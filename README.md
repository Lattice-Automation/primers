# primers

This is a small, straightforward tool for creating PCR primers. Its target use-case is DNA assembly.

Reasons to choose `primers` instead of [Primer3](https://github.com/primer3-org/primer3) include its:

- **features**: It is uniquely focused on DNA assembly flows like Gibson Assembly and Golden Gate cloning. You can design primers while adding sequence to the 5' ends of primers.
- **simplicity**: It is a small and simple Python CLI/library with a single dependency ([seqfold](https://github.com/Lattice-Automation/seqfold)). It is easier to install and use.
- **interface**: The Python library accepts and create primers for [Biopython `Seq`](https://biopython.org/wiki/Seq) classes. It outputs JSON for easy integration with other applications.
- **license**: It has a permissive, business-friendly license (MIT) instead of a copyleft GPL v2 license.

## Installation

```bash
pip install primers
```

## Usage

`primers` chooses pairs while optimizing for length, tm, GC ratio, secondary structure, and off-target binding. In the simplest case, you just pass the sequence you want to amplify:

```bash
$ primers create CTACTAATAGCACACACGGGGACTAGCATCTATCTCAGCTACGATCAGCATC
  dir    tm   ttm  gc     dg     p  seq
  FWD  63.6  63.6 0.5      0   2.6  CTACTAATAGCACACACGGG
  REV  63.2  63.2 0.5  -0.16  1.52  GATGCTGATCGTAGCTGAGATA
```

Additional sequence is added to the 5' end of primers via the `add_fwd/add_rev` args (`-f/-r` with CLI). By default, it will prepend the entire additional sequence. If you want it to choose the best subsequence to add to the 5' end (factoring in the features dicussed [below](#scoring)), allow it to choose from a range of indicies via the `add_fwd_len/add_rev_len` (`-fl/-rl` with CLI). Each primer has two tms: "tm", the melting temperature for the portion of the primer that binds to the template sequence and "tm_total", the melting temperature for the entire primer including the additional sequence added to primers' 5' end.

### Python

```python
from primers import create

# add enzyme recognition sequences to FWD and REV primers: BsaI, BpiI
fwd, rev = create("AATGAGACAATAGCACACACAGCTAGGTCAGCATACGAAA", add_fwd="GGTCTC", add_rev="GAAGAC")
print(fwd.fwd)      # True
print(fwd.seq)      # GGTCTCAATGAGACAATAGCACACACA; 5' to 3'
print(fwd.tm)       # 62.4; melting temp
print(fwd.tm_total) # 68.6; melting temp with added seq (GGTCTC)
print(fwd.dg)       # -1.86; minimum free energy of the secondary structure

# add from a range of sequence to the FWD primer: [5, 12] bp
fwd, rev = create("AATGAGACAATAGCACACACAGCTAGGTCAGCATACGAAA", add_fwd="GGATCGAGCTTGA", add_fwd_len=(5, 12))
print(fwd.seq)      # AGCTTGAAATGAGACAATAGCACACACAGC (AGCTTGA added from add_fwd)
print(fwd.tm)       # 62.2
print(fwd.tm_total) # 70.0
```

### CLI

```
$ primers create --help
usage: primers create [-h] [-f SEQ] [-fl INT INT] [-r SEQ] [-rl INT INT] [-t SEQ] [-j | --json | --no-json] SEQ

positional arguments:
  SEQ                   create primers to amplify this sequence

options:
  -h, --help            show this help message and exit
  -f SEQ                additional sequence to add to FWD primer (5' to 3')
  -fl INT INT           space separated min-max range for the length to add from '-f' (5' to 3')
  -r SEQ                additional sequence to add to REV primer (5' to 3')
  -rl INT INT           space separated min-max range for the length to add from '-r' (5' to 3')
  -t SEQ                sequence to check for off-target binding sites
  -j, --json, --no-json
                        write the primers to a JSON array
```

#### Table Output Format

By default, the primers are logged in table format in rows of `dir, tm, ttm, gc, dg, p, seq` where:

- dir: FWD or REV
- tm: the melting temperature of the annealing portion of the primer (Celsius)
- ttm: the total melting temperature of the primer with added seq (Celsius)
- gc: the GC ratio of the primer
- dg: the minimum free energy of the primer (kcal/mol)
- p: the primer's penalty score. Lower is better
- seq: the sequence of the primer in the 5' to the 3' direction

```txt
$ primers create -f GGTCTC -r GAAGAC AATGAGACAATAGCACACACAGCTAGGTCAGCATACGAAA
  dir    tm   ttm  gc     dg     p  seq
  FWD  60.8  67.0 0.5  -1.86  5.93  GGTCTCAATGAGACAATAGCACACAC
  REV  60.8  65.8 0.5      0   3.2  GAAGACTTTCGTATGCTGACCTAG
```

#### JSON Output Format

The `--json` flag prints primers in JSON format with more details on scoring. The example below is truncated for clarity:

```txt
$ primers create CTACTAATAGCACACACGGGGACTAGCATCTATCTCAGCTACGATCAGCATC --json| jq
[
  {
    "seq": "CTACTAATAGCACACACGGG",
    "len": 20,
    "tm": 63.6,
    "tm_total": 63.6,
    "gc": 0.5,
    "dg": 0,
    "fwd": true,
    "off_target_count": 0,
    "scoring": {
      "penalty": 2.6,
      "penalty_tm": 1.6,
      "penalty_tm_diff": 0,
      "penalty_gc": 0,
      "penalty_len": 1,
      "penalty_dg": 0,
      "penalty_off_target": 0
    }
  },
...
```

## Algorithm

Choosing PCR primers requires optimizing for a few different characteristics. Ideally, pairs of primers for PCR amplification would have similar tms, GC ratios close to 0.5, high minimum free energies (dg), and a lack off-target binding sites. In `primers`, like Primer3, choosing amongst those (sometimes competing) goals is accomplished with a linear function that penalizes undesirable characteristics. The primer pair with the lowest combined penalty is chosen.

### Scoring

The penalty for each possible primer, `p`, is calculated as:

```txt
PENALTY(p) =
    abs(p.tm - optimal_tm) * penalty_tm +     // penalize each deg of suboptimal melting temperature
    abs(p.gc - optimal_gc) * penalty_gc +     // penalize each percentage point of suboptimal GC ratio
    abs(len(p) - optimal_len) * penalty_len + // penalize each bp of suboptimal length
    abs(p.tm - p.pair.tm) * penalty_tm_diff + // penalize each deg of melting temperature diff between primers
    abs(p.dg) * penalty_dg +                  // penalize each kcal/mol of free energy in secondary structure
    p.offtarget_count * penalty_offtarget     // penalize each off-target binding site
```

Each of the optimal (`optimal_*`) and penalty (`penalty_*`) parameters is adjustable in the `primers.create()` function. The defaults are below:

```python
optimal_tm: float = 62.0
optimal_gc: float = 0.5
optimal_len: int = 22
penalty_tm: float = 1.0
penalty_gc: float = 0.2
penalty_len: float = 0.5
penalty_tm_diff: float = 1.0
penalty_dg: float = 2.0
penalty_offtarget: float = 20.0
```

#### Scoring Existing Primers

If you already have primers, and you want to see their features and penalty score, use the `primers score` command. The command below scores a FWD and REV primer against the sequence `-s` that they were created to amplify:

```txt
$ primers score GGTCTCAATGAGACAATA TTTCGTATGCTGACCTAG -s AATGAGACAATAGCACACACAGCTAGGTCAGCATACGAAATTT --json | jq
[
  {
    "seq": "GGTCTCAATGAGACAATA",
    "len": 18,
    "tm": 39.4,
    "tm_total": 55,
    "gc": 0.4,
    "dg": -1.86,
    "fwd": true,
    "off_target_count": 0,
    "scoring": {
      "penalty": 49.9,
      "penalty_tm": 22.6,
      "penalty_tm_diff": 19.6,
      "penalty_gc": 2,
      "penalty_len": 2,
      "penalty_dg": 3.7,
      "penalty_off_target": 0
    }
  },
  {
    "seq": "TTTCGTATGCTGACCTAG",
    "len": 18,
    "tm": 59,
    "tm_total": 59,
    "gc": 0.5,
    "dg": 0,
    "fwd": false,
    "off_target_count": 0,
    "scoring": {
      "penalty": 24.6,
      "penalty_tm": 3,
      "penalty_tm_diff": 19.6,
      "penalty_gc": 0,
      "penalty_len": 2,
      "penalty_dg": 0,
      "penalty_off_target": 0
    }
  }
]
```

### Off-target Binding Sites

Usually, off-target binding sites should be avoided. In `primers`, off-target binding sites are those with `<= 1` mismatch in the last 10 bair pairs of the primer's 3' end. This definition experimentally supported by:

> Wu, J. H., Hong, P. Y., & Liu, W. T. (2009). Quantitative effects of position and type of single mismatch on single base primer extension. Journal of microbiological methods, 77(3), 267-275

By default, primers are checked for off-targets within the `seq` parameter passed to `create(seq)`. But the primers can be checked against another sequence if it is passed to the optional `offtarget_check` argument (`-t` for CLI). This is useful when PCR'ing a subsequence of a larger DNA sequence like a plasmid.

```python
from primers import create

seq = "AATGAGACAATAGCACACACAGCTAGGTCAGCATACGAAA"
seq_parent = "ggaattacgtAATGAGACAATAGCACACACAGCTAGGTCAGCATACGAAAggaccagttacagga"

# primers are checked for offtargets in `seq_parent`
fwd, rev = create(seq, offtarget_check=seq_parent)
```
