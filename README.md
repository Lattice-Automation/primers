# primers

This is a tool for creating PCR primers. It has an emphasis on DNA assembly and makes it easy to add sequences to the end of PCR fragments. This is part of [overlap extension polymerase chain reaction](https://en.wikipedia.org/wiki/Overlap_extension_polymerase_chain_reaction) and preparing unstandardized DNA sequences for Gibson assembly and Golden gate cloning.

`primers` quickly creates pairs with optimized lengths, Tms, GC ratios, secondary structures (minimum free energies) and without off-target binding sites. Each returned primer has two tms: "tm", the melting temperature for the portion of the primer that binds to the template sequence, and "tm_total", the melting temperature for the entire primer with additional sequence added to its 5' end.

Unlike the most used alternative, Primer3, `primers` has a permissive MIT license and support for adding sequence to the 5' ends of primers.

## Installation

```bash
pip install primers
```

## Usage

### Python

```python
from primers import primers

# add enzyme recognition sequences to FWD and REV primers: BsaI, BpiI
fwd, rev = primers("AATGAGACAATAGCACACACAGCTAGGTCAGCATACGAAA", add_fwd="GGTCTC" add_rev="GAAGAC")
print(fwd.fwd)  # True
print(fwd.seq)  # GGTCTCAATGAGACAATAGCACACACA; 5' to 3'
print(fwd.tm)   # 62.4; melting temp
print(fwd.tm_total)  # 68.6; melting temp with added seq (GGTCTC)
print(fwd.dg)   # -1.86; minimum free energy of the secondary structure

# add from a range of sequence to the FWD primer: [5, 12] bp
add_fwd = "GGATCGAGCTTGA"
fwd, rev = primers("AATGAGACAATAGCACACACAGCTAGGTCAGCATACGAAA", add_fwd=add_fwd, add_fwd_len=(5, 12))
print(fwd.seq)  # AGCTTGAAATGAGACAATAGCACACACAGC
print(fwd.tm)   # 62.2
print(fwd.tm_total)  # 70.0
```

### CLI

```txt
$ primers AATGAGACAATAGCACACACAGCTAGGTCAGCATACGAAA -f GGTCTC -r GAAGAC
  dir    tm   ttm     dg   pen  seq
  FWD  62.4  68.6  -1.86  5.43  GGTCTCAATGAGACAATAGCACACACA
  REV  62.8  67.4      0   4.8  GAAGACTTTCGTATGCTGACCTAG
```

```txt
$ primers --help
usage: primers [-h] [-f SEQ] [-fl INT INT] [-r SEQ] [-rl INT INT] [-t SEQ] [--version] SEQ

Create PCR primers for a DNA sequence.

Logs the FWD and REV primer with columns:
    dir, tm, ttm, dg, pen, seq

Where:
    dir = FWD or REV.
    tm  = Melting temperature of the annealing/binding part of the primer (Celsius).
    ttm = The total melting temperature of the primer with added seq (Celsius).
    dg  = The minimum free energy of the primer's secondary structure (kcal/mol).
    pen = The primer's penalty score. Lower is better.
    seq = The sequence of the primer in the 5' to the 3' direction.

positional arguments:
  SEQ          DNA sequence

optional arguments:
  -h, --help   show this help message and exit
  -f SEQ       additional sequence to add to FWD primer (5' to 3')
  -fl INT INT  space separated min-max range for the length to add from '-f' (5' to 3')
  -r SEQ       additional sequence to add to REV primer (5' to 3')
  -rl INT INT  space separated min-max range for the length to add from '-r' (5' to 3')
  -t SEQ       sequence to check for offtargets binding sites
  --version    show program's version number and exit
```

## Algorithm

Creating primers for a DNA sequence is non-trivial because it's multi-objective optimization. Ideally pairs of primers for PCR amplification would have similar tms, GC ratios close to 0.5, high minimum free energies (dg), and a lack off-target binding sites. In `primers`, like Primer3, this is accomplished with a linear function that penalizes undesired characteristics. The primer pair with the lowest combined penalty is created.

### Scoring

The penalty for each possible primer, p, is calculated as:

```txt
PENALTY(p) =
    abs(p.tm - opt_tm) * penalty_tm +
    abs(p.gc - opt_gc) * penalty_gc +
    abs(len(p) - opt_len) * penalty_len +
    abs(p.tm - p.pair.tm) * penalty_tm_diff +
    abs(p.dg) * penalty_dg +
    p.offtargets * penalty_offtarget
```

Each of the optimal (`opt_*`) and penalty (`penalty_*`) parameters is adjustable through the `primers.primers()` function. The defaults are below.

```python
opt_tm: float = 62.0
opt_gc: float = 0.5
opt_len: int = 22
penalty_tm: float = 1.0
penalty_gc: float = 3.0
penalty_len: float = 1.0
penalty_tm_diff: float = 1.0
penalty_dg: float = 2.0
penalty_offtarget: float = 20.0
```

### Off-targets

Off-targets are defined as a subsequence within one mismatch of the last 10bp of a primer's 3' end. This is experimentally supported by:

> Wu, J. H., Hong, P. Y., & Liu, W. T. (2009). Quantitative effects of position and type of single mismatch on single base primer extension. Journal of microbiological methods, 77(3), 267-275

By default, primers are checked for off-targets within the `seq` parameter passed to `primers.primers(seq)`. But the primers can be checked against another sequence if it's passed through the `offtarget_check` argument. This is useful when PCR'ing a subsequence of a larger DNA sequence; for example: a plasmid.

```python
seq = "AATGAGACAATAGCACACACAGCTAGGTCAGCATACGAAA"
parent = "ggaattacgtAATGAGACAATAGCACACACAGCTAGGTCAGCATACGAAAggaccagttacagga"

# primers are checked for offtargets in `parent`
fwd, rev = primers(seq, offtarget_check=parent)
```
