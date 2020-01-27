# primers

This is a tool for making PCR primers for DNA sequences. `primers`' emphasis is on ease-of-use and DNA assembly. Adding additional sequence to 5' end of the FWD and REV primer, something common in creating fragments for Gibson assembly and Golden Gate cloning, is easy with `primers`. Finally, it has a permissive MIT license where other primer design tools don't.

## Installation

```bash
pip install primers
```

## Usage

### Python

```python
from primers import primers

# add recognition sequences to FWD and REV primers
fwd, rev = primers("AATGAGACAATAGCACACACAGCTAGGTCAGCATACGAAA", add_fwd="GGTCTC" add_rev="GAAGAC")
print(fwd.fwd)  # True
print(fwd.seq)  # GGTCTCAATGAGACAATAGCACACACA; 5' to 3'
print(fwd.tm)   # 62.4; melting temp
print(fwd.tm_total)  # 68.6; melting temp with added seq (GGTCTC)
print(fwd.dg)   # -1.86; minimum free energy of the secondary structure

# add from a range of sequence to the FWD primer
fwd, rev = primers("AATGAGACAATAGCACACACAGCTAGGTCAGCATACGAAA", add_fwd="GGATCGAGCTTGA", add_fwd_len=(5, 12))
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
usage: primers [-h] [-f SEQ] [-fl INT INT] [-r SEQ] [-rl INT INT] [--version] SEQ

Create PCR primers for a DNA sequence.

Logs the FWD and REV primer with columns:
dir, tm, ttm, dg, pen, seq

Where:
dir = FWD or REV.
tm  = Melting temperature of the annealing/binding part of the primer (Celsius).
ttm = The total melting temperature of the primer with added seq (Celsius).
dg  = The minimum free energy of the primer (kcal/mol).
pen = The primer's penalty score. Lower is better.
seq = The sequence of the primer in the 5' to the 3' direction.

positional arguments:
  SEQ                   DNA sequence

optional arguments:
  -h, --help            show this help message and exit
  -f SEQ, --fwd SEQ     additional sequence to add to FWD primer (5' to 3')
  -fl INT INT, --flen INT INT
                        space separated min-max range for the length to add from 'add_fwd' (5' to 3')
  -r SEQ, --rev SEQ     additional sequence to add to REV primer (5' to 3')
  -rl INT INT, --rlen INT INT
                        space separated min-max range for the length to add from 'add_rev' (5' to 3')
  --version             show program's version number and exit
```

## Overview

Selecting primers for a DNA sequence is non-trivial because it's a multi-objective optimization problem. Ideally, pairs of primers for PCR amplification would have similar, ideal tms, low gc%s, low free energies (dgs) and lack off-target binding sites.

### Scoring

In this module, the penalty for each possible primer, p, is calculated as:

```txt
PENALTY(p) =
    abs(p.tm - opt_tm) * penalty_tm +
    abs(p.gc - opt_gc) * penalty_gc +
    abs(len(p) - opt_len) * penalty_len +
    abs(p.tm - p.pair.tm) * penalty_tm_diff +
    abs(p.dg) * penalty_dg +
    p.offtargets * penalty_offtargets
```

Each of the optimal (`opt_*`) and penalty (`penalty_*`) parameters is adjustable through the `primers.primers()` function. The primer pair with the lowest combined penalty score is chosen.

Given this module's emphasis on DNA assembly, additional sequences added to the FWD and/or REV primer are considered in the PENALTY calculation.
