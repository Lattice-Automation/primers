"""Create primers to PCR amplify a DNA sequence.

Primers are used to PCR amplify DNA sequences. This `primers` python module
is made with synthetic biology and DNA assembly in mind, specifically.
It makes it easy to PCR amplify a DNA sequence while adding new sequences to its ends for:

1. adding homology for neighboring fragments (eg: Gibson Assembly)
2. adding restriction enzyme sites (eg: Golden Gate cloning).

Selecting primers for a DNA sequence is non-trivial because it's
a multi-objective optimization problem. Ideally, pairs of primers
for PCR amplification would have similar, ideal tms, low gc%s, low
free energies (dgs) and lack off-target binding sites.

In this module, the penalty for each possible primer, p, is calculated as:
    PENALTY(p) =
        abs(p.tm - opt_tm) * penalty_tm +
        abs(p.gc - opt_gc) * penalty_gc +
        abs(len(p) - opt_len) * penalty_len +
        abs(p.tm - p.pair.tm) * penalty_tm_diff +
        abs(p.dg) * penalty_dg +
        p.offtargets * penalty_offtarget

The primer pair with the lowest combined penalty score is chosen.

Given this module's emphasis on DNA assembly, additional sequences added to the FWD and/or REV primer
are considered in the PENALTY calculation.

An offtarget binding site of a primer is defined as a subsequence
within one mismatch of the last 10bp of the primer's 3' end. This is experimentally
supported by:

    Wu, J. H., Hong, P. Y., & Liu, W. T. (2009). Quantitative effects of
    position and type of single mismatch on single base primer extension.
    Journal of microbiological methods, 77(3), 267-275.
"""

import heapq
from logging import warning
from typing import Tuple, NamedTuple, List, Optional

from seqfold import gc_cache, dg_cache, tm_cache, Cache
from .offtargets import offtargets


LEN_MIN = 15  # min length of the annealing portion of primers
LEN_MAX = 32  # max length of the annealing portion of primers, based on IDT guidelines

PRIMER_FMT: str = "{:>5} {:>5} {:>5} {:>6} {:>5}  {}"
"""{fwd} {tm} {tm_total} {dg} {penalty} {seq}"""


class Primer(NamedTuple):
    """A single Primer for PCR amplification of a DNA sequence.

    Attributes:
        seq: The DNA sequence of the primer; 5' to 3'
        tm: The melting temperature of the primer (Celsius):
            for the binding region pre-addition of added sequence
        tm_total: The melting temperature of the total primer (Celsius):
            the tm of the primer with the binding and added sequence
        gc: The GC percentage of the primer
        dg: The minimum free energy of the primer
        fwd: Whether the primer anneals in the FWD
            direction of the template sequence
        penalty: The penalty score for this primer
    """

    seq: str
    tm: float
    tm_total: float
    gc: float
    dg: float
    fwd: bool
    offtargets: int
    penalty: float

    def __str__(self) -> str:
        """Create a string representation of the primer."""

        return PRIMER_FMT.format(
            "FWD" if self.fwd else "REV",
            self.tm,
            self.tm_total,
            round(self.dg, 2),
            round(self.penalty, 2),
            self.seq,
        )


class PrimerFactory(NamedTuple):
    """A factory for creating Primers with penalties.

    Holds the optimal values for a primer and the penalty for differences
    between primers' properties and those optimal values.
    
    Attributes:
        opt_tm: Optimal tm of a primer
        opt_gc: Optimal GC ratio of a primer
        opt_len: Optimal length of a primer
        penalty_tm: Penalty for a large tm difference
        penalty_tm_diff: Penalty for differences between primers in a pair
        penalty_dg: Penalty for very negative free energies
        penalty_offtarget: Penalty for offtargets
    """

    opt_tm: float
    opt_gc: float
    opt_len: int
    penalty_tm: float
    penalty_gc: float
    penalty_len: float
    penalty_tm_diff: float
    penalty_dg: float
    penalty_offtarget: float

    def build(
        self,
        seq: str,
        tm: float,
        tm_total: float,
        gc: float,
        dg: float,
        fwd: bool,
        offtargets: int,
    ) -> Primer:
        """Create a Primer with a scored penalty.
        
        Args:
            seq: Sequence of the primer, 5' to 3'
            tm: Tm of the created primer, Celsius
            tm_total: Tm of the created primer, with added sequence, Celsius
            gc: GC ratio of the created primer
            dg: Minimum free energy (kcal/mol) of the folded DNA sequence
            fwd: Whether this is a FWD primer
            offtargets: The number of offtarget binding sites in the template sequence
        
        Returns:
            Primer: A Primer with a penalty score
        """

        dg = min(dg, 0)
        penalty_tm = abs(tm - self.opt_tm) * self.penalty_tm
        penalty_gc = abs(gc - self.opt_gc) * self.penalty_gc
        penalty_len = abs(len(seq) - self.opt_len) * self.penalty_len
        penalty_dg = abs(dg) * self.penalty_dg
        penalty_offtarget = offtargets * self.penalty_offtarget
        penalty = penalty_tm + penalty_gc + penalty_len + penalty_dg + penalty_offtarget

        return Primer(
            seq=seq,
            tm=tm,
            tm_total=tm_total,
            gc=gc,
            dg=dg,
            fwd=fwd,
            offtargets=offtargets,
            penalty=penalty,
        )

    def build_pair(self, fwd: Primer, rev: Primer) -> Tuple[Primer, Primer]:
        """Create a pair of Primers with a tm_diff penalty added.
        
        Args:
            fwd: The FWD primer
            rev: The REV primer
        
        Returns:
            (Primer, Primer): A Primer pair with a tm_diff penalty applied
        """

        penalty_tm_diff = abs(fwd.tm - rev.tm) * self.penalty_tm_diff

        new_fwd = fwd._replace(penalty=fwd.penalty + penalty_tm_diff)
        new_rev = rev._replace(penalty=rev.penalty + penalty_tm_diff)

        return new_fwd, new_rev


def primers(
    seq: str,
    add_fwd: str = "",
    add_rev: str = "",
    add_fwd_len: Tuple[int, int] = (-1, -1),
    add_rev_len: Tuple[int, int] = (-1, -1),
    offtarget_check: str = "",
    opt_tm: float = 62.0,
    opt_gc: float = 0.5,
    opt_len: int = 22,
    penalty_tm: float = 1.0,
    penalty_gc: float = 3.0,
    penalty_len: float = 1.0,
    penalty_tm_diff: float = 1.0,
    penalty_dg: float = 2.0,
    penalty_offtarget: float = 20.0,
) -> Tuple[Primer, Primer]:
    """Create primers for PCR amplification of the sequence.
    
    Args:
        seq: The DNA sequence to amplify
    
    Keyword Args:
        add_fwd: Additional sequence to add to FWD primer (5' to 3')
        add_rev: Additional sequence to add to REV primer (5' to 3')
        add_fwd_len: Range (min, max) of number of bp to add from
            `add_fwd` (from 3' end)
        add_rev_len: Range (min, max) of number of bp to add from
            `add_rev` (from 3' end)
        offtarget_check: The sequence to check for offtarget binding sites
        opt_tm: The optimal tm of a primer based on IDT guidelines.
            Excluding added sequence
        opt_gc: The optimal GC ratio of a primer, based on IDT guidelines
        opt_len: The optimal length of a primer, excluding additional
            sequence added via `add_fwd` and `add_rev`
        penalty_tm: Penalty for tm differences from optimal
        penalty_gc: Penalty for differences between primers and optimal GC ratio
        penalty_len: Penalty for differences in primer length
        penalty_diff_tm: Penalty for tm differences between primers
        penalty_dg: Penalty for minimum free energy of a primer
        penalty_offtarget: Penalty for offtarget binding sites in the `seq`
    
    Returns:
        (Primer, Primer): Primers for PCR amplification
    """

    # parse input
    seq, offtarget_check = _parse(seq, offtarget_check)
    add_fwd, _ = _parse(add_fwd, "")
    add_rev, _ = _parse(add_rev, "")
    factory = PrimerFactory(
        opt_tm=opt_tm,
        opt_gc=opt_gc,
        opt_len=opt_len,
        penalty_tm=penalty_tm,
        penalty_gc=penalty_gc,
        penalty_len=penalty_len,
        penalty_tm_diff=penalty_tm_diff,
        penalty_dg=penalty_dg,
        penalty_offtarget=penalty_offtarget,
    )

    # set min/max if additional sequence was provided at FWD/REV
    add_fwd_min, add_fwd_max = _parse_add_len(add_fwd, add_fwd_len)
    add_rev_min, add_rev_max = _parse_add_len(add_rev, add_rev_len)

    # create the template sequence
    add_fwd = add_fwd[-add_fwd_max:]
    add_rev = _rc(add_rev)[:add_rev_max]
    seq_full = add_fwd + seq + add_rev

    if len(seq_full) < LEN_MAX:
        err = f"Template sequence length is too short: {len(seq_full)}bp < {LEN_MAX}bp"
        raise ValueError(err)

    # create two 2D arrays of primers in the FWD and REV directions
    opt_fwd_len = round(opt_len + (add_fwd_min + add_fwd_max) / 2)
    fwd_seq = seq_full[: add_fwd_max + LEN_MAX]
    fwd_primers = _primers(
        factory._replace(opt_len=opt_fwd_len),
        fwd_seq,
        offtarget_check,
        range(0, add_fwd_max - add_fwd_min + 1),
        range(LEN_MIN + add_fwd_max, LEN_MAX + add_fwd_max),
        True,
        add_fwd_max,
    )

    opt_rev_len = round(opt_len + (add_rev_min + add_rev_max) / 2)
    rev_seq = _rc(seq_full)
    rev_seq = rev_seq[: add_rev_max + LEN_MAX]
    rev_primers = _primers(
        factory._replace(opt_len=opt_rev_len),
        rev_seq,
        offtarget_check,
        range(0, add_rev_max - add_rev_min + 1),
        range(LEN_MIN + add_rev_max, LEN_MAX + add_rev_min),
        False,
        add_rev_max,
    )

    # choose the pair with the smallest total penalty
    return _choose(factory, fwd_primers, rev_primers)


def _primers(
    factory: PrimerFactory,
    seq: str,
    offtarget_check: str,
    start_range: range,
    end_range: range,
    fwd: bool,
    add_len: int,
) -> List[List[Optional[Primer]]]:
    """Return a matrix of primers for (i, j) where (i, j) are the start/end indexes
    
    Args:
        factory: Primer factory for creating the primers, assigning score
        seq: The sequence who primers are being created for
        offtarget_check: the sequence checked for offtargets
        start_range: A range of possible starting indicies
        end_range: A range of possible ending indicies
        fwd: Whether these primers are in the FWD direction
        add_len: The number of additional bp added to the sequence's end
    
    Returns:
        List[List[Primer]]: A 2D matrix of Primers with penalties. None if no
            primer created within a range. None if it's an invalid range (eg j < i) or
            not one allowed for by creation settings
    """

    gc = gc_cache(seq)
    tm = tm_cache(seq)
    dg = dg_cache(seq)
    ot = offtargets(seq, offtarget_check)

    assert len(gc) == len(tm) == len(dg)

    ps: List[List[Optional[Primer]]] = []
    for _ in range(len(gc)):
        ps.append([None] * len(gc))

    for s in start_range:
        for e in end_range:
            p_seq = seq[s : e + 1]
            p_tm = tm[s + add_len][e]
            p_tm_total = tm[s][e]
            p_gc = gc[s][e]
            p_dg = dg[s][e]
            p_ot = ot[e]

            ps[s][e] = factory.build(p_seq, p_tm, p_tm_total, p_gc, p_dg, fwd, p_ot)

    return ps


def _choose(
    factory: PrimerFactory,
    fwd_primers: List[List[Optional[Primer]]],
    rev_primers: List[List[Optional[Primer]]],
) -> Tuple[Primer, Primer]:
    """Choose the best combo of primers. One in FWD direction and one in the REV direction.
    
    The 10 best primers are chosen in the FWD and REV direction and their
    tm is compared to create NEW primers with the tm_diff penalty applied

    Args:
        factory: PrimerFactory for creating new primers with tm_diff penalties
        fwd_primers: FWD primer 2D array from `_primers`
        rev_primers: REV primer 2D array from `_primers`
    
    Returns:
        (Primer, Primer): The 'best' combo of primers for PCR
    """

    ranked_fwd: List[Tuple[float, Primer]] = []
    ranked_rev: List[Tuple[float, Primer]] = []

    for row in fwd_primers:
        for p in row:
            if not p:
                continue
            heapq.heappush(ranked_fwd, (p.penalty, p))

    for row in rev_primers:
        for p in row:
            if not p:
                continue
            heapq.heappush(ranked_rev, (p.penalty, p))

    if not ranked_fwd:
        raise RuntimeError("Failed to create any primers in the FWD direction")

    if not ranked_rev:
        raise RuntimeError("Failed to create any primers in the REV direction")

    min_penalty = float("inf")
    min_fwd, min_rev = None, None
    for _, fwd in heapq.nsmallest(10, ranked_fwd):
        for _, rev in heapq.nsmallest(10, ranked_rev):
            new_fwd, new_rev = factory.build_pair(fwd, rev)
            new_penalty = new_fwd.penalty + new_rev.penalty
            if new_penalty < min_penalty:
                min_penalty = new_penalty
                min_fwd, min_rev = fwd, rev

    if not min_fwd or not min_rev:
        raise RuntimeError("Failed to create a pair of PCR primers")

    return min_fwd, min_rev


def _parse(seq: str, offtarget_check: str) -> Tuple[str, str]:
    """Validate and parse the input sequence.

    Ensure it's just DNA, uppercase it, return. If it's a BioRecord, get sequence.

    Args:
        seq: Template DNA sequence
        offtarget_check: The sequence that's checked for offtargets
    
    Raises:
        ValueError: If invalid bases are in the DNA, anything other than {ATGC}
    
    Returns:
        str: The parsed DNA sequence to be amplified
    """

    if "SeqRecord" in str(type(seq)):
        seq = str(seq.seq)  # type: ignore

    seq = seq.upper()
    offtarget_check = offtarget_check or seq

    if "SeqRecord" in str(type(offtarget_check)):
        offtarget_check = str(offtarget_check.seq).upper()  # type: ignore

    diff = {c for c in seq if c not in "ATGC"}

    if diff:
        desc = ",".join(diff)
        raise ValueError(f"Invalid non-DNA bases found: {desc}")

    return seq, offtarget_check


def _parse_add_len(add: str, add_len: Tuple[int, int]) -> Tuple[int, int]:
    """Parse the additional length range

    Args:
        add: The additional sequence being added
        add_len: A tuple with a min and max number of bp to add, inclusive
    
    Returns:
        (int, int): A tuple with a min and max number of bp to add
    """

    if not add:
        # we're not adding anything
        return 0, 0

    if add_len == (-1, -1):
        # just add the whole sequence if <40bp
        # if greater than 40 bp, log a warning and add (30, 40) bp
        add_assume = min(40, len(add))
        if add_assume != len(add):
            warning(
                f"{len(add)}bp additional sequence added, but no `add_DIR_len` argument provided:"
                + f"\n\tAdding between {add_assume - 10} and {add_assume}bp"
            )
            return add_assume - 10, add_assume
        return add_assume, add_assume

    # one or both of them were set
    add_min, add_max = add_len
    if add_min < 0 and add_max >= 0:
        add_max = min(len(add), add_max)
        return 0, add_max
    if add_min >= 0 and add_max < 0:
        add_min = min(len(add), add_min)
        return add_min, len(add)

    if add_min > -1 and add_max > -1 and add_min > add_max:
        raise ValueError(f"add_DIR_len range has a min > max: {add_min} > {add_max}")

    add_min = max(add_min, 0)
    add_max = min(add_max, len(add))
    return add_min, add_max


def _rc(seq: str) -> str:
    """Return the reverse complement of a DNA sequence.

    Args:
        seq: The template sequence
    
    Returns:
        str: The reverse complement of the template sequence
    """

    rc = {"A": "T", "T": "A", "G": "C", "C": "G"}
    return "".join(rc[c] for c in reversed(seq))
