"""Create or score PCR primers"""

import heapq
from logging import warning
from typing import Any, Dict, Tuple, NamedTuple, List, Optional

from seqfold import gc_cache, dg_cache, tm_cache
from .off_targets import off_targets


LEN_MIN = 15  # min length of the annealing portion of primers
LEN_MAX = 32  # max length of the annealing portion of primers, based on IDT guidelines


class Scoring(NamedTuple):
    """A scoring for a single Primer."""

    penalty: float
    """The high-level penalty for this primer"""

    penalty_tm: float
    """Penalty for each degree of tm suboptimality (diff from optimal)"""

    penalty_tm_diff: float
    """Penalty for each degree of tm difference between primers in a pair"""

    penalty_gc: float
    """Penalty for each percentage point of GC suboptimality (diff from optional)"""

    penalty_len: float
    """Penalty for each base pair length of suboptimality (diff from optimal)"""

    penalty_dg: float
    """Penalty for every kcal/mol of free energy"""

    penalty_off_target: float
    """Penalty for each off-target binding site"""


class Primer(NamedTuple):
    """A single Primer for PCR amplification of a DNA sequence."""

    seq: str
    """The DNA sequence of the primer; 5' to 3'"""

    len: int
    """The length of the seq"""

    tm: float
    """The melting temperature of the primer (Celsius) for the
    binding region pre-addition of added sequence"""

    tm_total: float
    """The melting temperature of the total primer (Celsius):
    the tm of the primer with the binding and added sequence"""

    gc: float
    """The GC ratio of the primer"""

    dg: float
    """The minimum free energy of the primer (kcal/mol)"""

    fwd: bool
    """Whether the primer anneals in the FWD direction of the template sequence"""

    off_target_count: int
    """The count of off-targets in the primer"""

    scoring: Scoring
    """Scoring of this primer (contains penalty)"""

    @property
    def penalty(self) -> float:
        """Penalty of the primer."""
        return self.scoring.penalty

    def dict(self) -> Dict[str, Any]:
        j = self._asdict()
        j["scoring"] = self.scoring._asdict()
        return j


class PrimerFactory(NamedTuple):
    """A factory for creating Primers with penalties.

    Holds the optimal values for a primer and the penalty for differences
    between primers' properties and optimal values.
    """

    optimal_tm: float
    """Optimal tm of a primer"""

    optimal_gc: float
    """Optimal GC ratio of a primer"""

    optimal_len: int
    """Optimal length of a primer"""

    penalty_tm: float
    """Penalty for each degree of tm suboptimality (diff from optimal)"""

    penalty_tm_diff: float
    """Penalty for each degree of tm difference between primers in a pair"""

    penalty_gc: float
    """Penalty for each percentage point of GC suboptimality (diff from optional)"""

    penalty_len: float
    """Penalty for each base pair length of suboptimality (diff from optimal)"""

    penalty_dg: float
    """Penalty for every kcal/mol of free energy"""

    penalty_off_target: float
    """Penalty for each off-target binding site"""

    def build(
        self,
        seq: str,
        tm: float,
        tm_total: float,
        gc: float,
        dg: float,
        fwd: bool,
        off_target_count: int,
    ) -> Primer:
        """Create a Primer with a scored penalty.

        Args:
            seq: Sequence of the primer, 5' to 3'
            tm: Tm of the created primer, Celsius
            tm_total: Tm of the created primer, with added sequence, Celsius
            gc: GC ratio of the created primer
            dg: Minimum free energy (kcal/mol) of the folded DNA sequence
            fwd: Whether this is a FWD primer
            off_target_count: The number of offtarget binding sites in the template sequence

        Returns:
            Primer: A Primer with a penalty score
        """

        dg = min(dg, 0)
        penalty_tm = abs(tm - self.optimal_tm) * self.penalty_tm
        penalty_gc = abs(gc - self.optimal_gc) * self.penalty_gc * 100
        penalty_len = abs(len(seq) - self.optimal_len) * self.penalty_len
        penalty_dg = abs(dg) * self.penalty_dg
        penalty_off_target = off_target_count * self.penalty_off_target
        penalty = (
            penalty_tm + penalty_gc + penalty_len + penalty_dg + penalty_off_target
        )

        return Primer(
            seq=seq,
            len=len(seq),
            tm=tm,
            tm_total=tm_total,
            gc=round(gc, 2),
            dg=round(dg, 2),
            fwd=fwd,
            off_target_count=off_target_count,
            scoring=Scoring(
                penalty_tm=round(penalty_tm, 1),
                penalty_tm_diff=0,  # unknown at this point
                penalty_len=penalty_len,
                penalty_gc=round(penalty_gc, 1),
                penalty_dg=round(penalty_dg, 1),
                penalty_off_target=penalty_off_target,
                penalty=round(penalty, 1),
            ),
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

        new_fwd = fwd._replace(
            scoring=fwd.scoring._replace(
                penalty=round(fwd.scoring.penalty + penalty_tm_diff, 1),
                penalty_tm_diff=round(penalty_tm_diff, 1),
            )
        )
        new_rev = rev._replace(
            scoring=rev.scoring._replace(
                penalty=round(rev.scoring.penalty + penalty_tm_diff, 1),
                penalty_tm_diff=round(penalty_tm_diff, 1),
            )
        )

        return new_fwd, new_rev


def score(
    fwd: str,
    rev: str = "",
    seq: str = "",
    offtarget_check: str = "",
    optimal_tm: float = 62.0,
    optimal_gc: float = 0.5,
    optimal_len: int = 22,
    penalty_tm: float = 1.0,
    penalty_gc: float = 0.2,
    penalty_len: float = 0.5,
    penalty_tm_diff: float = 1.0,
    penalty_dg: float = 2.0,
    penalty_off_target: float = 20.0,
) -> Tuple[Primer, Optional[Primer]]:
    """Score primers from their sequence.

    Use-case: you already have primers and want to know their characteristics and scoring.

    Args:
        fwd: the sequence of the first primer

    Keyword Args:
        rev: optional sequence of the second primer
        seq: the sequence that the primers anneal to/amplify
        add_fwd: additional sequence to add to FWD primer (5' to 3')
        add_rev: additional sequence to add to REV primer (5' to 3')
        add_fwd_len: range (min, max) of number of bp to add from
            `add_fwd` (from 3' end)
        add_rev_len: range (min, max) of number of bp to add from
            `add_rev` (from 3' end)
        offtarget_check: the sequence to check for offtarget binding sites
        optimal_tm: the optimal tm of a primer based on IDT guidelines.
            Excluding added sequence
        optimal_gc: the optimal GC ratio of a primer, based on IDT guidelines
        optimal_len: the optimal length of a primer, excluding additional
            sequence added via `add_fwd` and `add_rev`
        penalty_tm: penalty for tm differences from optimal
        penalty_gc: penalty for each percentage point diff between primers
            and the optimal GC ratio
        penalty_len: penalty for differences in primer length
        penalty_diff_tm: penalty for tm differences between primers
        penalty_dg: penalty for minimum free energy of a primer
        penalty_off_target: penalty for offtarget binding sites in the `seq`
    """

    if len(fwd) < LEN_MIN:
        raise ValueError("`fwd` primer is too short: {} < {}", len(fwd), LEN_MIN)
    if rev and len(rev) < LEN_MIN:
        raise ValueError("`rev` primer is too short: {} < {}", len(rev), LEN_MIN)

    fwd = fwd.upper()
    rev = rev.upper()
    seq, off_target_check = _parse(seq, offtarget_check)
    add_fwd, _, add_rev = _binding_seq(fwd, rev=rev, seq=seq)

    factory = PrimerFactory(
        optimal_tm=optimal_tm,
        optimal_gc=optimal_gc,
        optimal_len=optimal_len,
        penalty_tm=penalty_tm,
        penalty_gc=penalty_gc,
        penalty_len=penalty_len,
        penalty_tm_diff=penalty_tm_diff,
        penalty_dg=penalty_dg,
        penalty_off_target=penalty_off_target,
    )

    fwd_primer = factory.build(
        fwd,
        tm=tm_cache(fwd)[add_fwd][-1],
        tm_total=tm_cache(fwd)[0][-1],
        dg=dg_cache(fwd)[0][-1],
        gc=gc_cache(fwd)[0][-1],
        fwd=True,
        off_target_count=off_targets(fwd, off_target_check)[len(fwd) - 1],
    )

    if not rev:
        return (fwd_primer, None)

    rev_primer = factory.build(
        rev,
        tm=tm_cache(rev)[add_rev][-1],
        tm_total=tm_cache(rev)[0][-1],
        dg=dg_cache(rev)[0][-1],
        gc=gc_cache(rev)[0][-1],
        fwd=False,
        off_target_count=off_targets(rev, off_target_check)[len(rev) - 1],
    )

    return factory.build_pair(fwd_primer, rev_primer)


def _binding_seq(fwd: str, rev: str = "", seq: str = "") -> Tuple[int, str, int]:
    """Attempt to find the binding region between the fwd and rev primers in seq"""

    if not seq:
        return (0, "", 0)

    fwd = fwd.upper()
    rev = rev.upper()
    seq = seq.upper()
    seq += seq + seq  # account for amplifications across the zero-index

    add_fwd = 0
    add_rev = 0

    # remove start of seq prior to binding site
    try:
        start_index_fwd = -10  # index on primer
        start_index_seq = seq.index(fwd[start_index_fwd:])

        while (
            start_index_seq
            and len(fwd) + start_index_fwd > 0
            and seq[start_index_seq] == fwd[start_index_fwd]
        ):
            start_index_fwd -= 1
            start_index_seq -= 1

        add_fwd = len(fwd) + start_index_fwd
        seq = seq[start_index_seq:]
    except Exception as err:
        raise ValueError("failed to find `fwd` binding site in `seq`: %s", err)

    if not rev:
        return (add_fwd, seq, add_rev)

    # remove end of seq after the reverse primer binding site
    try:
        end_index_rev = 10
        rev = _rc(rev)
        end_index_seq = seq.index(rev[:end_index_rev])

        while (
            end_index_rev < len(rev) - 1
            and rev[: end_index_rev + 2] in seq[end_index_seq:]
        ):
            end_index_rev += 1

        add_rev = len(rev) - end_index_rev - 1
        seq = seq[: end_index_seq + end_index_rev + 1]
    except Exception as err:
        raise ValueError("failed to find `rev` binding site in `seq`: %s", err)

    return add_fwd, seq, add_rev


def primers(
    seq: str,
    add_fwd: str = "",
    add_rev: str = "",
    add_fwd_len: Tuple[int, int] = (-1, -1),
    add_rev_len: Tuple[int, int] = (-1, -1),
    offtarget_check: str = "",
    optimal_tm: float = 62.0,
    optimal_gc: float = 0.5,
    optimal_len: int = 22,
    penalty_tm: float = 1.0,
    penalty_gc: float = 0.2,
    penalty_len: float = 0.5,
    penalty_tm_diff: float = 1.0,
    penalty_dg: float = 2.0,
    penalty_off_target: float = 20.0,
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
        optimal_tm: The optimal tm of a primer based on IDT guidelines.
            Excluding added sequence
        optimal_gc: The optimal GC ratio of a primer, based on IDT guidelines
        optimal_len: The optimal length of a primer, excluding additional
            sequence added via `add_fwd` and `add_rev`
        penalty_tm: Penalty for tm differences from optimal
        penalty_gc: Penalty for each percentage point diff between primers
            and the optimal GC ratio
        penalty_len: Penalty for differences in primer length
        penalty_diff_tm: Penalty for tm differences between primers
        penalty_dg: Penalty for minimum free energy of a primer
        penalty_off_target: Penalty for offtarget binding sites in the `seq`

    Returns:
        (Primer, Primer): Primers for PCR amplification
    """

    # parse input
    seq, offtarget_check = _parse(seq, offtarget_check)
    add_fwd, _ = _parse(add_fwd, "")
    add_rev, _ = _parse(add_rev, "")
    factory = PrimerFactory(
        optimal_tm=optimal_tm,
        optimal_gc=optimal_gc,
        optimal_len=optimal_len,
        penalty_tm=penalty_tm,
        penalty_gc=penalty_gc,
        penalty_len=penalty_len,
        penalty_tm_diff=penalty_tm_diff,
        penalty_dg=penalty_dg,
        penalty_off_target=penalty_off_target,
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
    optimal_fwd_len = round(optimal_len + (add_fwd_min + add_fwd_max) / 2)
    fwd_seq = seq_full[: add_fwd_max + LEN_MAX]
    fwd_primers = _primers(
        factory._replace(optimal_len=optimal_fwd_len),
        fwd_seq,
        offtarget_check,
        range(0, add_fwd_max - add_fwd_min + 1),
        range(LEN_MIN + add_fwd_max, LEN_MAX + add_fwd_max),
        True,
        add_fwd_max,
    )

    optimal_rev_len = round(optimal_len + (add_rev_min + add_rev_max) / 2)
    rev_seq = _rc(seq_full)
    rev_seq = rev_seq[: add_rev_max + LEN_MAX]
    rev_primers = _primers(
        factory._replace(optimal_len=optimal_rev_len),
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
    ot = off_targets(seq, offtarget_check)
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
            heapq.heappush(ranked_fwd, (p.scoring.penalty, p))

    for row in rev_primers:
        for p in row:
            if not p:
                continue
            heapq.heappush(ranked_rev, (p.scoring.penalty, p))

    if not ranked_fwd:
        raise RuntimeError("Failed to create any primers in the FWD direction")

    if not ranked_rev:
        raise RuntimeError("Failed to create any primers in the REV direction")

    min_penalty = float("inf")
    min_fwd, min_rev = None, None
    for _, fwd in heapq.nsmallest(10, ranked_fwd):
        for _, rev in heapq.nsmallest(10, ranked_rev):
            new_fwd, new_rev = factory.build_pair(fwd, rev)
            new_penalty = new_fwd.scoring.penalty + new_rev.scoring.penalty
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
    offtarget_check = (offtarget_check or seq).upper()

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
