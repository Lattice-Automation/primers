"""Find off-target binding sites.
"""

from collections import defaultdict
from typing import List, Dict


def off_targets(seq: str, check_seq: str) -> List[int]:
    """Return a list of off-target counts for primers that end is that index.

    For example, offtarget_cache[20] -> returns the number of offtarget binding
    sites whose last bp ends in the 20th index of `seq`

    Args:
        seq: The template sequence being checked
        check_seq: The sequence being checked for offtarget binding sites

    Returns:
        List[int]: A list of binding site counts for primers whose last
            bp is within that index of the list
    """

    check_map: Dict[str, int] = defaultdict(int)
    mutate_map = {"A": "TGC", "T": "AGC", "G": "ATC", "C": "ATG"}

    def mutate(s: str) -> List[str]:
        mutated_sites = [s]
        for i, c in enumerate(s):
            for m in mutate_map[c]:
                mutated_sites.append(s[:i] + m + s[i + 1 :])
        return mutated_sites

    for s in range(len(check_seq) - 9):
        for m in mutate(check_seq[s : s + 10]):
            check_map[m] += 1

    cache: List[int] = [0] * len(seq)
    for e in range(10, len(seq) + 1):
        ss = seq[e - 10 : e]

        # assumed to bind at least once
        cache[e - 1] = max(0, check_map[ss] + check_map[_rc(ss)] - 1)
    return cache


def _rc(seq: str) -> str:
    """Return the reverse complement of a DNA sequence.

    Args:
        seq: The template sequence

    Returns:
        str: The reverse complement of the template sequence
    """

    rc = {"A": "T", "T": "A", "G": "C", "C": "G"}
    return "".join(rc[c] for c in reversed(seq))
