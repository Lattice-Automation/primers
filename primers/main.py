"""Console entrypoint for creating PCR primers"""

import argparse
import sys
from typing import List

from . import __version__, primers, Primer
from .primers import PRIMER_FMT


def run():
    """Entry point for console_scripts.
    
    Create primers and log the results.
    """

    args = parse_args(sys.argv[1:])

    fwd, rev = primers(
        args.seq,
        add_fwd=args.fwd,
        add_fwd_len=tuple(args.flen),
        add_rev=args.rev,
        add_rev_len=tuple(args.rlen),
    )

    print(PRIMER_FMT.format("dir", "tm", "ttm", "dg", "pen", "seq"))
    print(fwd)
    print(rev)


def parse_args(args: List[str]) -> argparse.Namespace:
    """Parse command line parameters.

    Created and based on an example from pyscaffold:
    https://github.com/pyscaffold/pyscaffold/blob/master/src/pyscaffold/templates/skeleton.template

    Args:
        args ([str]): List of parameters as strings

    Returns:
        `argparse.Namespace`: command line parameters namespace
    """

    parser = argparse.ArgumentParser(
        description="""Create PCR primers for a DNA sequence.

Logs the FWD and REV primer with columns:
dir, tm, ttm, dg, pen, seq

Where:
dir = FWD or REV.
tm  = Melting temperature of the annealing/binding part of the primer (Celsius).
ttm = The total melting temperature of the primer with added seq (Celsius).
dg  = The minimum free energy of the primer (kcal/mol).
pen = The primer's penalty score. Lower is better.
seq = The sequence of the primer in the 5' to the 3' direction.
""",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument("seq", type=str, metavar="SEQ", help="DNA sequence")
    parser.add_argument(
        "-f",
        "--fwd",
        type=str,
        help="additional sequence to add to FWD primer (5' to 3')",
        default="",
        metavar="SEQ",
    )
    parser.add_argument(
        "-fl",
        "--flen",
        type=int,
        nargs=2,
        help="space separated min-max range for the length to add from 'add_fwd' (5' to 3')",
        default=[-1, -1],
        metavar="INT",
    )
    parser.add_argument(
        "-r",
        "--rev",
        type=str,
        help="additional sequence to add to REV primer (5' to 3')",
        default="",
        metavar="SEQ",
    )
    parser.add_argument(
        "-rl",
        "--rlen",
        type=int,
        nargs=2,
        help="space separated min-max range for the length to add from 'add_rev' (5' to 3')",
        default=[-1, -1],
        metavar="INT",
    )
    parser.add_argument(
        "--version", action="version", version="seqfold {ver}".format(ver=__version__)
    )

    return parser.parse_args(args)


if __name__ == "__main__":
    run()
