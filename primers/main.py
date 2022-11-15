"""Console entrypoint for creating PCR primers"""

import argparse
import json
import sys
from typing import List, Optional

from . import __version__, primers, score, Primer

"""{fwd} {tm} {tm_total} {gc} {dg} {penalty} {seq}"""


def run():
    """Entry point for the CLI"""

    args = parse_args(sys.argv[1:])
    args.run(args)


def parse_args(args: List[str]) -> argparse.Namespace:
    """Create argparser"""

    parser = argparse.ArgumentParser(
        description="""Create or score PCR primers""",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--version", action="version", version="seqfold {ver}".format(ver=__version__)
    )
    sub_parsers = parser.add_subparsers()

    # Add a sub-command for primer scoring
    parser_score = sub_parsers.add_parser("score", help="score existing primers")
    parser_score.set_defaults(run=run_score)
    parser_score.add_argument("primer", type=str, nargs="+", help="primer sequences")
    parser_score.add_argument(
        "-s",
        type=str,
        help="the sequence that was amplified",
        default="",
        metavar="SEQ",
    )
    parser_score.add_argument(
        "-t",
        type=str,
        help="sequence to check for off-target binding sites",
        default="",
        metavar="SEQ",
    )
    parser_score.add_argument(
        "-j",
        "--json",
        action=argparse.BooleanOptionalAction,
        help="write the primers to a JSON array",
    )

    # Add a sub-command for primer creation
    parser_create = sub_parsers.add_parser(
        "create", help="create primers to amplify a sequence"
    )
    parser_create.set_defaults(run=run_create)
    parser_create.add_argument(
        "seq", type=str, metavar="SEQ", help="create primers to amplify this sequence"
    )
    parser_create.add_argument(
        "-f",
        type=str,
        help="additional sequence to add to FWD primer (5' to 3')",
        default="",
        metavar="SEQ",
    )
    parser_create.add_argument(
        "-fl",
        type=int,
        nargs=2,
        help="space separated min-max range for the length to add from '-f' (5' to 3')",
        default=[-1, -1],
        metavar="INT",
    )
    parser_create.add_argument(
        "-r",
        type=str,
        help="additional sequence to add to REV primer (5' to 3')",
        default="",
        metavar="SEQ",
    )
    parser_create.add_argument(
        "-rl",
        type=int,
        nargs=2,
        help="space separated min-max range for the length to add from '-r' (5' to 3')",
        default=[-1, -1],
        metavar="INT",
    )
    parser_create.add_argument(
        "-t",
        type=str,
        help="sequence to check for off-target binding sites",
        default="",
        metavar="SEQ",
    )
    parser_create.add_argument(
        "-j",
        "--json",
        action=argparse.BooleanOptionalAction,
        help="write the primers to a JSON array",
    )

    return parser.parse_args(args)


def run_create(args: argparse.Namespace):
    fwd, rev = primers(
        args.seq,
        add_fwd=args.f,
        add_fwd_len=tuple(args.fl),  # type: ignore
        add_rev=args.r,
        add_rev_len=tuple(args.rl),  # type: ignore
        offtarget_check=args.t,
    )

    print_output(fwd, rev, args.json)


def run_score(args: argparse.Namespace):
    seq_fwd = args.primer[0]
    seq_rev = ""
    if len(args.primer) > 1:
        seq_rev = args.primer[1]
    elif len(args.primer) > 2:
        print(
            "cannot pass more than 2 primers to 'primers score' at this time. Open a GH issue if that's a desired feature"
        )
        exit(1)

    fwd, rev = score(seq_fwd, rev=seq_rev, seq=args.s, offtarget_check=args.t)

    print_output(fwd, rev, args.json)


def print_output(fwd: Primer, rev: Optional[Primer], print_json=False):
    if print_json:
        if rev:
            print(json.dumps([fwd.dict(), rev.dict()]))
        else:
            print(json.dumps([fwd.dict()]))
    else:
        table_fmt = "{:>5} {:>5} {:>5} {:>5} {:>5} {:>5}  {}"
        print(table_fmt.format("dir", "tm", "ttm", "gc", "dg", "p", "seq"))

        for p in [fwd, rev]:
            if not p:
                continue

            print(
                table_fmt.format(
                    "FWD" if p.fwd else "REV",
                    p.tm,
                    p.tm_total,
                    p.gc,
                    p.dg,
                    p.scoring.penalty,
                    p.seq,
                )
            )


if __name__ == "__main__":
    run()
