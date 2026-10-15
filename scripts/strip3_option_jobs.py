#!/usr/bin/env python3
"""Generate first-move option jobs for re-deriving one board's SG with the
old solver's `masks` mode (option decomposition).

A position's options are exactly {position minus the closed neighbourhood of
one occupied cell}.  Feeding every option of a board as its own `masks` job
keeps the old solver's shared memo (work is reused across jobs) and
parallelises over jobs even when a single board search would not.

Board families use the same coordinate definitions as
tests/test_strip3_ports.py:

  A(t)       full 3xt board minus the middle cell of the last column;
  B(n, j)    full 3xn board minus the closed neighbourhood of (1, j)
             (the edge-orbit representative of scripts/strip3_ports.py).

Outputs a jobs file for `src/cfast/strip3 masks` (one "width mask" line per
option, masks in hex) and a JSON span index mapping board key to the
[start, stop) line range of its options.

Run:
  python3 scripts/strip3_option_jobs.py A 29 30 31 32 \
      --out data/strip3_a_options_in.txt --index data/strip3_a_options_index.json
  python3 scripts/strip3_option_jobs.py B 33 --columns 1 6 17 \
      --out /tmp/b33_jobs.txt --index /tmp/b33_index.json
"""

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests"))
import test_strip3_ports as base


def encode(cells):
    """Column-major mask (same encoding as the strip3 interface)."""
    mask = 0
    for r, c in cells:
        mask |= 1 << (3 * (c - 1) + (r - 1))
    return mask


def closed_cells(cell):
    r, c = cell
    return {(rr, cc) for rr, cc in
            ((r, c), (r - 1, c), (r + 1, c), (r, c - 1), (r, c + 1))
            if 1 <= rr <= 3 and cc >= 1}


def options_of(cells):
    """Options in increasing mask-bit order (lsb first)."""
    cells = frozenset(cells)
    order = sorted(cells, key=lambda rc: 3 * (rc[1] - 1) + (rc[0] - 1))
    return [encode(cells - closed_cells(cell)) for cell in order]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("family", choices=("A", "B"))
    parser.add_argument("n", type=int,
                        help="A: the widths t (one or more); B: the single board width")
    parser.add_argument("columns", nargs="*", type=int,
                        help="A: more widths; B: the hole columns j")
    parser.add_argument("--columns", dest="columns_opt", nargs="+", type=int,
                        help="B only: hole columns (alternative to positional)")
    parser.add_argument("--out", type=Path, required=True, help="jobs file for `strip3 masks`")
    parser.add_argument("--index", type=Path, required=True, help="span index JSON")
    args = parser.parse_args()

    columns = args.columns_opt if args.columns_opt else args.columns
    if args.family == "B" and not columns:
        parser.error("B requires at least one hole column j")

    jobs, index = [], {}
    if args.family == "A":
        for t in [args.n] + columns:
            cells = base.cells_of("A", t)
            options = options_of(cells)
            index[str(t)] = [len(jobs), len(jobs) + len(options)]
            jobs.extend((t, option) for option in options)
            print(f"A({t}): cells={len(cells)} options={len(options)}")
    else:
        n = args.n
        for j in columns:
            cells = base.cells_of("B", n, j)
            options = options_of(cells)
            index[str(j)] = [len(jobs), len(jobs) + len(options)]
            jobs.extend((n, option) for option in options)
            print(f"B({n},{j}): cells={len(cells)} options={len(options)}")

    with args.out.open("w") as handle:
        for width, option in jobs:
            handle.write(f"{width} {option:x}\n")
    args.index.write_text(json.dumps(index))
    print(f"total jobs: {len(jobs)} -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
