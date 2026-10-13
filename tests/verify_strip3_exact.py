#!/usr/bin/env python3
"""Derive the exact Grundy value of the full 3xn strip from a sweep file.

Inputs:
  * a sweep result produced by `scripts/strip3_ports.py sweep --n N` (exact
    values of the n//2 edge orbits B(n,j)); the remaining edge values follow
    from the proven reflection lemmas (b(n,j) = b(n,n+1-j) horizontally and
    top row = bottom row vertically);
  * the a-sequence: `a(t)` = value of the right-end-middle-gap strip, taken
    from data/strip3_survey.jsonl for t<=28 and from data/strip3_a_ext.json
    for larger t;
  * the middle-move split formula a(i-1) xor a(n-i) for a move in the middle
    row at column i (notes/08).

The verifier reconstructs every option value independently, checks the sweep
masks against coordinates, and prints the mex.  With --record it appends the
value to data/raw/strip3_3xn.tsv and data/grid_values.csv if consistent.

Run: python3 tests/verify_strip3_exact.py --n 30 --sweep /tmp/sweep_30.json
"""

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests"))
import test_strip3_ports as base


def encode(cells):
    """Column-major encoding (interface definition, re-typed here)."""
    mask = 0
    for r, c in cells:
        mask |= 1 << (3 * (c - 1) + (r - 1))
    return mask


def load_a():
    survey = [json.loads(line) for line in
              (ROOT / "data/strip3_survey.jsonl").read_text().splitlines()
              if line.strip()]
    a = {0: 0}
    for record in survey:
        a[record["n"]] = record["end_gap"]
    ext = ROOT / "data/strip3_a_ext.json"
    if ext.exists():
        for key, value in json.loads(ext.read_text())["values"].items():
            a[int(key)] = value
    return a


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--sweep", type=Path, required=True)
    parser.add_argument("--record", action="store_true",
                        help="append to data/raw/strip3_3xn.tsv and data/grid_values.csv")
    args = parser.parse_args()
    n = args.n

    sweep = json.loads(args.sweep.read_text())
    if sweep.get("kind") != "strip3_sweep" or sweep.get("n") != n:
        raise ValueError("wrong sweep file")
    values = {}
    for entry in sweep["results"]:
        if entry["value"] is None:
            raise ValueError(f"column {entry['column']} is UNKNOWN; cannot derive f({n})")
        values[entry["column"]] = entry["value"]
        expect = encode(base.cells_of("B", n, entry["column"]))
        if int(entry["mask"], 16) != expect:
            raise ValueError(f"mask mismatch for column {entry['column']}")

    orbits = (n + 1) // 2
    if sorted(values) != list(range(1, orbits + 1)):
        raise ValueError(f"orbit coverage is {sorted(values)}, need 1..{orbits}")

    # extend to every top-row value by the horizontal reflection lemma
    top = {}
    for j in range(1, n + 1):
        top[j] = values[j] if j <= orbits else values[n + 1 - j]

    a = load_a()
    missing = [t for i in range(1, n + 1) for t in (i - 1, n - i) if t not in a]
    if missing:
        raise ValueError(f"missing a-values {sorted(set(missing))}")

    middle = [a[i - 1] ^ a[n - i] for i in range(1, n + 1)]
    options = middle + [top[j] for j in range(1, n + 1)] + [top[j] for j in range(1, n + 1)]
    present = set(options)
    g = 0
    while g in present:
        g += 1
    print(f"f({n}) = {g}   [middle {len(middle)} values, both edge rows {n} each]")
    print(f"  zero options: {sum(v == 0 for v in options)}; "
          f"first-move winners: "
          f"{'none (P-position)' if g == 0 else 'some'}")

    if args.record:
        tsv = ROOT / "data/raw/strip3_3xn.tsv"
        lines = [l for l in tsv.read_text().splitlines() if l.strip()]
        rows = {}
        for line in lines:
            parts = line.split()
            if len(parts) == 3 and parts[0].isdigit():
                rows[int(parts[1])] = line
        if n in rows:
            print(f"  {n} already recorded: {rows[n]}")
        else:
            with tsv.open("a") as handle:
                handle.write(f"3\t{n}\t{g}\n")
            with (ROOT / "data/grid_values.csv").open("a") as handle:
                handle.write(f"3,{n},{g}\n")
            print(f"  recorded 3,{n},{g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
