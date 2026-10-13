#!/usr/bin/env python3
"""Independent check of the 3x30 edge-orbit target file.

If data/strip3_target_30.json contains exact values for all 15 edge orbits
(j = 1..15, the top row; the bottom row is its vertical reflection and
j -> 31-j is the horizontal reflection), then, together with the verified
middle-move exclusions for n = 30, the mex argument gives the exact value
G(3x30) = 0.  This script re-derives every premise from the raw survey data
and the analysis file, and refuses to conclude anything if a premise is
missing.  A zero column would instead mean an N-position with that first
move; the script then checks the recorded option list (if any) for the
P-certificate shape (all options nonzero).

Run: python3 tests/verify_strip3_target.py \
        [--target data/strip3_target_30.json] [--n 30]
"""

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests"))
import test_strip3_ports as base  # coordinate helpers only


def load_a_d():
    survey = [json.loads(line) for line in
              (ROOT / "data/strip3_survey.jsonl").read_text().splitlines()
              if line.strip()]
    a = [0] + [record["end_gap"] for record in survey]
    d = [0] + [record["both_gaps"] for record in survey]
    return survey, a, d


def middle_option_values(n, a, d):
    """The 30 middle-move option values of F_n, with the zero-tail premise
    used exactly where the analysis file uses it (unknown a(q))."""
    known = len(a) - 1          # a(0..known) are exact
    forbidden = {t + 1 for t in range(1, len(d)) if d[t] == 0}
    values = []
    for i in range(1, n + 1):
        p, q = i - 1, n - i
        if p >= len(a) and q < len(a):
            p, q = q, p          # the screen always keeps the unknown on the right
        if p < len(a) and q < len(a):
            values.append(a[p] ^ a[q])
        elif q >= len(a) and p < len(a):
            # the screen excludes this via a witness or a zero-tail gap
            gap = q - p
            witnesses = [w for w in range(len(a))
                         if a[w] == a[p] and abs(q - w) in forbidden]
            if gap in forbidden:
                values.append(None)  # proven nonzero, exact value unknown
            elif witnesses:
                values.append(None)  # proven nonzero: a(q) != a(p)
            else:
                values.append("UNRESOLVED")
        else:
            values.append("UNRESOLVED")
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--target", type=Path,
                        default=ROOT / "data/strip3_target_30.json")
    parser.add_argument("--n", type=int, default=30)
    args = parser.parse_args()
    n = args.n

    if not args.target.exists():
        print(f"{args.target} not present; nothing to verify")
        return 1
    target = json.loads(args.target.read_text())
    if target.get("kind") != "strip3_target" or target.get("n") != n:
        raise ValueError("target file has the wrong kind or n")

    survey, a, d = load_a_d()
    if any(d[t] == 0 and a[k] == a[k + t + 1]
           for t in range(1, len(d)) for k in range(len(a) - t - 1)):
        raise ValueError("zero-tail theorem violated by the stored sequences")

    analysis = json.loads((ROOT / "data/strip3_analysis.json").read_text())
    screen = {s["n"]: s for s in analysis["middle_screens"]}[n]
    if screen["status"] != "no_middle_win":
        raise ValueError("middle screen for this n is not complete")
    if screen["unresolved_split_lengths"]:
        raise ValueError("middle screen still has unresolved splits")
    values = middle_option_values(n, a, d)
    if any(v == "UNRESOLVED" for v in values):
        raise ValueError("an unknown a(q) is not covered by the screen")
    if any(v == 0 for v in values):
        raise ValueError("a middle move has value 0; f(30) != 0")

    # masks and orbit coverage from the target file
    by_column = {}
    for entry in target.get("results", []):
        if entry["value"] is not None:
            by_column[entry["column"]] = entry

    def mask_ok(column, mask_hex):
        n = target["n"]
        left = (1 << (3 * n)) - 1
        for c in (column - 1, column, column + 1):
            if 1 <= c <= n:
                left &= ~(1 << (3 * (c - 1)))
        left &= ~(1 << (3 * (column - 1) + 1))
        return int(mask_hex, 16) == left

    for j, entry in by_column.items():
        if not mask_ok(j, entry["mask"]):
            raise ValueError(f"mask for column {j} does not match the geometry")

    orbits = range(1, n // 2 + 1)
    covered = [j for j in orbits if j in by_column or (n + 1 - j) in by_column]
    zeros = [j for j in by_column if by_column[j]["value"] == 0]
    print(f"middle moves: all {n} excluded (screen, {len(values)} values, "
          f"{sum(v is None for v in values)} by the zero-tail premise)")
    print(f"edge orbits covered: {len(covered)}/{len(orbits)} "
          f"({sorted(by_column)})")
    if zeros:
        for j in zeros:
            entry = by_column[j]
            children = entry.get("children") or []
            if children and all(c["value"] not in (None, 0) for c in children):
                print(f"B({n},{j}) = 0 with {len(children)} nonzero options: "
                      f"3x{n} is an N-position, winning first move (1,{j}) "
                      f"(or its reflections)")
            else:
                print(f"B({n},{j}) = 0 but the option certificate is "
                      f"incomplete; run the target children for column {j}")
        return 0
    if len(covered) == len(list(orbits)):
        print(f"all {len(covered)} edge orbits are nonzero and every middle "
              f"move is excluded:")
        print(f"    every first move of 3x{n} has a nonzero successor, so "
              f"mex = 0 and G(3x{n}) = 0 (P-position)")
        return 0
    print("target incomplete; no outcome claimed")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
