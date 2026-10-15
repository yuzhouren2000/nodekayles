#!/usr/bin/env python3
"""Verify the exact value of the full 3x33 strip from the language-survey list.

Inputs:
  * data/strip3_language_survey_33.jsonl  (per-cell first-move values emitted by
    `src/cfast/strip3_language survey 33`; column-major v = 3c+r);
  * the a-sequence a(t) = value of the right-end-middle-gap strip, from
    data/strip3_survey.jsonl (t<=28) and data/strip3_a_ext.json (t=29..32);
  * the proven middle-move split formula a(i-1) xor a(n-i) for a move in the
    middle row at column i (notes/08).

Independently re-derives:
  * every middle-row option value from the a-sequence alone (old solver data);
  * board reflection symmetry of the full option list;
  * the mex of all 3n option values -> f(n), and the winning first moves.

Run: python3 tests/verify_strip3_language_33.py [--n 33] [--survey FILE]
"""

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_N = 33


def load_a():
    a = {0: 0}
    for line in (ROOT / "data/strip3_survey.jsonl").read_text().splitlines():
        if line.strip():
            record = json.loads(line)
            a[record["n"]] = record["end_gap"]
    ext = json.loads((ROOT / "data/strip3_a_ext.json").read_text())["values"]
    for key, value in ext.items():
        a[int(key)] = value
    return a


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--n", type=int, default=DEFAULT_N)
    parser.add_argument("--survey", type=Path,
                        default=ROOT / "data/strip3_language_survey_33.jsonl")
    args = parser.parse_args()
    N = args.n

    survey = json.loads(args.survey.read_text())
    if survey["n"] != N:
        raise SystemExit(f"survey is for n={survey['n']}, expected {N}")
    options = survey["options"]
    assert len(options) == 3 * N, "option count"

    a = load_a()
    problems = []

    # 1) middle-row options recomputed from the old solver's a-sequence
    checked, pending = 0, []
    for c in range(N):
        want = a.get(c)
        other = a.get(N - 1 - c)
        got = options[3 * c + 1]
        if want is None or other is None:
            pending.append(c + 1)
            continue
        expect = want ^ other
        if got != expect:
            problems.append(f"middle column {c + 1}: survey {got}, "
                            f"a({c})^a({N - 1 - c})={expect}")
        else:
            checked += 1

    # 2) board reflection symmetry of the whole option list
    for c in range(N):
        for r in range(3):
            here = options[3 * c + r]
            mirror = options[3 * (N - 1 - c) + (2 - r)]
            if here != mirror:
                problems.append(f"reflection mismatch at (r{r + 1},c{c + 1}): "
                                f"{here} vs {mirror}")

    # 3) mex over all option values
    present = set(options)
    g = 0
    while g in present:
        g += 1

    winners = [(v // 3 + 1, v % 3 + 1) for v in range(3 * N) if options[v] == 0]
    print(f"survey mex = {g}; recorded grundy = {survey['grundy']}")
    if g != survey["grundy"]:
        problems.append(f"mex {g} != recorded grundy {survey['grundy']}")
    print(f"middle options cross-checked against the a-sequence: {checked}/{N}"
          + (f"; pending a-values for columns {pending}" if pending else ""))
    print(f"winning first moves (value-0 options): {winners}")
    print(f"option value set: {sorted(set(options))}")

    if problems:
        print("PROBLEMS:")
        for p in problems:
            print("  " + p)
        return 1
    print(f"f({N}) = {g} verified from the survey certificate")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
