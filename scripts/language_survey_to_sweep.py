#!/usr/bin/env python3
"""Convert a language-solver survey list into the standard strip3_sweep format.

`strip3_language survey N` prints one JSON line with the value of every first
move of the empty 3xN board (column-major cell order v = 3c+r).  The verifier
tests/verify_strip3_exact.py consumes a "strip3_sweep" file holding the n//2
distinct edge orbits B(n,j); the top-row option at column j is exactly that
value, and the sweep masks are re-checked against base.cells_of("B", n, j).

Run: python3 scripts/language_survey_to_sweep.py --survey FILE --output FILE
"""

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests"))
import test_strip3_ports as base


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--survey", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source", default="src/cfast/strip3_language survey")
    args = parser.parse_args()

    survey = json.loads(args.survey.read_text())
    n = survey["n"]
    options = survey["options"]
    if len(options) != 3 * n:
        raise SystemExit("option count does not match n")

    orbits = (n + 1) // 2
    results = []
    for j in range(1, orbits + 1):
        value = options[3 * (j - 1)]          # top row of column j
        mask = 0
        for r, c in base.cells_of("B", n, j):
            mask |= 1 << (3 * (c - 1) + (r - 1))
        results.append({"column": j, "mask": f"{mask:x}", "value": value})

    out = {
        "version": 1,
        "kind": "strip3_sweep",
        "n": n,
        "columns": list(range(1, orbits + 1)),
        "ht_bits": None,
        "load": None,
        "prealloc": None,
        "threads": None,
        "timeout_s": None,
        "status": "ok",
        "returncode": 0,
        "elapsed_s": None,
        "union_memo": None,
        "memo_at_exit": None,
        "binary_sha256": None,
        "results": results,
        "note": (f"edge orbits extracted from {args.source} (language solver); "
                 "grundy and middle options are re-derived by verify_strip3_exact"),
    }
    args.output.write_text(json.dumps(out, indent=1) + "\n")
    print(f"wrote {args.output}: {orbits} orbits, n={n}")


if __name__ == "__main__":
    raise SystemExit(main())
