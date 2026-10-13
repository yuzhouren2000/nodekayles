"""Extract first-response structure from the exact solver on a small board.

For each orbit of first moves, evaluate every legal reply and independently
search the remaining induced graph for an arbitrary fixed pairing. SAT is
checked against SG=0; UNKNOWN aborts instead of becoming a negative finding.

Run: python3 scripts/response_structure.py 5 8 --output data/responses_5x8.json.gz
Indices and printed coordinates are 0-based. This is a finite computation,
not a proof of a strategy invariant under arbitrarily many future moves.
"""

import argparse
import gzip
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nodekayles import Solver
from fixed_pairing import bits, rotation_pairing, solve_pairing


def representatives(rows, cols):
    for r in range((rows + 1) // 2):
        for c in range((cols + 1) // 2):
            if rows == cols and r > c:
                continue
            yield r * cols + c


def analyze_responses(rows, cols):
    solver = Solver()
    nb = solver.neighbours(rows, cols)
    full = (1 << (rows * cols)) - 1
    closed = [(1 << v) | nb[v] for v in range(rows * cols)]
    value = solver.grid_value(rows, cols)
    table = []
    matching_cache = {}
    for first in representatives(rows, cols):
        rest = full & ~closed[first]
        replies = []
        pairings = []
        rotation_replies = []
        for reply in bits(rest):
            tail = rest & ~closed[reply]
            g = solver.value(rows, cols, tail)
            if tail not in matching_cache:
                matching_cache[tail] = solve_pairing(rows, cols, tail)
            result = matching_cache[tail]
            if result["status"] == "UNKNOWN":
                raise RuntimeError("pairing search limit reached; survey is incomplete")
            if result["status"] == "SAT":
                if g != 0:
                    raise AssertionError("a fixed pairing must imply SG=0")
                pairings.append({"reply": reply, "certificate": result})
            rotation = rotation_pairing(rows, cols, tail)
            if rotation is not None:
                rotation_replies.append({"reply": reply, "pairs": rotation})
            replies.append({"reply": reply, "grundy": g})
        table.append({
            "first": first,
            "after_first_grundy": solver.value(rows, cols, rest),
            "replies": replies,
            "fixed_pairings": pairings,
            "rotation_replies": rotation_replies,
        })
        wins = [divmod(t["reply"], cols) for t in replies if t["grundy"] == 0]
        paired = [divmod(t["reply"], cols) for t in pairings]
        print(f"first {divmod(first, cols)}: {len(wins)} winning replies; fixed pairing {paired}", flush=True)
    return {
        "format": "nodekayles-response-survey-v1",
        "rows": rows,
        "cols": cols,
        "grundy": value,
        "coordinate_convention": "0-based row-major vertex indices",
        "first_move_orbits": table,
        "pairing_results": {str(mask): result for mask, result in matching_cache.items()},
        "solver_stats": solver.stats(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rows", type=int)
    parser.add_argument("cols", type=int)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.rows <= 0 or args.cols <= 0:
        parser.error("board dimensions must be positive")
    result = analyze_responses(args.rows, args.cols)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    encoded = (json.dumps(result, separators=(",", ":")) + "\n").encode()
    if args.output.suffix == ".gz":
        args.output.write_bytes(gzip.compress(encoded, mtime=0))
    else:
        args.output.write_bytes(encoded)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
