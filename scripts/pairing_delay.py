"""Prove that fixed pairing cannot be forced within two response rounds.

After a selected first attack a, for every legal response b the remaining
graph has no fixed pairing, and some legal next attack c makes EVERY legal
response d leave a graph without a fixed pairing. Each negative subproblem
contains a standalone refutation; no SG values are used.

Run: python3 scripts/pairing_delay.py 5 8 --first 2 1 \
        --output data/proofs/pairing_delay_5x8.json.gz
The first move uses 0-based coordinates. Gzip is optional, inferred by suffix.
"""

import argparse
import gzip
import json
from pathlib import Path

from fixed_pairing import bits, solve_pairing


def build_obstruction(rows, cols, first):
    if rows <= 0 or cols <= 0 or not 0 <= first < rows * cols:
        raise ValueError("invalid board or first move")
    full = (1 << (rows * cols)) - 1
    closed = []
    for v in range(rows * cols):
        r, c = divmod(v, cols)
        nb = 1 << v
        if r:
            nb |= 1 << (v - cols)
        if r + 1 < rows:
            nb |= 1 << (v + cols)
        if c:
            nb |= 1 << (v - 1)
        if c + 1 < cols:
            nb |= 1 << (v + 1)
        closed.append(nb)
    cache = {}
    needed = set()

    def fixed(mask):
        if mask not in cache:
            cache[mask] = solve_pairing(rows, cols, mask)
        result = cache[mask]
        if result["status"] == "UNKNOWN":
            raise RuntimeError("pairing search limit reached; obstruction is incomplete")
        return result["status"] == "SAT"

    if fixed(full):
        return {"status": "NO_OBSTRUCTION", "reason": "initial graph already has a pairing"}
    needed.add(full)
    after_first = full & ~closed[first]
    branches = []
    for reply in bits(after_first):
        after_reply = after_first & ~closed[reply]
        if fixed(after_reply):
            return {"status": "NO_OBSTRUCTION", "reason": "pairing after first response", "reply": reply}
        needed.add(after_reply)
        for attack in bits(after_reply):
            after_attack = after_reply & ~closed[attack]
            tails = []
            for second_reply in bits(after_attack):
                tail = after_attack & ~closed[second_reply]
                if fixed(tail):
                    break
                tails.append({"reply": second_reply, "refutation": str(tail)})
            else:
                needed.update(int(t["refutation"]) for t in tails)
                branches.append({
                    "reply": reply,
                    "after_reply_refutation": str(after_reply),
                    "blocking_attack": attack,
                    "second_replies": tails,
                })
                break
        else:
            return {"status": "NO_OBSTRUCTION", "reason": "no blocking second attack", "reply": reply}
    return {
        "format": "nodekayles-pairing-delay-v1",
        "status": "OBSTRUCTION",
        "rows": rows,
        "cols": cols,
        "first": first,
        "initial_refutation": str(full),
        "branches": branches,
        "refutations": {str(mask): cache[mask] for mask in sorted(needed)},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rows", type=int)
    parser.add_argument("cols", type=int)
    parser.add_argument("--first", type=int, nargs=2, required=True, metavar=("ROW", "COL"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    r, c = args.first
    if not 0 <= r < args.rows or not 0 <= c < args.cols:
        parser.error("first move is outside the board")
    result = build_obstruction(args.rows, args.cols, r * args.cols + c)
    encoded = (json.dumps(result, separators=(",", ":")) + "\n").encode()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.output.suffix == ".gz":
        args.output.write_bytes(gzip.compress(encoded, mtime=0))
    else:
        args.output.write_bytes(encoded)
    print(f"{args.rows}x{args.cols}, first {args.first}: {result['status']}")
    if result["status"] == "OBSTRUCTION":
        print(f"{len(result['branches'])} first replies; {len(result['refutations'])} shared refutations")
    print(f"wrote {args.output}")
    return 0 if result["status"] == "OBSTRUCTION" else 2


if __name__ == "__main__":
    raise SystemExit(main())
