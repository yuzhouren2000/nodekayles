"""Search arbitrary fixed-pairing certificates, with checkable refutations.

A candidate pair A={u,v} must contain nonadjacent playable vertices. Two
selected pairs A,B must be disjoint, and each must intersect the other's
closed neighbourhood union in either zero or two vertices. A compatible
perfect matching is exactly a certificate satisfying (C1)-(C3).

The search branches on the uncovered vertex with the fewest allowed pairs.
UNSAT includes every branch and a dead-vertex leaf, as a proof DAG. The
independent verifier in tests/verify_fixed_pairing.py reconstructs the
constraints and checks branch coverage without running this search.

Run: python3 scripts/fixed_pairing.py 5 8 --output data/proofs/fixed_5x8.json
No external SAT library is required. All vertex indices are row-major, 0-based.
"""

import argparse
import json
import time
from pathlib import Path

FORMAT = "nodekayles-fixed-pairing-v1"


def bits(mask):
    while mask:
        bit = mask & -mask
        mask ^= bit
        yield bit.bit_length() - 1


class SearchLimit(Exception):
    pass


def rotation_pairing(rows, cols, mask):
    """A sufficient structural certificate; None says only rotation fails.

    The 180-degree map must preserve the playable set and pair no cell with
    itself or an orthogonal neighbour. It is then an induced-graph automorphism.
    """
    if rows <= 0 or cols <= 0 or mask < 0 or mask >> (rows * cols):
        raise ValueError("invalid board or playable mask")
    pairs = []
    for a in bits(mask):
        b = rows * cols - 1 - a
        if not mask >> b & 1 or a == b:
            return None
        ar, ac = divmod(a, cols)
        br, bc = divmod(b, cols)
        if abs(ar - br) + abs(ac - bc) == 1:
            return None
        if a < b:
            pairs.append((a, b))
    return pairs


def solve_pairing(rows, cols, mask=None, *, node_limit=100_000, timeout=30):
    """Return SAT + pairs, UNSAT + proof, or UNKNOWN when a limit is reached.

    mask describes the induced graph of currently playable cells, not stones.
    None means the full board. An empty position has the empty certificate.
    """
    if rows <= 0 or cols <= 0:
        raise ValueError("board dimensions must be positive")
    full = (1 << (rows * cols)) - 1
    if mask is None:
        mask = full
    if mask < 0 or mask & ~full:
        raise ValueError("playable mask extends outside the board")
    if node_limit < 0 or timeout < 0:
        raise ValueError("search limits must be nonnegative")
    start = time.monotonic()
    result = {
        "format": FORMAT,
        "rows": rows,
        "cols": cols,
        "playable": list(bits(mask)),
    }
    visits = 0
    nodes = []
    pair_count = 0

    def finish(status, **payload):
        return {
            **result,
            "status": status,
            **payload,
            "stats": {
                "candidate_pairs": pair_count,
                "search_nodes": visits,
                "proof_nodes": len(nodes) if status == "UNSAT" else 0,
                "seconds": round(time.monotonic() - start, 6),
            },
        }

    if mask.bit_count() % 2:
        return finish("UNSAT", proof={"kind": "odd-order"})
    if not mask:
        return finish("SAT", pairs=[])

    def check_time():
        if time.monotonic() - start >= timeout:
            raise SearchLimit("time limit")

    try:
        check_time()
        closed = [0] * (rows * cols)
        for v in bits(mask):
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
            closed[v] = nb & mask
        vertices = result["playable"]
        pairs = [
            (a, b)
            for i, a in enumerate(vertices)
            for b in vertices[i + 1 :]
            if not (closed[a] >> b & 1)
        ]
        pair_count = len(pairs)
        pair_masks = [(1 << a) | (1 << b) for a, b in pairs]
        unions = [closed[a] | closed[b] for a, b in pairs]
        incident = [0] * (rows * cols)
        compatible = [0] * pair_count
        for i, (a, b) in enumerate(pairs):
            incident[a] |= 1 << i
            incident[b] |= 1 << i
        for i in range(pair_count):
            check_time()
            for j in range(i + 1, pair_count):
                if pair_masks[i] & pair_masks[j]:
                    continue
                if unions[i] & pair_masks[j] not in (0, pair_masks[j]):
                    continue
                if unions[j] & pair_masks[i] not in (0, pair_masks[i]):
                    continue
                compatible[i] |= 1 << j
                compatible[j] |= 1 << i

        dead = {}

        def dfs(remaining, allowed):
            nonlocal visits
            check_time()
            if visits >= node_limit:
                raise SearchLimit("node limit")
            visits += 1
            if not remaining:
                return [], None
            state = remaining, allowed
            if state in dead:
                return None, dead[state]
            vertex, choices = None, None
            for v in bits(remaining):
                candidates = incident[v] & allowed
                if choices is None or candidates.bit_count() < choices.bit_count():
                    vertex, choices = v, candidates
                if not candidates:
                    break
            branches = []
            for i in bits(choices):
                certificate, child = dfs(
                    remaining ^ pair_masks[i], allowed & compatible[i]
                )
                if certificate is not None:
                    return [pairs[i], *certificate], None
                a, b = pairs[i]
                partner = b if a == vertex else a
                branches.append([partner, child])
            node_id = len(nodes)
            nodes.append({"vertex": vertex, "branches": branches})
            dead[state] = node_id
            return None, node_id

        certificate, root = dfs(mask, (1 << pair_count) - 1)
        if certificate is not None:
            return finish("SAT", pairs=certificate)
        return finish("UNSAT", proof={"kind": "branching", "root": root, "nodes": nodes})
    except SearchLimit as exc:
        return finish("UNKNOWN", reason=str(exc))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rows", type=int)
    parser.add_argument("cols", type=int)
    parser.add_argument("--mask", type=lambda x: int(x, 0), help="playable bitmask")
    parser.add_argument("--node-limit", type=int, default=100_000)
    parser.add_argument("--timeout", type=float, default=30)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = solve_pairing(
        args.rows, args.cols, args.mask,
        node_limit=args.node_limit, timeout=args.timeout,
    )
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"{args.rows}x{args.cols}: {result['status']} {result['stats']}")
    if args.output:
        print(f"wrote {args.output}")
    return 2 if result["status"] == "UNKNOWN" else 0


if __name__ == "__main__":
    raise SystemExit(main())
