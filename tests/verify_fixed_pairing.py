"""Independently verify fixed-pairing certificates and exhaustive refutations.

This module imports neither the search nor the SG solver. Geometry and pair
compatibility are reconstructed with sets. A refutation must branch on EVERY
currently allowed partner of its chosen vertex; an empty branch list proves
that vertex cannot be covered. Repeated proof nodes must have the same state.

Run: python3 tests/verify_fixed_pairing.py data/proofs/fixed_*.json
"""

import json
import sys
from pathlib import Path


def neighbourhoods(rows, cols, vertices):
    out = {}
    for v in vertices:
        r, c = divmod(v, cols)
        neighbours = {v}
        for dr, dc in ((-1, 0), (1, 0), (0, -1), (0, 1)):
            rr, cc = r + dr, c + dc
            if 0 <= rr < rows and 0 <= cc < cols:
                neighbours.add(rr * cols + cc)
        out[v] = neighbours & vertices
    return out


def verify_pairs(rows, cols, vertices, pairs):
    """Check the actual response and invariant after every possible first move."""
    vertices = set(vertices)
    nb = neighbourhoods(rows, cols, vertices)
    phi = {}
    for pair in pairs:
        if len(pair) != 2:
            raise ValueError("a pair must have two vertices")
        a, b = pair
        if a == b or a not in vertices or b not in vertices or a in phi or b in phi:
            raise ValueError("pairs do not form a matching of playable vertices")
        phi[a], phi[b] = b, a
    if set(phi) != vertices:
        raise ValueError("some playable vertices are unpaired")
    for a in vertices:
        b = phi[a]
        if b in nb[a]:
            raise ValueError("the prescribed response is illegal")
        removed = nb[a] | nb[b]
        if {phi[v] for v in removed} != removed:
            raise ValueError("the two moves do not restore the pairing invariant")
    return True


def verify_result(data):
    if data.get("format") != "nodekayles-fixed-pairing-v1":
        raise ValueError("unknown proof format")
    rows, cols = data["rows"], data["cols"]
    if type(rows) is not int or type(cols) is not int or rows <= 0 or cols <= 0:
        raise ValueError("invalid board dimensions")
    raw_vertices = data["playable"]
    if any(type(v) is not int or not 0 <= v < rows * cols for v in raw_vertices):
        raise ValueError("invalid playable vertex")
    vertices = frozenset(raw_vertices)
    if len(vertices) != len(raw_vertices):
        raise ValueError("duplicate playable vertex")
    status = data["status"]
    if status == "SAT":
        return verify_pairs(rows, cols, vertices, data["pairs"])
    if status != "UNSAT":
        raise ValueError("UNKNOWN is not a certificate or a refutation")
    proof = data["proof"]
    if proof["kind"] == "odd-order":
        if len(vertices) % 2 != 1:
            raise ValueError("odd-order obstruction on an even-order graph")
        return True
    if proof["kind"] != "branching":
        raise ValueError("unknown refutation kind")
    nb = neighbourhoods(rows, cols, vertices)
    ordered = sorted(vertices)
    candidates = frozenset(
        (a, b)
        for i, a in enumerate(ordered)
        for b in ordered[i + 1 :]
        if b not in nb[a]
    )
    unions = {pair: nb[pair[0]] | nb[pair[1]] for pair in candidates}

    def compatible(a, b):
        aset, bset = set(a), set(b)
        return (
            not aset & bset
            and len(unions[a] & bset) in (0, 2)
            and len(unions[b] & aset) in (0, 2)
        )

    nodes = proof["nodes"]
    verified = {}
    active = set()

    def check(node_id, remaining, allowed):
        if type(node_id) is not int or not 0 <= node_id < len(nodes):
            raise ValueError("invalid proof node reference")
        state = remaining, allowed
        if node_id in active:
            raise ValueError("cycle in refutation")
        if node_id in verified:
            if verified[node_id] != state:
                raise ValueError("shared node is used for different states")
            return
        if not remaining:
            raise ValueError("cannot refute the empty matching problem")
        node = nodes[node_id]
        v = node["vertex"]
        if type(v) is not int or v not in remaining:
            raise ValueError("branch vertex is not uncovered")
        choices = {
            (b if a == v else a): (a, b)
            for a, b in allowed
            if v == a or v == b
        }
        branches = node["branches"]
        partners = [branch[0] for branch in branches]
        if len(set(partners)) != len(partners) or set(partners) != set(choices):
            raise ValueError("refutation does not cover every allowed partner")
        active.add(node_id)
        for partner, child in branches:
            pair = choices[partner]
            next_allowed = frozenset(p for p in allowed if compatible(pair, p))
            check(child, remaining - set(pair), next_allowed)
        active.remove(node_id)
        verified[node_id] = state

    check(proof["root"], vertices, candidates)
    return True


def main():
    if not sys.argv[1:]:
        print(__doc__)
        return 2
    for name in sys.argv[1:]:
        data = json.loads(Path(name).read_text())
        verify_result(data)
        print(f"verified {name}: {data['status']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
