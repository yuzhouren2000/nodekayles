#!/usr/bin/env python3
"""Independent raw-graph checks of composition and the concrete replacement."""

from copy import deepcopy
from functools import cache
from itertools import combinations
import json
from pathlib import Path
import random
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from kayles_boundary_algebra import BoundaryAlgebra, FORK, PORTS, SQUARE
from strip3_history_language import decode
from verify_kayles_boundary import verify


@cache
def raw_sg(graph):
    """Vertex-set recursion; deliberately no component or twin reductions."""
    if not graph:
        return 0
    edges = dict(graph)
    outcomes = set()
    for v, neighbours in graph:
        remaining = {u for u in edges if u != v and u not in neighbours}
        child = tuple((u, tuple(w for w in edges[u] if w in remaining))
                      for u in sorted(remaining))
        outcomes.add(raw_sg(child))
    result = 0
    while result in outcomes:
        result += 1
    return result


def glue(a, la, b, lb):
    """Build the actual union graph using vertex names, never interface forms."""
    def names(graph, labels, side):
        return [f"port{label}" if label else f"{side}{v}"
                for v, label in enumerate(labels)]
    na, nb = names(a, la, "a"), names(b, lb, "b")
    edges = {v: set() for v in na + nb}
    for graph, ns in ((a, na), (b, nb)):
        for v, adjacency in enumerate(graph):
            edges[ns[v]].update(ns[u] for u in range(len(graph)) if adjacency >> u & 1)
    return tuple((v, tuple(sorted(edges[v]))) for v in sorted(edges))


def graphs(n):
    edges = tuple(combinations(range(n), 2))
    for code in range(1 << len(edges)):
        a = [0] * n
        for i, (u, v) in enumerate(edges):
            if code >> i & 1:
                a[u] |= 1 << v
                a[v] |= 1 << u
        yield tuple(a)


def labels(n, k):
    return tuple(range(1, k + 1)) + (0,) * (n - k)


def grid_graph(mask):
    cells = [i for i in range(mask.bit_length()) if mask >> i & 1]
    coords = [(i % 3, i // 3) for i in cells]
    graph = tuple(sum(1 << u for u, (s, d) in enumerate(coords)
                      if abs(r - s) + abs(c - d) == 1) for r, c in coords)
    return graph, cells


def rooted_isomorphic(a, b):
    if len(a) != len(b):
        return False
    from itertools import permutations
    for tail in permutations(range(1, len(a))):
        p = (0,) + tail
        if all(bool(a[i] >> j & 1) == bool(b[p[i]] >> p[j] & 1)
               for i in range(len(a)) for j in range(len(a))):
            return True
    return False


def reachable_examples():
    for mask, model in ((2484, SQUARE), (189, FORK)):
        code = decode(mask)
        n, history, target = code.witness()
        position = {(r, c) for c in range(n) for r in range(3)}
        played = set()
        for r, c in history:
            assert all(abs(r - s) + abs(c - d) != 1 for s, d in played)
            assert (r, c) in position
            played.add((r, c))
            position -= {(r, c), (r - 1, c), (r + 1, c), (r, c - 1), (r, c + 1)}
        desired = {(i % 3, i // 3) for i in range(target.bit_length()) if target >> i & 1}
        assert desired <= position
        assert all(abs(r - s) + abs(c - d) != 1
                   for r, c in desired for s, d in position - desired)
        graph, _ = grid_graph(mask)
        assert rooted_isomorphic(model, graph)
    return 2


def certificate_checks():
    proof = json.loads((ROOT / "data/proofs/kayles_boundary_replacement.json").read_text())
    result = verify(proof)
    # Reject a omitted internal safe move, wrong detached nimber and forged root.
    bad = []
    item = deepcopy(proof)
    item["nodes"][-1][3].pop()
    bad.append(item)
    item = deepcopy(proof)
    item["nodes"][-1][3][0][1][1] ^= 1
    bad.append(item)
    item = deepcopy(proof)
    item["states"][1][0] -= 1
    bad.append(item)
    for item in bad:
        try:
            verify(item)
        except (AssertionError, KeyError):
            continue
        raise AssertionError("damaged certificate accepted")
    return result


def main():
    algebra = BoundaryAlgebra()
    square, fork = algebra.state(SQUARE, PORTS), algebra.state(FORK, PORTS)
    assert square == fork
    invalid = (((1,), (1,)), ((2, 0), (1, 0)), ((2, 1), (1, 1)),
               ((0,), (-1,)), ((0,), (1, 2)), ((0,), (True,)))
    for graph, boundary in invalid:
        try:
            algebra.state(graph, boundary)
        except ValueError:
            continue
        raise AssertionError("invalid graph or boundary accepted")
    checks = 0
    # Exhaust all small simple graphs, including disconnected contexts and
    # differing adjacency between the two boundary vertices.
    for k in (1, 2, 3):
        samples = [(a, labels(n, k)) for n in range(k, 5) for a in graphs(n)]
        for a, la in samples:
            sa = algebra.state(a, la)
            for b, lb in samples:
                actual = raw_sg(glue(a, la, b, lb))
                assert algebra.glued_value(sa, algebra.state(b, lb)) == actual
                checks += 1
    random_checks = 0
    rng = random.Random(7318)
    for _ in range(1200):
        k = rng.randrange(1, 4)
        a, b = [], []
        for target in (a, b):
            n = rng.randrange(k, 7)
            target.extend([0] * n)
            for u, v in combinations(range(n), 2):
                if rng.random() < 0.35:
                    target[u] |= 1 << v
                    target[v] |= 1 << u
        la, lb = labels(len(a), k), labels(len(b), k)
        assert algebra.glued_value(algebra.state(a, la), algebra.state(b, lb)) == raw_sg(
            glue(a, la, b, lb))
        random_checks += 1
    replacements = 0
    for n in range(2, 6):
        for context in graphs(n):
            lc = labels(n, 2)
            expected = algebra.glued_value(square, algebra.state(context, lc))
            assert raw_sg(glue(SQUARE, PORTS, context, lc)) == expected
            assert raw_sg(glue(FORK, PORTS, context, lc)) == expected
            replacements += 1
    # The two distant terminal leaves are NOT a valid boundary for this rule.
    distant = (1, 0, 0, 0, 0, 2)
    assert algebra.state(SQUARE, distant) != algebra.state(FORK, distant)
    witness = None
    for n in range(2, 5):
        for context in graphs(n):
            lc = labels(n, 2)
            x, y = (raw_sg(glue(g, distant, context, lc)) for g in (SQUARE, FORK))
            if x != y:
                witness = {"context": context, "labels": lc, "values": [x, y]}
                break
        if witness:
            break
    assert witness
    # Retain all three final-column ports for the thick-board interface. The
    # two 64-entry-table counterexamples must remain distinguishable.
    h, _ = grid_graph(0x7ff4d37)
    j, _ = grid_graph(0x7ff4df4)
    lh = tuple(1 + i % 3 if i // 3 == 8 else 0
               for i in range(27) if 0x7ff4d37 >> i & 1)
    lj = tuple(1 + i % 3 if i // 3 == 8 else 0
               for i in range(27) if 0x7ff4df4 >> i & 1)
    assert algebra.state(h, lh) != algebra.state(j, lj)
    # Evaluate that actual full-column extension, rather than checking only
    # that the forms differ. Context is two full columns sharing the first.
    column, _ = grid_graph((1 << 6) - 1)
    lc = (1, 2, 3, 0, 0, 0)
    assert algebra.glued_value(algebra.state(h, lh), algebra.state(column, lc)) == 0
    assert algebra.glued_value(algebra.state(j, lj), algebra.state(column, lc)) == 4
    certificate = certificate_checks()
    audit = {"certificate": certificate, "exhaustive_products": checks,
             "random_products": random_checks, "replacement_contexts": replacements,
             "actual_reachable_examples": reachable_examples(),
             "distant_port_counterexample": witness,
             "notes13_extension_values": [0, 4], "invalid_inputs_rejected": len(invalid)}
    (ROOT / "data/proofs/kayles_boundary_audit.json").write_text(
        json.dumps(audit, indent=2) + "\n")
    print(json.dumps(audit, ensure_ascii=False))


if __name__ == "__main__":
    main()
