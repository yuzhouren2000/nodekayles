#!/usr/bin/env python3
"""Independent SG, certificate and exterior checks for neutral cores."""

from copy import deepcopy
from functools import cache
from itertools import combinations, product
import json
from pathlib import Path
import random
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from neutral_core import (neutral_core_certificate, paired_core_certificate,
                          reflection_involution, leaf_support_reduction,
                          bipartite_port_obstruction, guarded_twin_reduction)
from kayles_boundary_transfer import column_graph
from verify_neutral_core import verify


def raw_value(graph, initial=None):
    """Actual vertex-subset game; no structural theorem is used."""
    @cache
    def value(mask):
        if not mask:
            return 0
        reached = pending = mask & -mask
        while pending:
            bit = pending & -pending
            pending ^= bit
            added = graph[bit.bit_length() - 1] & mask & ~reached
            reached |= added
            pending |= added
        if reached != mask:
            return value(reached) ^ value(mask ^ reached)
        outcomes = {value(mask & ~(graph[v] | 1 << v))
                    for v in range(len(graph)) if mask >> v & 1}
        result = 0
        while result in outcomes:
            result += 1
        return result
    return value((1 << len(graph)) - 1 if initial is None else initial)


def graph_from_edges(n, edges):
    graph = [0] * n
    for u, v in edges:
        graph[u] |= 1 << v
        graph[v] |= 1 << u
    return tuple(graph)


def attach_two_path(graph, support, odd_leaves=1):
    """Add a new root adjacent to support and its odd pendant leaf group."""
    result = list(graph) + [0] * (odd_leaves + 1)
    root = len(graph)
    result[root] |= 1 << support
    result[support] |= 1 << root
    for v in range(root + 1, len(result)):
        result[root] |= 1 << v
        result[v] |= 1 << root
    return tuple(result), root


def duplicate_with_leaves(graph, twin, odd_leaves=1):
    n = len(graph)
    result = list(graph) + [graph[twin]] + [0] * odd_leaves
    for v in range(n):
        if graph[twin] >> v & 1:
            result[v] |= 1 << n
    for v in range(n + 1, len(result)):
        result[n] |= 1 << v
        result[v] |= 1 << n
    return tuple(result), n


def pair_extension(core, ports, pair_count):
    n = len(core)
    base = tuple(core) + (0,) * (2 * pair_count)
    sigma = tuple(range(n)) + tuple(n + (i ^ 1) for i in range(2 * pair_count))
    orbits = [((v, n + 2 * i), (v, n + 2 * i + 1))
              for v in ports for i in range(pair_count)]
    for i, j in combinations(range(pair_count), 2):
        a, b = n + 2 * i, n + 2 * j
        orbits += [((a, b), (a + 1, b + 1)), ((a, b + 1), (a + 1, b))]
    return base, sigma, orbits


def with_orbits(base, orbits, code):
    graph = list(base)
    for i, orbit in enumerate(orbits):
        if code >> i & 1:
            for u, v in orbit:
                graph[u] |= 1 << v
                graph[v] |= 1 << u
    return tuple(graph)


def is_automorphism(graph, sigma):
    return all(sum(1 << sigma[u] for u in range(len(graph)) if neighbours >> u & 1)
               == graph[sigma[v]] for v, neighbours in enumerate(graph))


def main():
    counts = {"small_core_port_queries": 0, "accepted_core_certificates": 0,
              "exhaustive_extensions": 0, "random_extensions": 0,
              "two_path_identities": 0, "nonautomorphic_extensions": 0,
              "bipartite_obstruction_witnesses": 0,
              "guarded_twin_identities": 0}
    examples = []
    # Every graph on <=4 vertices, every set of proposed ports. Multiple
    # direct deletions and every subsequent actual move must be closed.
    for n in range(1, 5):
        edges = tuple(combinations(range(n), 2))
        for code in range(1 << len(edges)):
            graph = graph_from_edges(n, (e for i, e in enumerate(edges) if code >> i & 1))
            for twin in range(n):
                for leaves in (1, 3):
                    enlarged, port = duplicate_with_leaves(graph, twin, leaves)
                    target = raw_value(graph) ^ 1
                    assert raw_value(enlarged) == target
                    certificate = neutral_core_certificate(enlarged, (port,))
                    assert certificate is not None and verify(certificate) == target
                    reduced, offset, removed = guarded_twin_reduction(enlarged)
                    assert raw_value(reduced) ^ offset == target
                    counts["guarded_twin_identities"] += 1
            for port in range(n):
                obstruction = bipartite_port_obstruction(graph, port)
                if obstruction is None:
                    continue
                current = (1 << n) - 1
                for v in obstruction["legal_history"]:
                    assert v != port and current >> v & 1
                    current &= ~(graph[v] | 1 << v)
                    assert current >> port & 1
                assert current == obstruction["terminal_mask"]
                target = raw_value(graph, current)
                deleted = raw_value(graph, current & ~(1 << port))
                assert target == obstruction["terminal_grundy"]
                assert deleted == obstruction["deleted_port_grundy"] != target
                assert verify({"format": "node-kayles-neutral-obstruction-v1",
                               "adjacency": graph, "witness": obstruction}) == {
                                   "terminal_grundy": target, "deleted_port_grundy": deleted}
                counts["bipartite_obstruction_witnesses"] += 1
            for port_mask in range(1 << n):
                ports = tuple(v for v in range(n) if port_mask >> v & 1)
                certificate = neutral_core_certificate(graph, ports)
                counts["small_core_port_queries"] += 1
                if certificate is None:
                    continue
                assert verify(certificate) == raw_value(graph)
                counts["accepted_core_certificates"] += 1

    cores = [((2, 1), (0,)),
             (graph_from_edges(4, ((0, 1), (1, 2), (2, 3))), (1, 2)),
             (graph_from_edges(5, ((0, 1), (0, 2), (2, 3), (2, 4))), (0,)),
             (graph_from_edges(3, combinations(range(3), 2)), (0,)),
             (graph_from_edges(6, ((0, 1), (1, 2), (0, 3), (1, 4), (2, 5))), (0, 1, 2))]
    for graph, ports in cores:
        certificate = neutral_core_certificate(graph, ports)
        assert certificate is not None
        target = verify(certificate)
        base, sigma, orbits = pair_extension(graph, ports, 2)
        for code in range(1 << len(orbits)):
            full = with_orbits(base, orbits, code)
            result = paired_core_certificate(full, sigma)
            assert result is not None and result["grundy"] == target
            assert verify(result) == target
            assert raw_value(full) == target
            counts["exhaustive_extensions"] += 1

    # Literal P4 pairing is compatible but fails to preserve individual edges.
    # Every fixed-core vertex may see both vertices of either pair.
    exterior = ((0, 1), (1, 2), (2, 3))
    for graph, ports in cores:
        n = len(graph)
        base = graph_from_edges(n + 4, [(u, v) for u in range(n)
                                      for v in range(u + 1, n) if graph[u] >> v & 1]
                                + [(u + n, v + n) for u, v in exterior])
        sigma = tuple(range(n)) + (n + 2, n + 3, n, n + 1)
        orbits = [((v, n + i), (v, n + i + 2)) for v in ports for i in range(2)]
        for code in range(1 << len(orbits)):
            full = with_orbits(base, orbits, code)
            assert not is_automorphism(full, sigma)
            result = paired_core_certificate(full, sigma)
            assert result is not None and result["grundy"] == raw_value(graph)
            assert verify(result) == raw_value(graph)
            assert raw_value(full) == raw_value(graph)
            counts["nonautomorphic_extensions"] += 1

    # The local two-point-tail identity permits arbitrary core edges, not
    # merely trees, bipartite graphs or symmetric cores.
    for n in range(1, 5):
        possible = tuple(combinations(range(n), 2))
        for code in range(1 << len(possible)):
            inner = graph_from_edges(n + 1, [(u, v) for i, (u, v) in enumerate(possible)
                                             if code >> i & 1] + [(0, n)])
            for leaves in (1, 3):
                enlarged, port = attach_two_path(inner, 0, leaves)
                certificate = neutral_core_certificate(enlarged, (port,))
                assert certificate is not None
                assert verify(certificate) == (raw_value(inner) ^ 1)
                assert raw_value(enlarged) == (raw_value(inner) ^ 1)
                reduced, offset, removed = leaf_support_reduction(enlarged)
                assert raw_value(enlarged) == raw_value(reduced) ^ offset
                counts["two_path_identities"] += 1

    # An unbounded family: the corner-defect ladder has value n. Its second
    # lower vertex has a genuine leaf at the missing-corner end.
    ladders = []
    for width in range(2, 10):
        graph, cells = column_graph((2,) + (3,) * (width - 1))
        support = cells.index((1, 1))
        assert raw_value(graph) == width
        enlarged, root = attach_two_path(graph, support)
        certificate = neutral_core_certificate(enlarged, (root,))
        assert certificate is not None and verify(certificate) == width ^ 1
        assert raw_value(enlarged) == width ^ 1
        ladders.append({"ladder_width": width, "vertices": len(enlarged),
                        "neutral_root": root, "grundy": width ^ 1,
                        "proof_states": len(certificate["states"])})

    # These fixed cores are actual reachable three-row components. Their
    # external ports still need a geometry check before embedding an extension
    # in a three-row strip; an abstract extension is not such an embedding.
    for columns, target in (((7, 3), 3), ((7, 6, 4, 6, 4, 6, 6), 2),
                            ((7, 6, 4, 6, 4, 6, 7, 7), 5),
                            ((7, 6, 4, 6, 4, 6, 7, 7, 6), 4)):
        graph, cells = column_graph(columns)
        certificate = neutral_core_certificate(graph, (1,))
        assert certificate is not None and verify(certificate) == target
        assert raw_value(graph) == target
        from strip3_history_language import decode
        mask = sum(col << (3 * c) for c, col in enumerate(columns))
        parent_width, history, component = decode(mask).witness()
        available = {(r, c) for r in range(3) for c in range(parent_width)}
        for r, c in history:
            assert (r, c) in available
            available -= {(r, c), (r - 1, c), (r + 1, c), (r, c - 1), (r, c + 1)}
        wanted = {(r, c) for r, c in available if component >> (3 * c + r) & 1}
        assert len(wanted) == len(graph)
        first_column = min(c for r, c in wanted)
        assert {(r, c - first_column) for r, c in wanted} == set(cells)
        reached = {next(iter(wanted))}
        pending = list(reached)
        while pending:
            r, c = pending.pop()
            for cell in ((r - 1, c), (r + 1, c), (r, c - 1), (r, c + 1)):
                if cell in available and cell not in reached:
                    reached.add(cell); pending.append(cell)
        assert reached == wanted
        examples.append({"columns": columns, "neutral_root": 1, "grundy": target,
                         "proof_states": len(certificate["states"]),
                         "parent_width": parent_width,
                         "legal_history": tuple((r + 1, c + 1) for r, c in history)})

    # A constant local rewrite applies to arbitrarily wide connected dense
    # three-row regions. The prefix is unrestricted; it need not be symmetric.
    tail_cases = 0
    square_cases = 0
    for n in range(1, 5):
        for prefix in product(range(8), repeat=n - 1):
            graph, _ = column_graph(prefix + (4, 6, 4, 4))
            reduced = column_graph(prefix + (4, 6))[0]
            assert raw_value(graph) == raw_value(reduced) ^ 1
            tail_cases += 1
            square, cells = column_graph(prefix + (4, 6, 6, 4))
            tree, _ = column_graph(prefix + (4, 6, 2))
            # Check the specific root used in the arbitrary-prefix proof;
            # the recognizer may prefer a different valid root in the prefix.
            k = len(prefix)
            root = cells.index((2, k + 2))
            twin = cells.index((1, k + 1))
            leaf = cells.index((2, k + 3))
            assert square[leaf] == 1 << root
            assert square[root] & ~(1 << leaf) == square[twin]
            retained = [v for v in range(len(square)) if v not in (root, leaf)]
            literal_reduction = tuple(sum(1 << j for j, u in enumerate(retained)
                                           if square[v] >> u & 1) for v in retained)
            assert literal_reduction == tree
            assert raw_value(square) == raw_value(tree) ^ 1
            reduction = guarded_twin_reduction(square)
            assert reduction is not None
            assert raw_value(reduction[0]) ^ reduction[1] == raw_value(square)
            square_cases += 1
    from guarded_symmetry import remaining_columns
    dense_cases = []
    for width in range(7, 18):
        history = ((1, width - 5), (1, width - 3), (2, width - 1))
        columns = remaining_columns(width, history)
        available = {(r, c) for r in range(1, 4) for c in range(1, width + 1)}
        for r, c in history:
            assert (r, c) in available
            available -= {(r, c), (r - 1, c), (r + 1, c), (r, c - 1), (r, c + 1)}
        assert columns == tuple(sum(1 << (r - 1) for r in range(1, 4)
                                    if (r, c) in available)
                                for c in range(1, width + 1))
        assert columns == (7,) * (width - 7) + (6, 4, 6, 4, 4, 0, 5)
        assert sum(col.bit_count() for col in columns) == 3 * width - 12
        if width <= 11:
            graph, _ = column_graph(columns)
            reduced, _ = column_graph(columns[:-4])
            assert raw_value(graph) == raw_value(reduced) ^ 1
        dense_cases.append({"width": width, "history": history,
                            "playable_cells": 3 * width - 12})

    dense_square_cases = []
    for width in range(5, 18):
        history = ((1, width - 3), (1, width))
        available = {(r, c) for r in range(1, 4) for c in range(1, width + 1)}
        for r, c in history:
            assert (r, c) in available
            available -= {(r, c), (r - 1, c), (r + 1, c), (r, c - 1), (r, c + 1)}
        columns = tuple(sum(1 << (r - 1) for r in range(1, 4) if (r, c) in available)
                        for c in range(1, width + 1))
        assert columns == (7,) * (width - 5) + (6, 4, 6, 6, 4)
        assert len(available) == 3 * width - 7
        if width <= 11:
            graph, _ = column_graph(columns)
            reduced, _ = column_graph((7,) * (width - 5) + (6, 4, 6, 2))
            assert raw_value(graph) == raw_value(reduced) ^ 1
        dense_square_cases.append({"width": width, "history": history,
                                   "playable_cells": 3 * width - 7})

    rng = random.Random(31020)
    for _ in range(200):
        k = rng.randrange(1, 7)
        graph = graph_from_edges(k + 1, [(u, v) for u, v in combinations(range(k), 2)
                                        if rng.randrange(3) == 0] + [(0, k)])
        graph, port = attach_two_path(graph, 0)
        base, sigma, orbits = pair_extension(graph, (port,), rng.randrange(4))
        full = with_orbits(base, orbits, rng.randrange(1 << len(orbits)))
        result = paired_core_certificate(full, sigma)
        assert result is not None and result["grundy"] == raw_value(full)
        assert verify(result) == raw_value(full)
        counts["random_extensions"] += 1

    # A one-time equality is insufficient: P5 rooted at its second vertex
    # changes from 3 to 4 when two identical two-leaf stars are attached.
    one_time = (6, 1, 9, 20, 8)
    assert raw_value(one_time) == raw_value(one_time, 30) == 3
    assert neutral_core_certificate(one_time, (0,)) is None
    failure = (102, 1, 9, 20, 8, 641, 1281, 32, 64, 32, 64)
    sigma = tuple(range(5)) + (6, 5, 8, 7, 10, 9)
    assert is_automorphism(failure, sigma) and raw_value(failure) == 4
    assert paired_core_certificate(failure, sigma) is None
    # Omitting the existing support leaf breaks the two-point-tail formula.
    no_support_leaf = (0,)
    enlarged, root = attach_two_path(no_support_leaf, 0)
    assert raw_value(enlarged) == 2 != (raw_value(no_support_leaf) ^ 1)
    assert neutral_core_certificate(enlarged, (root,)) is None
    assert leaf_support_reduction(enlarged) is None
    # Odd original leaf count is necessary, not sufficient: P5's second
    # vertex passes this linear obstruction but fails hereditary neutrality.
    assert bipartite_port_obstruction(one_time, 0) is None
    full_strip_obstructions = []
    for width in range(1, 25):
        graph, _ = column_graph((7,) * width)
        for port in range(len(graph)):
            obstruction = bipartite_port_obstruction(graph, port)
            assert obstruction is not None
            current = (1 << len(graph)) - 1
            for v in obstruction["legal_history"]:
                assert current >> v & 1
                current &= ~(graph[v] | 1 << v)
            assert current >> port & 1 and current == obstruction["terminal_mask"]
            verify({"format": "node-kayles-neutral-obstruction-v1",
                    "adjacency": graph, "witness": obstruction})
        full_strip_obstructions.append({"width": width, "rejected_ports": len(graph)})
    obstruction_graph, _ = column_graph((7,) * 12)
    obstruction_proof = {"format": "node-kayles-neutral-obstruction-v1",
                         "adjacency": obstruction_graph,
                         "witness": bipartite_port_obstruction(obstruction_graph, 16)}
    obstruction_tampered = []
    bad = deepcopy(obstruction_proof); bad["witness"]["legal_history"] = (16,); obstruction_tampered.append(bad)
    bad = deepcopy(obstruction_proof); bad["witness"]["terminal_grundy"] ^= 1; obstruction_tampered.append(bad)
    bad = deepcopy(obstruction_proof); bad["witness"]["terminal_mask"] ^= 1; obstruction_tampered.append(bad)
    for bad in obstruction_tampered:
        try:
            verify(bad)
        except ValueError:
            continue
        raise AssertionError("tampered obstruction accepted")
    # Removing all neutral ports may change a proposed value on an ordinary
    # full odd strip. Recognition must reject it, not infer an empty-board SG.
    for axis, widths in (("horizontal", (3, 5, 7)), ("vertical", (1, 2, 3))):
        for width in widths:
            reflected = reflection_involution((7,) * width, axis)
            assert reflected is not None and paired_core_certificate(*reflected) is None

    # Independent checker catches SG changes, missing descendants, forged
    # root values, and duplicate proof rows.
    valid = neutral_core_certificate((6, 1, 25, 4, 4), (0,))
    tampered = []
    bad = deepcopy(valid); bad["grundy"] ^= 1; tampered.append(bad)
    bad = deepcopy(valid); bad["states"] = bad["states"][1:]; tampered.append(bad)
    bad = deepcopy(valid); bad["states"] += (bad["states"][0],); tampered.append(bad)
    bad = deepcopy(valid); bad["states"] = tuple((m, g ^ 1) if m == 31 else (m, g)
                                                for m, g in bad["states"]); tampered.append(bad)
    bad = deepcopy(valid); bad["ports"] = (2,); tampered.append(bad)
    for bad in tampered:
        try:
            verify(bad)
        except ValueError:
            continue
        raise AssertionError("tampered proof accepted")

    paired = paired_core_certificate((14, 1, 1, 1), (0, 1, 3, 2))
    assert paired is not None and verify(paired) == 1
    paired_tampered = []
    bad = deepcopy(paired); bad["involution"] = (0, 1, 2, 3); paired_tampered.append(bad)
    bad = deepcopy(paired); bad["core_ports"] = (); paired_tampered.append(bad)
    bad = deepcopy(paired); bad["core"]["ports"] = (); paired_tampered.append(bad)
    bad = deepcopy(paired); bad["grundy"] ^= 1; paired_tampered.append(bad)
    bad = deepcopy(paired); bad["adjacency"] = (6, 5, 3, 0); paired_tampered.append(bad)
    for bad in paired_tampered:
        try:
            verify(bad)
        except ValueError:
            continue
        raise AssertionError("tampered paired proof accepted")

    audit = {"format": "node-kayles-neutral-core-audit-v1", **counts,
             "unbounded_family": ladders, "reachable_three_row_cores": examples,
             "arbitrary_three_row_tail_rewrites": tail_cases,
             "arbitrary_three_row_square_rewrites": square_cases,
             "dense_three_move_histories": dense_cases,
             "dense_two_move_square_histories": dense_square_cases,
             "one_time_equality_counterexample": {"core": one_time, "port": 0,
                  "core_value": 3, "extended_value": 4, "extended_graph": failure},
             "empty_strip_rejections": 6, "tamper_rejections": len(tampered),
             "paired_tamper_rejections": len(paired_tampered),
             "obstruction_tamper_rejections": len(obstruction_tampered),
             "full_strip_neutral_port_obstructions": full_strip_obstructions}
    (ROOT / "data/proofs/neutral_core_audit.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2) + "\n")
    (ROOT / "data/proofs/neutral_core_sg5.json").write_text(json.dumps(
        neutral_core_certificate(column_graph((7, 6, 4, 6, 4, 6, 7, 7))[0], (1,)),
        indent=2) + "\n")
    paired_example = paired_core_certificate(*reflection_involution((7, 7, 6, 4, 6, 4, 6, 7, 7)))
    assert paired_example is not None and verify(paired_example) == 1
    (ROOT / "data/proofs/paired_neutral_core_9.json").write_text(
        json.dumps(paired_example, indent=2) + "\n")
    (ROOT / "data/proofs/bipartite_neutral_obstruction_12.json").write_text(
        json.dumps(obstruction_proof, indent=2) + "\n")
    print(json.dumps(audit, ensure_ascii=False))


if __name__ == "__main__":
    main()
