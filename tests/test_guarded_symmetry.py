#!/usr/bin/env python3
"""Independent vertex-set SG and legal-history checks of guarded symmetry."""

from functools import cache
from itertools import combinations, product
import json
from pathlib import Path
import random
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from guarded_symmetry import (central_pair_history, four_move_p_history,
                              guarded_certificate, reflection_certificate,
                              remaining_columns)


def raw_value(graph):
    """All moves on actual vertex subsets, without any symmetry reductions."""
    @cache
    def value(mask):
        if not mask:
            return 0
        # Independent connected-component splitting; no theorem recognition.
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
    return value((1 << len(graph)) - 1)


def induced(graph, sigma, removed):
    keep = [v for v in range(len(graph)) if not removed >> v & 1]
    index = {v: i for i, v in enumerate(keep)}
    child = tuple(sum(1 << index[u] for u in keep if graph[v] >> u & 1)
                  for v in keep)
    return child, tuple(index[sigma[v]] for v in keep)


def guarded_graph(leaf_counts, pair_count, isolates=0):
    """All possible invariant edges, including edges across swapped halves."""
    k = len(leaf_counts)
    graph = [0] * (k + sum(leaf_counts) + 2 * pair_count + isolates)
    sigma = list(range(len(graph)))
    cursor = k
    for centre, count in enumerate(leaf_counts):
        for _ in range(count):
            graph[centre] |= 1 << cursor
            graph[cursor] |= 1 << centre
            cursor += 1
    pairs = []
    for _ in range(pair_count):
        pairs.append((cursor, cursor + 1))
        sigma[cursor], sigma[cursor + 1] = cursor + 1, cursor
        cursor += 2
    edge_orbits = [((u, v),) for u, v in combinations(range(k), 2)]
    edge_orbits += [((v, a), (v, b)) for v in range(k) for a, b in pairs]
    for (a, b), (c, d) in combinations(pairs, 2):
        edge_orbits += [((a, c), (b, d)), ((a, d), (b, c))]
    return graph, tuple(sigma), edge_orbits


def add_edges(base, orbits, code):
    graph = list(base)
    for i, orbit in enumerate(orbits):
        if code >> i & 1:
            for u, v in orbit:
                graph[u] |= 1 << v
                graph[v] |= 1 << u
    return tuple(graph)


def coordinate_graph(cells):
    cells = tuple(sorted(cells))
    return tuple(sum(1 << j for j, u in enumerate(cells)
                     if abs(v[0] - u[0]) + abs(v[1] - u[1]) == 1)
                 for v in cells)


def actual_history(width, history):
    """Use coordinate sets, independently of bitwise column deletion."""
    available = {(r, c) for r in (1, 2, 3) for c in range(1, width + 1)}
    for r, c in history:
        assert (r, c) in available
        available -= {(r, c), (r - 1, c), (r + 1, c), (r, c - 1), (r, c + 1)}
    columns = remaining_columns(width, history)
    assert available == {(r, c) for c, col in enumerate(columns, 1)
                         for r in (1, 2, 3) if col >> (r - 1) & 1}
    return columns, coordinate_graph(available)


def main():
    counts = {"exhaustive_graphs": 0, "move_closure_checks": 0,
              "random_graphs": 0, "central_histories": 0,
              "symmetric_history_variants": 0, "four_move_histories": 0}
    # No assumption that paired vertices live in separate components. Both
    # parallel and crossing edge orbits between pairs are independently free.
    configs = (((), 4), ((1,), 3), ((1, 1), 2), ((1, 1, 1), 1),
               ((3,), 2), ((1, 3), 1), ((3, 3), 1))
    for leaves, pairs in configs:
        base, sigma, orbits = guarded_graph(leaves, pairs)
        for code in range(1 << len(orbits)):
            graph = add_edges(base, orbits, code)
            cert = guarded_certificate(graph, sigma)
            assert cert is not None and cert.grundy == len(leaves) % 2
            assert raw_value(graph) == cert.grundy
            for v in range(len(graph)):
                removed = graph[v] | 1 << v
                if sigma[v] == v:
                    child, child_sigma = induced(graph, sigma, removed)
                    target = cert.grundy ^ 1
                else:
                    partner = sigma[v]
                    assert not removed >> partner & 1
                    removed |= graph[partner] | 1 << partner
                    child, child_sigma = induced(graph, sigma, removed)
                    target = cert.grundy
                successor = guarded_certificate(child, child_sigma)
                assert successor is not None and successor.grundy == target
                counts["move_closure_checks"] += 1
            counts["exhaustive_graphs"] += 1

    rng = random.Random(7320)
    for _ in range(200):
        k, pairs, isolates = rng.randrange(1, 4), rng.randrange(4), rng.randrange(3)
        leaves = [1] * k
        if rng.randrange(2):
            leaves[rng.randrange(k)] = 3
        base, sigma, orbits = guarded_graph(leaves, pairs, isolates)
        graph = add_edges(base, orbits, rng.randrange(1 << len(orbits)))
        cert = guarded_certificate(graph, sigma)
        assert cert is not None and cert.grundy == (k + isolates) % 2
        assert raw_value(graph) == cert.grundy
        counts["random_graphs"] += 1

    # Every small invariant induced graph with an explicitly guarded centre
    # column, even when its two sides contain arbitrary holes and components.
    masks = 0
    for half in range(4):
        for prefix in product(range(8), repeat=half):
            columns = prefix + (6,) + prefix[::-1]
            if half and prefix[-1] & 2:
                continue  # the middle fixed vertex must be a pendant leaf
            cert = reflection_certificate(columns)
            assert cert is not None and cert.grundy == 1
            cells = {(r, c) for c, col in enumerate(columns)
                     for r in range(3) if col >> r & 1}
            assert raw_value(coordinate_graph(cells)) == 1
            masks += 1
    counts["symmetric_masks"] = masks

    for width in range(3, 18, 2):
        columns, graph = actual_history(width, central_pair_history(width))
        assert reflection_certificate(columns).grundy == 1
        # Raw recursion is bounded here; larger cases use the theorem only.
        if width <= 13:
            assert raw_value(graph) == 1
        expected_cells = 3 * width - (5 if width == 3 else 7)
        assert sum(col.bit_count() for col in columns) == expected_cells
        centre = (width + 1) // 2
        if width <= 13:
            for move in ((2, centre), (3, centre)):
                _, child = actual_history(width, central_pair_history(width) + (move,))
                assert raw_value(child) == 0
        counts["central_histories"] += 1

    # Further legal history may be arbitrary in each half, then reflected.
    for half in range(1, 6):
        width = 2 * half + 1
        for choices in product((0, 1, 3), repeat=half - 1):
            side = choices + (1,)
            word = side + (0,) + side[::-1]
            history = tuple((row, c) for c, row in enumerate(word, 1) if row)
            try:
                columns, graph = actual_history(width, history)
            except (AssertionError, ValueError):
                continue
            cert = reflection_certificate(columns)
            assert cert is not None and cert.grundy == 1
            assert raw_value(graph) == 1
            counts["symmetric_history_variants"] += 1

    small_core = coordinate_graph({(r, c) for r in range(1, 4) for c in range(1, 4)
                                   if (r, c) != (2, 3)})
    core_options = [raw_value(induced(small_core, tuple(range(8)), ns | 1 << v)[0])
                    for v, ns in enumerate(small_core)]
    assert core_options == [0] * 8 and raw_value(small_core) == 1

    family_examples = []
    for width in range(8, 22):
        history = four_move_p_history(width)
        columns, graph = actual_history(width, history)
        even = width - width % 2
        assert columns[:4] == (7, 7, 5, 0) and columns[even - 1] == 0
        assert columns[even:] == ((5,) if width % 2 else ())
        assert reflection_certificate(columns[4:even - 1]).grundy == 1
        if width <= 19:
            assert raw_value(graph) == 0
        assert sum(col.bit_count() for col in columns) == 3 * width - (12 if width < 10 else 16) - width % 2
        family_examples.append({"width": width, "history": history, "grundy": 0})
        counts["four_move_histories"] += 1

    # Each omitted condition has a concrete counterexample.
    no_guard = (6, 1, 1)  # P3, reflection fixing its centre
    assert raw_value(no_guard) == 2
    assert guarded_certificate(no_guard, (0, 2, 1)) is None
    even_leaves, sigma, _ = guarded_graph((2,), 0)
    assert guarded_certificate(even_leaves, sigma) is None and raw_value(even_leaves) == 2
    adjacent_pair = (2, 1, 8, 4)  # fixed K2 plus a swapped adjacent K2
    assert guarded_certificate(adjacent_pair, (0, 1, 3, 2)) is None
    assert raw_value(adjacent_pair) == 0
    not_pendant = (14, 13, 3, 3)  # fixed edge, both endpoints touch the pair
    assert guarded_certificate(not_pendant, (0, 1, 3, 2)) is None
    assert raw_value(not_pendant) == 2

    invalid = (((2, 0), (0, 1)), ((0,), (1,)), ((2, 1), (1, 1)),
               ((2, 1), (True, 0)), ((2, 5, 2), (1, 0, 2)))
    for graph, sigma in invalid:
        try:
            guarded_certificate(graph, sigma)
        except ValueError:
            continue
        raise AssertionError("invalid symmetry accepted")
    for function, widths in ((central_pair_history, (0, 1, 4, True)),
                             (four_move_p_history, (0, 6, 7, True))):
        for width in widths:
            try:
                function(width)
            except ValueError:
                continue
            raise AssertionError("invalid family width accepted")

    assert central_pair_history(1_000_001) == ((1, 500_000), (1, 500_002))
    assert four_move_p_history(1_000_000_000) == (
        (2, 4), (1, 500_000_001), (1, 500_000_003), (2, 1_000_000_000))
    audit = {"format": "node-kayles-guarded-symmetry-audit-v1", **counts,
             "small_core": {"columns": [7, 7, 5], "option_values": core_options,
                            "grundy": 1},
             "four_move_p_examples": family_examples,
             "counterexamples": 4, "invalid_input_checks": len(invalid) + 8}
    (ROOT / "data/proofs/guarded_symmetry_audit.json").write_text(
        json.dumps(audit, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(audit, ensure_ascii=False))


if __name__ == "__main__":
    main()
