#!/usr/bin/env python3
"""Check interface operations against independent vertex-set graph recursion."""

from itertools import combinations, product
import json
from pathlib import Path
import random
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from kayles_boundary_algebra import FORK, SQUARE
from kayles_boundary_transfer import BoundaryTransfer, ThreeRowTransfer
from test_kayles_boundary import glue, graphs, grid_graph, raw_sg


def raw(graph):
    return raw_sg(glue(graph, (0,) * len(graph), (), ()))


def random_graph(rng, n):
    graph = [0] * n
    for u, v in combinations(range(n), 2):
        if rng.random() < 0.4:
            graph[u] |= 1 << v
            graph[v] |= 1 << u
    return tuple(graph)


def random_labels(rng, n):
    ports = rng.sample((1, 2, 3, 4), rng.randrange(min(n, 4) + 1))
    result = ports + [0] * (n - len(ports))
    rng.shuffle(result)
    return tuple(result)


def named_join(*items):
    """Independently identify ports in several actual named graphs."""
    edges = {}
    for side, (graph, labels) in enumerate(items):
        names = [f"p{label}" if label else f"s{side}v{i}"
                 for i, label in enumerate(labels)]
        for name in names:
            edges.setdefault(name, set())
        for v, neighbours in enumerate(graph):
            edges[names[v]].update(names[u] for u in range(len(graph))
                                   if neighbours >> u & 1)
    return tuple((v, tuple(sorted(edges[v]))) for v in sorted(edges))


def graph_for_word(word):
    mask = sum(w << (3 * c) for c, w in enumerate(word))
    graph, cells = grid_graph(mask)
    labels = tuple(i % 3 + 1 if i // 3 == len(word) - 1 else 0 for i in cells)
    return graph, labels


def main():
    rng = random.Random(7319)
    algebra = BoundaryTransfer()
    # Exhaust partial-overlap interfaces, including no shared ports and
    # disconnected graphs. Port 2 occurs on either side independently.
    samples = [(graph, boundary) for n in range(1, 4) for graph in graphs(n)
               for boundary in product((0, 1, 2), repeat=n)
               if len([i for i in boundary if i]) == len(set(boundary) - {0})]
    joins = 0
    for a, la in samples:
        sa = algebra.state(a, la)
        for b, lb in samples:
            joined = algebra.join(sa, algebra.state(b, lb))
            assert algebra.value(joined) == raw_sg(glue(a, la, b, lb))
            assert algebra.join(algebra.state(b, lb), sa) == joined
            joins += 1

    contextual = 0
    forgets = 0
    renames = 0
    associates = 0
    for _ in range(800):
        items = [(random_graph(rng, n), random_labels(rng, n))
                 for n in (rng.randrange(1, 6), rng.randrange(1, 6), rng.randrange(1, 5))]
        states = [algebra.state(g, labels) for g, labels in items]
        joined = algebra.join(states[0], states[1])
        context = states[2]
        assert algebra.value(algebra.join(joined, context)) == raw_sg(named_join(*items))
        contextual += 1
        # Associativity is checked semantically in a third graph context.
        other_order = algebra.join(states[0], algebra.join(states[1], context))
        assert algebra.value(other_order) == algebra.value(algebra.join(joined, context))
        associates += 1

        a, la = items[0]
        removed = tuple(i for i in (1, 2, 3, 4) if rng.randrange(2))
        lf = tuple(0 if i in removed else i for i in la)
        forgotten = algebra.forget(states[0], removed)
        b, lb = items[2]
        # A smaller-boundary context must not reconnect the forgotten ports.
        smaller = tuple(0 if i in removed else i for i in lb)
        assert algebra.value(algebra.join(forgotten, algebra.state(b, smaller))) == raw_sg(
            glue(a, lf, b, smaller))
        assert algebra.value(forgotten) == raw(a)
        forgets += 1

        targets = rng.sample((11, 12, 13, 14), 4)
        mapping = dict(zip((1, 2, 3, 4), targets))
        renamed = algebra.rename(joined, mapping)
        new_items = [(g, tuple(mapping[i] if i else 0 for i in labels)) for g, labels in items]
        assert algebra.value(algebra.join(renamed, algebra.state(*new_items[2]))) == raw_sg(
            named_join(*new_items))
        renames += 1

    # Detached nimbers must survive every operation, including empty states.
    st = algebra.state((2, 1), (1, 2))
    assert algebra.value(algebra.join(st, (0, 9))) == algebra.value(st) ^ 9
    assert algebra.forget((st[0], 9), (1, 2)) == (0, algebra.value(st) ^ 9)
    assert algebra.value(algebra.rename((st[0], 9), {1: 3, 2: 4})) == algebra.value(st) ^ 9
    assert algebra.state((), ()) == (0, 0)

    model = ThreeRowTransfer()
    word_checks = 0
    word_contexts = 0
    for n in range(5):
        for word in product(range(8), repeat=n):
            state = model.state(word)
            graph, boundary = graph_for_word(word)
            assert model.algebra.value(state) == raw(graph)
            word_checks += 1
            if n <= 3:
                m = rng.randrange(1, 6)
                context = random_graph(rng, m)
                labels = random_labels(rng, m)
                actual = raw_sg(glue(graph, boundary, context, labels))
                combined = model.algebra.join(state, model.algebra.state(context, labels))
                assert model.algebra.value(combined) == actual
                word_contexts += 1
    for _ in range(400):
        word = tuple(rng.randrange(8) for _ in range(rng.randrange(5, 8)))
        graph, _ = graph_for_word(word)
        assert model.value(word) == raw(graph)
        word_checks += 1

    # Historical independent C results, with a separate raw recursion through 8.
    baseline = {}
    for line in (ROOT / "data/raw/strip3_3xn.tsv").read_text().splitlines():
        if line and not line.startswith("#"):
            m, n, g = map(int, line.split())
            if m == 3:
                baseline[n] = g
    full = ThreeRowTransfer()
    state, previous = (0, 0), 0
    values = []
    counts = []
    for n in range(1, 13):
        state = full.append(state, previous, 7)
        previous = 7
        value = full.algebra.value(state)
        assert value == baseline[n]
        if n <= 8:
            graph, _ = graph_for_word((7,) * n)
            assert value == raw(graph)
        values.append(value)
        counts.append(len(full.algebra.nodes))

    obstruction = []
    for mask in (0x7ff4d37, 0x7ff4df4):
        word = tuple(mask >> (3 * i) & 7 for i in range(9))
        value = model.value(word + (7,))
        graph, _ = graph_for_word(word + (7,))
        assert value == raw(graph)
        obstruction.append(value)
    assert obstruction == [0, 4]

    # Premises of the infinite-interface theorem, checked on actual graphs.
    # Corner-defect ladders carry n, and a full ladder context adds length
    # by ordinary addition, not xor with a context-independent nimber.
    ladders = ThreeRowTransfer()
    ladder_extensions = 0
    for n in range(2, 9):
        state = ladders.state((2,) + (3,) * (n - 1))
        left_graph, left_labels = graph_for_word((2,) + (3,) * (n - 1))
        for k in range(4):
            mask = sum(3 << (3 * c) for c in range(k + 1))
            ctx_graph, cells = grid_graph(mask)
            ctx_labels = tuple(i % 3 + 1 if i // 3 == 0 else 0 for i in cells)
            actual = raw_sg(glue(left_graph, left_labels, ctx_graph, ctx_labels))
            assert actual == n + k
            assert ladders.algebra.value(ladders.algebra.join(
                state, ladders.algebra.state(ctx_graph, ctx_labels))) == actual
            ladder_extensions += 1
    initial, extended = [2, 3], [3, 4]
    assert initial[0] ^ initial[1] == 1
    assert extended[0] ^ extended[1] == 7

    # A genuine failure of the distant two-port square/tree replacement, even
    # with separate path contexts on the two sides, rather than a closing edge.
    distant = (1, 0, 0, 0, 0, 2)
    left = ((2, 1), (1, 0))
    right = ((2, 1), (2, 0))
    end_values = [raw_sg(named_join((g, distant), left, right)) for g in (SQUARE, FORK)]
    assert end_values == [0, 1]

    # Three copies of an arbitrary rooted branch cannot be reduced to one.
    branch = ((2, 5, 2), (1, 0, 0))
    context = ((6, 25, 25, 6, 6), (1, 0, 0, 0, 0))  # K2,3
    bouquet_values = [raw_sg(named_join(context, *([branch] * k))) for k in (1, 3)]
    assert bouquet_values == [4, 2]

    invalid = [lambda: algebra.join((-1, 0), (0, 0)),
               lambda: algebra.join((0, -1), (0, 0)),
               lambda: algebra.forget(st, (0,)),
               lambda: algebra.rename(st, {1: 3}),
               lambda: algebra.rename(st, {1: 3, 2: 3}),
               lambda: algebra.rename(st, ((1, 3), (1, 4), (2, 5))),
               lambda: algebra.rename(st, {1: 0, 2: 4}),
               lambda: model.state((True,)),
               lambda: model.state((8,)),
               lambda: model.append((0, 0), 7, 7)]
    for operation in invalid:
        try:
            operation()
        except ValueError:
            continue
        raise AssertionError("invalid operation accepted")

    audit = {"exhaustive_partial_joins": joins, "contextual_joins": contextual,
             "forget_contexts": forgets, "rename_contexts": renames,
             "associativity_checks": associates, "strip_words": word_checks,
             "strip_contexts": word_contexts, "full_widths": list(range(1, 13)),
             "full_values": values, "full_accumulated_nodes": counts,
             "raw_full_widths": list(range(1, 9)), "notes13_extension_values": obstruction,
             "raw_corner_ladder_extensions": ladder_extensions,
             "nimber_translation_counterexample": {"widths": [2, 3],
                 "initial_values": initial, "extension_columns": 1,
                 "extended_values": extended, "initial_xor": 1, "extended_xor": 7},
             "distant_separate_path_counterexample": {"values": end_values,
                 "left_context": left, "right_context": right},
             "arbitrary_branch_parity_counterexample": {"values": bouquet_values,
                 "branch": branch, "context": context, "copy_counts": [1, 3]},
             "invalid_inputs_rejected": len(invalid)}
    (ROOT / "data/proofs/kayles_boundary_transfer_audit.json").write_text(
        json.dumps(audit, indent=2) + "\n")
    print(json.dumps(audit, ensure_ascii=False))


if __name__ == "__main__":
    main()
