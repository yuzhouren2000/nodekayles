#!/usr/bin/env python3
"""Independent checks of reachability witnesses and isomorphism filtering."""

from itertools import combinations, permutations
import json
from pathlib import Path
import random
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from search_reachable_boundary import (boundary_graph, labelled_isomorphic,
                                      reachable_columns)
from strip3_history_language import decode
from kayles_boundary_transfer import column_graph
from test_kayles_boundary import glue, raw_sg
from strip3_defect_algebra import full_head
from verify_comb_bridge import original_game


def brute_isomorphic(a, la, b, lb):
    if len(a) != len(b):
        return False
    return any(all(la[v] == lb[p[v]] for v in range(len(a)))
               and all(bool(a[v] >> u & 1) == bool(b[p[v]] >> p[u] & 1)
                       for v in range(len(a)) for u in range(len(a)))
               for p in permutations(range(len(a))))


def main():
    checked = 0
    for width in range(2, 9):
        for columns in reachable_columns(width):
            mask = sum(col << (3 * c) for c, col in enumerate(columns))
            parent_width, history, wanted_mask = decode(mask).witness()
            available = {(r, c) for c in range(parent_width) for r in range(3)}
            for r, c in history:
                assert (r, c) in available
                available -= {(r, c), (r - 1, c), (r + 1, c), (r, c - 1), (r, c + 1)}
            wanted = {(r, c) for r, c in available if wanted_mask >> (3 * c + r) & 1}
            assert len(wanted) == mask.bit_count()
            reached, pending = {next(iter(wanted))}, [next(iter(wanted))]
            while pending:
                r, c = pending.pop()
                for cell in ((r - 1, c), (r + 1, c), (r, c - 1), (r, c + 1)):
                    if cell in available and cell not in reached:
                        reached.add(cell)
                        pending.append(cell)
            assert reached == wanted
            for boundary in ("right", "both"):
                graph, labels = boundary_graph(columns, boundary)
                for label in labels:
                    assert labels.count(label) == 1 or not label
            checked += 1
    rng = random.Random(22104)
    isomorphisms = 0
    for n in range(1, 6):
        for _ in range(100):
            graphs, labels = [], []
            for side in range(2):
                graph = [0] * n
                for u, v in combinations(range(n), 2):
                    if rng.randrange(2):
                        graph[u] |= 1 << v
                        graph[v] |= 1 << u
                graphs.append(tuple(graph))
                labels.append((1,) + (0,) * (n - 1))
            assert labelled_isomorphic(graphs[0], labels[0], graphs[1], labels[1]) == (
                brute_isomorphic(graphs[0], labels[0], graphs[1], labels[1]))
            p = list(range(n))
            rng.shuffle(p)
            renamed = [0] * n
            renamed_labels = [0] * n
            for v in range(n):
                renamed[p[v]] = sum(1 << p[u] for u in range(n) if graphs[0][v] >> u & 1)
                renamed_labels[p[v]] = labels[0][v]
            assert labelled_isomorphic(graphs[0], labels[0], renamed, renamed_labels)
            isomorphisms += 2
    # Even a three-number ROOT summary loses genuinely relevant internal moves.
    context, _ = column_graph((7, 7, 2))
    context_labels = (0,) * (len(context) - 1) + (1,)
    triples, extended, items = [], [], []
    edge_values = []
    thick_core = full_head(7)[1]
    for columns in ((4,), (4, 6, 4, 6, 4)):
        graph, _ = column_graph(columns)
        labels = (0,) * (len(graph) - 1) + (1,)
        named = glue(graph, labels, (), ())
        edges = dict(named)
        def removed(vertices):
            return tuple((v, tuple(u for u in ns if u not in vertices))
                         for v, ns in named if v not in vertices)
        triple = (raw_sg(named), raw_sg(removed({"port1"})),
                  raw_sg(removed({"port1"} | set(edges["port1"]))))
        target = raw_sg(glue(graph, labels, context, context_labels))
        ordering = [len(graph) - 1] + list(range(len(graph) - 1))
        rooted_graph = tuple(sum(1 << j for j, u in enumerate(ordering)
                                 if graph[v] >> u & 1) for v in ordering)
        edge_values.append(original_game(rooted_graph, thick_core, 0, 0)[0])
        triples.append(triple)
        extended.append(target)
        items.append({"columns": columns, "adjacency": graph, "labels": labels,
                      "root_triple": triple, "glued_value": target})
    assert triples == [(1, 0, 0), (1, 0, 0)] and extended == [2, 4]
    assert edge_values == [3, 6]
    reachable_joined_words = [(7,) * 6 + (6, 4),
                              (7,) * 6 + (6, 4, 6, 4, 6, 4)]
    for columns, expected in zip(reachable_joined_words, edge_values):
        width = len(columns)
        available = {(r, c) for c in range(width) for r in range(3)}
        for c in range(7, width, 2):
            assert (0, c) in available
            available -= {(0, c), (1, c), (0, c - 1), (0, c + 1)}
        actual = tuple(sum(1 << r for r in range(3) if (r, c) in available)
                       for c in range(width))
        assert actual == columns
        graph, _ = column_graph(columns)
        assert raw_sg(glue(graph, (0,) * len(graph), (), ())) == expected
    proof = {"format": "root-triple-obstruction-v1", "graphs": items,
             "context": {"columns": (7, 7, 2), "adjacency": context,
                         "labels": context_labels},
             "reachable_edge_context": {"head_width": 7, "adjacency": thick_core,
                                         "root": 0, "joined_values": edge_values,
                                         "joined_columns": reachable_joined_words}}
    (ROOT / "data/proofs/root_triple_obstruction.json").write_text(
        json.dumps(proof, indent=2) + "\n")
    print({"reachable_witnesses": checked, "isomorphism_checks": isomorphisms,
           "root_identification_counterexample": extended,
           "reachable_edge_root_summary_counterexample": edge_values})


if __name__ == "__main__":
    main()
