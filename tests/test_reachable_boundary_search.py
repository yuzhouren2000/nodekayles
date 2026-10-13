#!/usr/bin/env python3
"""Independent checks of reachability witnesses and isomorphism filtering."""

from itertools import combinations, permutations
from pathlib import Path
import random
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from search_reachable_boundary import (boundary_graph, labelled_isomorphic,
                                      reachable_columns)
from strip3_history_language import decode


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
    print({"reachable_witnesses": checked, "isomorphism_checks": isomorphisms})


if __name__ == "__main__":
    main()
