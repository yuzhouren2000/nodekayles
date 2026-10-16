"""Coordinate-only checks for the tree classification and interval evaluator.

No C imports, game-search cache, or precomputed SG table. The finite local
check covers left/right distances 0..3 from the marked column; the proof that
this covers arbitrary strip widths is in notes/38_caterpillar_components.md.
"""
from collections import Counter
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import time


def adjacency(mask):
    points = {(v//3, v%3) for v in range(mask.bit_length()) if mask >> v & 1}
    return {v: {u for u in ((v[0]-1, v[1]), (v[0]+1, v[1]),
                          (v[0], v[1]-1), (v[0], v[1]+1)) if u in points} for v in points}


def caterpillar_spine(graph):
    """Return attached leaf counts, or None unless this is a caterpillar."""
    if not graph:
        return None
    reached, pending = set(), [next(iter(graph))]
    while pending:
        v = pending.pop()
        if v not in reached:
            reached.add(v); pending.extend(graph[v]-reached)
    if len(reached) != len(graph) or sum(map(len, graph.values())) != 2*(len(graph)-1):
        return None
    spine = {v for v in graph if len(graph[v]) >= 2}
    if not spine:
        return []  # K1 or K2, SG=1
    if any(len(graph[v] & spine) > 2 for v in spine):
        return None
    start = min(v for v in spine if len(graph[v] & spine) < 2)
    result, previous, v = [], None, start
    while v is not None:
        result.append(len(graph[v]-spine))
        following = (graph[v] & spine) - {previous}
        previous, v = v, next(iter(following)) if following else None
    assert len(result) == len(spine)
    return result


def interval_value(leaves):
    if leaves == []:
        return 1
    p = len(leaves)
    values = {}

    def part(a, b):
        return values[a, b] if a <= b else 0

    for length in range(1, p+1):
        for a in range(p-length+1):
            b = a+length-1
            options = set()
            for i in range(a, b+1):
                options.add(part(a, i-2) ^ part(i+2, b)
                            ^ (leaves[i-1] % 2 if i > a else 0)
                            ^ (leaves[i+1] % 2 if i < b else 0))
                if leaves[i]:
                    options.add(part(a, i-1) ^ part(i+1, b) ^ ((leaves[i]-1) % 2))
            g = 0
            while g in options:
                g += 1
            values[a, b] = g
    return values[0, p-1]


def local_check():
    counts = Counter()
    for left in range(4):
        for right in range(4):
            n, c = left+right+1, left
            all_points = {(x, y) for x in range(n) for y in range(3)}
            neighbours = {v: {u for u in ((v[0]-1, v[1]), (v[0]+1, v[1]),
                                         (v[0], v[1]-1), (v[0], v[1]+1)) if u in all_points}
                          for v in all_points}

            def histories(prefix=()):
                if len(prefix) == n:
                    yield {(x, y) for x, a in enumerate(prefix) for y in range(3) if a >> y & 1}
                else:
                    for a in (0, 1, 2, 4, 5):
                        if not prefix or not a & prefix[-1]:
                            yield from histories((*prefix, a))

            for history in histories():
                counts["local_histories"] += 1
                deleted = set(history)
                for v in history:
                    deleted.update(neighbours[v])
                live = all_points-deleted
                for r in range(3):
                    v = c, r
                    if v not in live:
                        continue
                    counts["live_centres"] += 1
                    nonleaf = [u for u in neighbours[v] & live if len(neighbours[u] & live) >= 2]
                    if len(nonleaf) >= 3:
                        counts["three_nonleaf_centres"] += 1
                        assert any({(x, y), (x+1, y), (x, y+1), (x+1, y+1)} <= live
                                   for x in (c-1, c) for y in (r-1, r)), (left, right, history, v)
                # Independent global corollary: every bounded face is a unit
                # square. For a grid graph this says beta equals square count.
                graph = {v: neighbours[v] & live for v in live}
                unseen, components = set(live), 0
                while unseen:
                    components += 1
                    todo = [unseen.pop()]
                    while todo:
                        reached = graph[todo.pop()] & unseen
                        unseen.difference_update(reached); todo.extend(reached)
                squares = sum({(x, y), (x+1, y), (x, y+1), (x+1, y+1)} <= live
                              for x in range(n-1) for y in range(2))
                assert sum(map(len, graph.values()))//2-len(live)+components == squares
    return dict(counts)


def audit_corpus(path):
    records = json.loads(gzip.decompress(path.read_bytes()))["nodes"]
    geometries = {int(row["mask"], 16) for row in records}
    hits, branches = {}, {}
    for mask in geometries:
        graph = adjacency(mask)
        spine = caterpillar_spine(graph)
        if spine is not None:
            hits[mask] = interval_value(spine)
            if max(map(len, graph.values())) > 2:
                branches[mask] = spine
    for row in records:
        mask = int(row["mask"], 16)
        if mask in hits:
            assert (hits[mask] == row["k"]) is row["equal"]
    return dict(geometries=len(geometries), trees=len(hits), branching_trees=len(branches),
                branching_queries=sum(int(row["mask"], 16) in branches for row in records),
                max_vertices=max(map(int.bit_count, branches)), max_spine=max(map(len, branches.values())),
                maximum_replaced_tree_cost=max(int(row["cost"]) for row in records
                                               if int(row["mask"], 16) in branches))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    root = Path(__file__).resolve().parents[1]
    corpus = root / "data/proofs/strip3_equality/strip3_equality_target_tree_min_26.json.gz"
    sources = [Path(__file__).resolve(), root / "tests/test_strip3_caterpillar.py",
               root / "src/cfast/strip3_caterpillar.c", root / "src/cfast/strip3_target.c", corpus]
    report = dict(format="strip3-caterpillar-audit-v1", local=local_check(), corpus=audit_corpus(corpus),
                  provenance=dict(commit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root,
                                                                 text=True).strip(),
                                  sha256={str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
                                          for p in sources}),
                  solver_processes=0, wall_s=round(time.monotonic()-started, 6))
    if args.output:
        args.output.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report))
