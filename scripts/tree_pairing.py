"""Recognise every fixed-pairing certificate on a finite tree.

A tree admits a fixed response pairing iff it consists of two isomorphic
trees joined by one edge whose endpoints correspond to adjacent vertices
under the isomorphism. See notes/11_tree_pairing.md for the complete proof.

Examples:
  python3 scripts/tree_pairing.py --spider 1 3 5 --output /tmp/tree.json
  python3 scripts/tree_pairing.py --graph graph.json --output /tmp/tree.json

The input graph JSON contains order and edges; vertices are 0..order-1.
UNSAT excludes fixed pairing only, not other winning strategies.
"""

import argparse
import json
from pathlib import Path

FORMAT = "nodekayles-tree-pairing-v1"


def from_edges(order, edges):
    if type(order) is not int or order < 1:
        raise ValueError("order must be a positive integer")
    adjacency = [set() for _ in range(order)]
    for edge in edges:
        if len(edge) != 2:
            raise ValueError("an edge must have two endpoints")
        a, b = edge
        if any(type(v) is not int or not 0 <= v < order for v in (a, b)):
            raise ValueError("invalid vertex")
        if a == b or b in adjacency[a]:
            raise ValueError("loop or duplicate edge")
        adjacency[a].add(b)
        adjacency[b].add(a)
    return adjacency


def validate_tree(adjacency):
    n = len(adjacency)
    if n < 1:
        raise ValueError("a tree must be nonempty")
    for v, neighbours in enumerate(adjacency):
        for u in neighbours:
            if type(u) is not int or not 0 <= u < n or u == v:
                raise ValueError("invalid adjacency")
            if v not in adjacency[u]:
                raise ValueError("adjacency is not symmetric")
    if sum(map(len, adjacency)) != 2 * (n - 1):
        raise ValueError("graph does not have tree edge count")
    seen = {0}
    stack = [0]
    while stack:
        for u in adjacency[stack.pop()]:
            if u not in seen:
                seen.add(u)
                stack.append(u)
    if len(seen) != n:
        raise ValueError("graph is disconnected")


def spider(arms):
    """Three paths incident with vertex 0; arm lengths count edges."""
    if len(arms) != 3 or any(type(t) is not int or t < 1 for t in arms):
        raise ValueError("three positive arm lengths are required")
    adjacency = [set()]
    for length in arms:
        previous = 0
        for _ in range(length):
            v = len(adjacency)
            adjacency.append({previous})
            adjacency[previous].add(v)
            previous = v
    return adjacency


def double_along_edge(adjacency, edge):
    """Two copies of a tree, connected from a in copy 0 to b in copy 1."""
    validate_tree(adjacency)
    a, b = edge
    if a not in range(len(adjacency)) or b not in adjacency[a]:
        raise ValueError("the marked endpoints must be adjacent")
    n = len(adjacency)
    result = [set(ns) for ns in adjacency]
    result.extend({v + n for v in ns} for ns in adjacency)
    result[a].add(b + n)
    result[b + n].add(a)
    return result


def _rooted_types(adjacency, vertices, root, intern):
    """Iterative AHU rooted-tree types, interned across both halves."""
    parent = {root: None}
    order = [root]
    children = {}
    for v in order:
        children[v] = sorted(u for u in adjacency[v]
                             if u in vertices and u != parent[v])
        for u in children[v]:
            parent[u] = v
            order.append(u)
    types = {}
    for v in reversed(order):
        key = tuple(sorted(types[u] for u in children[v]))
        if key not in intern:
            intern[key] = len(intern)
        types[v] = intern[key]
    return children, types


def classify_tree(adjacency):
    """Return a structural SAT certificate or a theorem-based UNSAT result.

    No matching enumeration or SG search is used. Rooted-tree isomorphism is
    tried at each neighbour of one bridge endpoint. Sorting adjacency lists
    makes the output deterministic; traversal is iterative for long paths.
    """
    validate_tree(adjacency)
    n = len(adjacency)
    result = {"format": FORMAT, "order": n,
              "edges": [[a, b] for a, ns in enumerate(adjacency)
                        for b in sorted(ns) if a < b]}
    if n % 2:
        return {**result, "status": "UNSAT", "reason": "odd-order"}

    parent = {0: None}
    order = [0]
    for v in order:
        for u in sorted(adjacency[v]):
            if u != parent[v]:
                parent[u] = v
                order.append(u)
    sizes = {v: 1 for v in order}
    for v in reversed(order[1:]):
        sizes[parent[v]] += sizes[v]
    balanced = [v for v in order[1:] if 2 * sizes[v] == n]
    if not balanced:
        return {**result, "status": "UNSAT", "reason": "no-balanced-edge"}
    assert len(balanced) == 1
    a = balanced[0]
    b = parent[a]
    left = {a}
    stack = [a]
    while stack:
        for u in adjacency[stack.pop()]:
            if u != b and u not in left:
                left.add(u)
                stack.append(u)
    right = set(range(n)) - left
    intern = {}
    left_children, left_types = _rooted_types(adjacency, left, a, intern)
    for root in sorted(adjacency[b] & right):
        right_children, right_types = _rooted_types(adjacency, right, root, intern)
        if left_types[a] != right_types[root]:
            continue
        mapping = []
        pending = [(a, root)]
        while pending:
            v, w = pending.pop()
            mapping.append([v, w])
            lc = sorted(left_children[v], key=lambda x: (left_types[x], x))
            rc = sorted(right_children[w], key=lambda x: (right_types[x], x))
            pending.extend(zip(lc, rc))
        return {**result, "status": "SAT", "bridge": [a, b],
                "pairs": sorted(mapping)}
    return {**result, "status": "UNSAT", "reason": "incompatible-rooted-halves",
            "bridge": [a, b]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--spider", type=int, nargs=3, metavar=("A", "B", "C"))
    group.add_argument("--graph", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.spider:
        adjacency = spider(args.spider)
    else:
        graph = json.loads(args.graph.read_text())
        adjacency = from_edges(graph["order"], graph["edges"])
    result = classify_tree(adjacency)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"tree of order {len(adjacency)}: {result['status']}; wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
