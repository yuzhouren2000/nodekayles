"""Check a Node Kayles equality certificate WITHOUT a solver or SG oracle.

Use coordinates for deletion, components and rectangle isomorphisms. Leaves
use only the path recurrence, cycle rule, safe involution or vertex bound.
This establishes certificate validity/cost, not minimum-DAG optimality.
"""
from collections import Counter
from functools import cache
import gzip
import json
from pathlib import Path
import sys


@cache
def cells(mask):
    return frozenset((v//3, v%3) for v in range(mask.bit_length()) if mask >> v & 1)


def mask_of(points):
    return sum(1 << (3*c+r) for c, r in points)


def components(points):
    unseen = set(points)
    while unseen:
        start = unseen.pop(); reached = {start}; todo = [start]
        while todo:
            c, r = todo.pop()
            for v in ((c-1, r), (c+1, r), (c, r-1), (c, r+1)):
                if v in unseen:
                    unseen.remove(v); reached.add(v); todo.append(v)
        yield reached


@cache
def canonical(mask):
    points = cells(mask)
    assert points
    c0 = min(c for c, _ in points); r0 = min(r for _, r in points)
    pts = {(c-c0, r-r0) for c, r in points}
    w = max(c for c, _ in pts)+1; h = max(r for _, r in pts)+1
    bases = [pts]
    if w == h == 3: bases.append({(r, c) for c, r in pts})
    return min(mask_of({(w-1-c if fc else c, h-1-r if fr else r) for c, r in base})
               for base in bases for fc in (False, True) for fr in (False, True))


@cache
def normal_parts(mask):
    counts = Counter(canonical(mask_of(part)) for part in components(cells(mask)))
    return tuple(sorted((m for m, count in counts.items() if count%2), key=lambda m: (-m.bit_count(), m)))


def after_move(mask, v):
    assert mask >> v & 1
    c, r = divmod(v, 3)
    return mask_of({(x, y) for x, y in cells(mask) if abs(x-c)+abs(y-r) > 1})


@cache
def options(mask):
    return {normal_parts(after_move(mask, v)) for v in range(mask.bit_length()) if mask >> v & 1}


@cache
def path_value(n):
    values = {path_value(max(i-1, 0)) ^ path_value(max(n-i-2, 0)) for i in range(n)}
    g = 0
    while g in values: g += 1
    return g


def axiom_value(mask, rule):
    points = cells(mask); assert len(tuple(components(points))) == 1
    degrees = [sum(abs(c-d)+abs(r-s) == 1 for d, s in points) for c, r in points]
    if rule == 2:
        assert max(degrees) <= 2 and (min(degrees) < 2 or len(points) == 1)
        return path_value(len(points))
    if rule == 3:
        assert min(degrees) == max(degrees) == 2
        return int(path_value(len(points)-3) == 0)
    assert rule == 4
    w = max(c for c, _ in points)+1; h = max(r for _, r in points)+1
    image = {(w-1-c, h-1-r) for c, r in points}
    assert image == points
    assert all((c, r) != (w-1-c, h-1-r) and abs(2*c-w+1)+abs(2*r-h+1) > 1 for c, r in points)
    return 0


def verify(report):
    assert report["format"] == "strip3-equality-proof-v1"
    n, g = report["n"], report["SG"]
    assert type(n) is int and 1 <= n <= 42 and type(g) is int and 0 <= g <= 3*n
    records = report["nodes"]
    ids = [r["id"] for r in records]
    assert sorted(ids) == list(range(1, len(records)+1))
    nodes = {r["id"]: r for r in records}
    active, costs = set(), {}
    sum_count = edges = 0

    def query(node_id, mask=None, k=None, truth=None):
        assert type(node_id) is int and node_id in nodes
        node = nodes[node_id]
        if mask is not None: assert int(node["mask"], 16) == mask
        if k is not None: assert node["k"] == k
        if truth is not None: assert node["equal"] is truth
        return check(node_id)

    def check_sum(record, expected_parts, k, truth):
        nonlocal sum_count, edges
        parts = tuple(int(m, 16) for m in record["parts"])
        assert parts == expected_parts
        assert record["k"] == k and record["equal"] is truth
        children = record["children"]; assert len(children) == len(parts)
        pivot = record.get("pivot", 0)
        assert type(pivot) is int and (0 <= pivot < len(parts) if parts else pivot == 0)
        sum_count += 1; edges += len(children)
        cost, offset = 1, 0
        for i in range(len(parts)):
            if i == pivot: continue
            child = nodes[children[i]]
            cost += query(children[i], parts[i], truth=True)
            offset ^= child["k"]
        if parts:
            cost += query(children[pivot], parts[pivot], k^offset, truth)
        else:
            assert truth is (k == 0)
        return cost

    def check(node_id):
        nonlocal edges
        if node_id in costs: return costs[node_id]
        assert node_id not in active, "cyclic proof"
        active.add(node_id)
        node = nodes[node_id]
        mask, k, truth, rule = int(node["mask"], 16), node["k"], node["equal"], node["rule"]
        assert mask > 0 and mask.bit_length() <= 126 and canonical(mask) == mask
        assert len(tuple(components(cells(mask)))) == 1
        assert type(k) is int and 0 <= k <= 127 and type(truth) is bool
        cost = 1
        if rule == 1:
            assert k > mask.bit_count() and truth is False
        elif rule in (2, 3, 4):
            assert truth is (axiom_value(mask, rule) == k)
        elif rule == 5:
            assert truth is True
            expected = options(mask)
            observed = [tuple(int(m, 16) for m in s["parts"]) for s in node["options"]]
            assert len(observed) == len(set(observed)) and set(observed) == expected, "missing universal move"
            for record, parts in zip(node["options"], observed):
                cost += check_sum(record, parts, k, False)
            lower = node["lower"]; assert len(lower) == k
            edges += len(lower)
            for j, child in enumerate(lower): cost += query(child, mask, j, False)
        elif rule == 6:
            assert truth is False
            v = node["move"]; assert type(v) is int and 0 <= v < 126 and mask >> v & 1
            cost += check_sum(node["witness"], normal_parts(after_move(mask, v)), k, True)
        elif rule == 7:
            assert truth is False and type(node["lower_k"]) is int and 0 <= node["lower_k"] < k
            cost += query(node["child"], mask, node["lower_k"], True); edges += 1
        else:
            raise AssertionError("unknown rule")
        assert str(cost) == node["cost"], "incorrect unfolded cost"
        active.remove(node_id); costs[node_id] = cost
        return cost

    cost = query(report["root"], (1 << (3*n))-1, g, True)
    assert len(costs) == len(records), "unreachable certificate records"
    assert str(cost) == report["tree_cost"]
    return dict(n=n, SG=g, queries=len(records), component_masks=len({r["mask"] for r in records}),
                sum_rules=sum_count, edges=edges, tree_cost=str(cost))


def main():
    for path in sys.argv[1:]:
        raw = Path(path).read_bytes()
        report = json.loads(gzip.decompress(raw) if path.endswith(".gz") else raw)
        print(json.dumps(dict(file=path, **verify(report))))


if __name__ == "__main__":
    main()
