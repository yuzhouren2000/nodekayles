#!/usr/bin/env python3
"""Check a complete boundary DAG certificate by independent vertex-set deletion.

The verifier does not import the interface compiler or its SG solver. It checks
all options in each boundary-attached induced subset of each supplied graph,
including the SG of every detached component; a matching root is then covered
by the arbitrary-context replacement theorem in notes/18.
"""

from functools import cache
import json
from pathlib import Path
import sys


def freeze(value):
    return tuple(map(freeze, value)) if isinstance(value, list) else value


def verify(proof):
    assert proof["format"] == "node-kayles-boundary-dag-v1"
    nodes = tuple(map(freeze, proof["nodes"]))
    assert nodes[0] == ((), (), (), ())
    assert len(set(nodes)) == len(nodes), "duplicate interned nodes"
    lookup = {node: i for i, node in enumerate(nodes)}
    for i, (live, direct, moves, internal) in enumerate(nodes):
        assert live == tuple(sorted(set(live))) and all(x > 0 for x in live)
        assert tuple(label for label, _ in direct) == live
        assert tuple(label for label, _ in moves) == live
        assert internal == tuple(sorted(set(internal)))
        for label, state in direct + moves:
            child, nimber = state
            assert 0 <= child < i and type(nimber) is int and nimber >= 0
            assert label not in nodes[child][0]
            assert set(nodes[child][0]) <= set(live)
        for label, (child, _) in direct:
            assert nodes[child][0] == tuple(x for x in live if x != label)
        for removed, (child, nimber) in internal:
            assert removed == tuple(sorted(set(removed)))
            assert 0 <= child < i and type(nimber) is int and nimber >= 0
            assert nodes[child][0] == tuple(x for x in live if x not in removed)
    checked = 0
    actual = []
    for item in proof["graphs"]:
        adjacency, labels = item["adjacency"], item["labels"]
        n = len(adjacency)
        assert n > 0 and len(labels) == n
        positive = [x for x in labels if x]
        assert len(positive) == len(set(positive))
        assert all(type(x) is int and x >= 0 for x in labels)
        edges = {v: frozenset(u for u in range(n) if adjacency[v] >> u & 1)
                 for v in range(n)}
        for v in range(n):
            assert 0 <= adjacency[v] < 1 << n and v not in edges[v]
            assert all(v in edges[u] for u in edges[v])

        def component(position, start):
            seen, pending = {start}, [start]
            while pending:
                v = pending.pop()
                new = (edges[v] & position) - seen
                seen.update(new)
                pending.extend(new)
            return frozenset(seen)

        @cache
        def raw_sg(position):
            if not position:
                return 0
            outcomes = {raw_sg(position - edges[v] - {v}) for v in position}
            result = 0
            while result in outcomes:
                result += 1
            return result

        @cache
        def state(position):
            nonlocal checked
            checked += 1
            boundary = {v for v in position if labels[v]}
            if not boundary:
                return 0, raw_sg(position)
            attached = frozenset().union(*(component(position, v) for v in boundary))
            if attached != position:
                child, nimber = state(attached)
                return child, nimber ^ raw_sg(position - attached)
            direct, moves, internal = [], [], set()
            for v in sorted(position):
                after = state(position - edges[v] - {v})
                if labels[v]:
                    direct.append((labels[v], state(position - {v})))
                    moves.append((labels[v], after))
                else:
                    killed = tuple(sorted(labels[u] for u in edges[v] & boundary))
                    internal.add((killed, after))
            key = (tuple(sorted(labels[v] for v in boundary)), tuple(sorted(direct)),
                   tuple(sorted(moves)), tuple(sorted(internal)))
            assert key in lookup, "certificate omits or alters an actual option"
            return lookup[key], 0

        actual.append(state(frozenset(range(n))))
    assert actual == list(map(tuple, proof["states"]))
    assert len(actual) >= 2 and len(set(actual)) == 1, "roots do not match"
    used = {state[0] for state in actual}
    pending = list(used)
    while pending:
        node = nodes[pending.pop()]
        children = {state[0] for _, state in node[1] + node[2] + node[3]}
        new = children - used
        used.update(new)
        pending.extend(new)
    assert used == set(range(len(nodes))), "unreferenced certificate nodes"
    return {"nodes": len(nodes), "induced_states_checked": checked,
            "common_state": actual[0]}


if __name__ == "__main__":
    source = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(
        "data/proofs/kayles_boundary_replacement.json")
    print(json.dumps(verify(json.loads(source.read_text())), ensure_ascii=False))
