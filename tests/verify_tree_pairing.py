"""Standalone verifier for structural tree-pairing certificates/refutations.

Imports neither the constructor, matching search nor any SG solver. SAT
checks the actual response invariant plus the two-halves witness. UNSAT
reconstructs every possible balanced cut and checks rooted isomorphism using
parenthesis words, independently of the constructor's interned tree types.

Run: python3 tests/verify_tree_pairing.py data/proofs/tree_pairing_examples.json
"""

import json
import sys
from pathlib import Path


def graph_from_record(record):
    if record.get("format") != "nodekayles-tree-pairing-v1":
        raise ValueError("unknown certificate format")
    n = record["order"]
    if type(n) is not int or n < 1:
        raise ValueError("invalid tree order")
    adjacency = {v: set() for v in range(n)}
    edges = set()
    for edge in record["edges"]:
        if len(edge) != 2:
            raise ValueError("invalid edge")
        a, b = edge
        if any(type(v) is not int or v not in adjacency for v in (a, b)):
            raise ValueError("invalid vertex")
        if a == b or frozenset(edge) in edges:
            raise ValueError("loop or duplicate edge")
        edges.add(frozenset(edge))
        adjacency[a].add(b)
        adjacency[b].add(a)
    if len(edges) != n - 1 or len(component(adjacency, 0)) != n:
        raise ValueError("record is not a tree")
    return adjacency


def component(adjacency, start, cut=None):
    seen = {start}
    todo = [start]
    for v in todo:
        for u in adjacency[v]:
            if cut is not None and {u, v} == set(cut):
                continue
            if u not in seen:
                seen.add(u)
                todo.append(u)
    return seen


def rooted_word(adjacency, vertices, root):
    """Canonical balanced-parenthesis word for a rooted tree."""
    parent = {root: None}
    traversal = [root]
    children = {}
    for v in traversal:
        children[v] = adjacency[v].intersection(vertices) - {parent[v]}
        for u in children[v]:
            parent[u] = v
            traversal.append(u)
    words = {}
    for v in reversed(traversal):
        words[v] = "(" + "".join(sorted(words.pop(u) for u in children[v])) + ")"
    return words[root]


def verify_pairs(adjacency, pairs):
    """Simulate every attack and prescribed response with set deletions."""
    phi = {}
    for pair in pairs:
        if len(pair) != 2:
            raise ValueError("invalid pair")
        a, b = pair
        if any(type(v) is not int or v not in adjacency for v in (a, b)):
            raise ValueError("invalid paired vertex")
        if a == b or a in phi or b in phi or b in adjacency[a]:
            raise ValueError("pairs overlap or prescribe an illegal response")
        phi[a], phi[b] = b, a
    if set(phi) != set(adjacency):
        raise ValueError("pairing does not cover every vertex")
    for attack in adjacency:
        after_attack = set(adjacency) - {attack} - adjacency[attack]
        reply = phi[attack]
        if reply not in after_attack:
            raise ValueError("response is not legal")
        after_reply = after_attack - {reply} - adjacency[reply]
        if {phi[v] for v in after_reply} != after_reply:
            raise ValueError("two moves do not restore the pairing invariant")
    return phi


def verify_result(record):
    adjacency = graph_from_record(record)
    n = len(adjacency)
    status = record["status"]
    if status == "SAT":
        phi = verify_pairs(adjacency, record["pairs"])
        a, b = record["bridge"]
        if any(type(v) is not int or v not in adjacency for v in (a, b)):
            raise ValueError("invalid bridge")
        if b not in adjacency[a]:
            raise ValueError("bridge is not an edge")
        left = component(adjacency, a, (a, b))
        right = set(adjacency) - left
        if len(left) != len(right) or {phi[v] for v in left} != right:
            raise ValueError("pairing does not identify the two bridge halves")
        if phi[a] not in adjacency[b] or phi[b] not in adjacency[a]:
            raise ValueError("bridge endpoints do not correspond to adjacent vertices")
        for v in left:
            if {phi[u] for u in adjacency[v] & left} != adjacency[phi[v]] & right:
                raise ValueError("halves are not isomorphic under the proposed pairing")
        return True
    if status != "UNSAT":
        raise ValueError("UNKNOWN is not a certificate")
    reason = record["reason"]
    if n % 2:
        if reason != "odd-order":
            raise ValueError("wrong odd-order reason")
        return True
    balanced = []
    for a in adjacency:
        for b in adjacency[a]:
            if a < b:
                left = component(adjacency, a, (a, b))
                if 2 * len(left) == n:
                    balanced.append((a, b, left))
    if not balanced:
        if reason != "no-balanced-edge":
            raise ValueError("wrong balanced-edge reason")
        return True
    if reason != "incompatible-rooted-halves" or len(balanced) != 1:
        raise ValueError("wrong rooted-halves reason")
    a, b, left = balanced[0]
    if set(record["bridge"]) != {a, b}:
        raise ValueError("wrong balanced bridge")
    right = set(adjacency) - left
    # An isomorphism must map a to a neighbour of b in the other half.
    # Check every such neighbour; the word characterises rooted isomorphism.
    target = rooted_word(adjacency, left, a)
    for root in adjacency[b] & right:
        if target == rooted_word(adjacency, right, root):
            raise ValueError("a structural pairing exists; false refutation")
    return True


def main():
    if not sys.argv[1:]:
        print(__doc__)
        return 2
    for filename in sys.argv[1:]:
        data = json.loads(Path(filename).read_text())
        records = data["examples"] if "examples" in data else [data]
        for record in records:
            verify_result(record)
        print(f"verified {len(records)} structural tree records in {filename}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
