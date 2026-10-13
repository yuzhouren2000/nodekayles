#!/usr/bin/env python3
"""Exact recursive interfaces for Node Kayles under vertex identification.

An interface retains every internal option while compressing detached components
by xor. Equal forms can replace one another in any context on the labelled
boundary. Forms are sufficient, not asserted minimal or finite in number.
"""

import argparse
from functools import cache
import json
from pathlib import Path


def mex(values):
    values = set(values)
    result = 0
    while result in values:
        result += 1
    return result


def validate(graph, labels):
    if not graph or len(labels) != len(graph):
        raise ValueError("nonempty graph and one label per vertex required")
    active = [label for label in labels if label]
    if any(type(label) is not int or label < 0 for label in labels):
        raise ValueError("boundary labels must be nonnegative integers")
    if len(set(active)) != len(active):
        raise ValueError("positive boundary labels must be unique")
    full = (1 << len(graph)) - 1
    for v, neighbours in enumerate(graph):
        if type(neighbours) is not int or neighbours < 0 or neighbours & ~full:
            raise ValueError("invalid adjacency bitset")
        if neighbours >> v & 1:
            raise ValueError("loops are not allowed")
        for u in range(len(graph)):
            if bool(neighbours >> u & 1) != bool(graph[u] >> v & 1):
                raise ValueError("adjacency must be symmetric")


def induced(graph, mask):
    vertices = [v for v in range(len(graph)) if mask >> v & 1]
    return tuple(sum(1 << j for j, u in enumerate(vertices)
                     if graph[v] >> u & 1) for v in vertices)


def reach(graph, seeds):
    seen = pending = seeds
    while pending:
        bit = pending & -pending
        pending ^= bit
        new = graph[bit.bit_length() - 1] & ~seen
        seen |= new
        pending |= new
    return seen


@cache
def grundy(graph):
    if not graph:
        return 0
    full = (1 << len(graph)) - 1
    first = reach(graph, 1)
    if first != full:
        return grundy(induced(graph, first)) ^ grundy(induced(graph, full ^ first))
    return mex(grundy(induced(graph, full & ~(neighbours | (1 << v))))
               for v, neighbours in enumerate(graph))


# p-a-b-c-d-a is a square with pendant q at d; roots are p,a.
SQUARE = (2, 13, 18, 50, 12, 8)
# p-a-z-d-q is a path with an additional pendant t at z.
FORK = (2, 5, 26, 4, 36, 16)
PORTS = (1, 2, 0, 0, 0, 0)


class BoundaryAlgebra:
    """Hash-cons a complete boundary game DAG, then evaluate its products.

    A state is (form_id, detached_nimber). Form zero has no live boundary;
    its state is just a nimber. Every other node is
      (live_labels, direct_deletions, boundary_moves, internal_moves).
    Entries in both boundary tables are (label, child_state). An internal
    entry is (removed_boundary_labels, child_state).
    """

    def __init__(self):
        self.nodes = [((), (), (), ())]
        self.intern = {self.nodes[0]: 0}
        self._form = cache(self._form)
        self.delete = cache(self.delete)
        self.product = cache(self.product)

    def state(self, graph, labels):
        graph, labels = tuple(graph), tuple(labels)
        validate(graph, labels)
        return self._form(graph, labels)

    def _form(self, graph, labels):
        full = (1 << len(graph)) - 1
        boundary = sum(1 << i for i, label in enumerate(labels) if label)
        if not boundary:
            return 0, grundy(graph)
        attached = reach(graph, boundary)
        if attached != full:
            attached_labels = tuple(labels[i] for i in range(len(graph))
                                    if attached >> i & 1)
            child, offset = self._form(induced(graph, attached), attached_labels)
            return child, offset ^ grundy(induced(graph, full ^ attached))

        def child(mask):
            return self._form(induced(graph, mask),
                              tuple(labels[i] for i in range(len(graph)) if mask >> i & 1))

        direct, moves, internal = [], [], set()
        for v, neighbours in enumerate(graph):
            after = child(full & ~(neighbours | (1 << v)))
            if labels[v]:
                direct.append((labels[v], child(full & ~(1 << v))))
                moves.append((labels[v], after))
            else:
                removed = tuple(sorted(labels[u] for u in range(len(graph))
                                       if labels[u] and neighbours >> u & 1))
                internal.add((removed, after))
        node = (tuple(sorted(label for label in labels if label)),
                tuple(sorted(direct)), tuple(sorted(moves)), tuple(sorted(internal)))
        if node not in self.intern:
            self.intern[node] = len(self.nodes)
            self.nodes.append(node)
        return self.intern[node], 0

    def delete(self, form, removed):
        """Directly delete these labels, accumulating detached nimbers."""
        current, offset = form, 0
        for label in sorted(set(removed)):
            transition = dict(self.nodes[current][1]).get(label)
            if transition is not None:
                current, extra = transition
                offset ^= extra
        return current, offset

    def product(self, left, right):
        """SG after identifying all equally labelled boundary vertices."""
        if left > right:
            return self.product(right, left)
        a, b = self.nodes[left], self.nodes[right]
        if a[0] != b[0]:
            raise ValueError("product needs exactly the same live boundary labels")
        if not a[0]:
            return 0
        options = set()
        # An internal move only directly removes its boundary neighbours on
        # the other side; it cannot hit the other side's internal vertices.
        for source, context in ((a, right), (b, left)):
            for removed, (child, offset) in source[3]:
                other, extra = self.delete(context, removed)
                options.add(offset ^ extra ^ self.product(child, other))
        ma, mb = dict(a[2]), dict(b[2])
        for label in a[0]:
            ca, xa = ma[label]
            cb, xb = mb[label]
            live = set(self.nodes[ca][0]) & set(self.nodes[cb][0])
            da, ya = self.delete(ca, tuple(set(self.nodes[ca][0]) - live))
            db, yb = self.delete(cb, tuple(set(self.nodes[cb][0]) - live))
            options.add(xa ^ xb ^ ya ^ yb ^ self.product(da, db))
        return mex(options)

    def glued_value(self, left_state, right_state):
        a, x = left_state
        b, y = right_state
        return x ^ y ^ self.product(a, b)

    def proof(self, graphs):
        states = [self.state(graph, labels) for graph, labels in graphs]
        return {"format": "node-kayles-boundary-dag-v1", "states": states,
                "graphs": [{"adjacency": graph, "labels": labels}
                           for graph, labels in graphs], "nodes": self.nodes}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    algebra = BoundaryAlgebra()
    proof = algebra.proof(((SQUARE, PORTS), (FORK, PORTS)))
    assert proof["states"][0] == proof["states"][1]
    proof["reachable_three_row_examples"] = [
        {"width": 4, "mask": 2484, "root_cell": 2, "kind": "square"},
        {"width": 3, "mask": 189, "root_cell": 0, "kind": "fork"}]
    output = json.dumps(proof, indent=2) + "\n"
    if args.output:
        args.output.write_text(output)
    else:
        print(output, end="")


if __name__ == "__main__":
    main()
