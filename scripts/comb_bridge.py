#!/usr/bin/env python3
"""Exact arbitrary-core / alternating-comb / arbitrary-core bridge algebra.

Both roots are vertex zero in their own cores. Length zero means a direct
root-to-root edge. The eventual period-four proof is in notes/22; it uses
certified one-ended core tables, never a bare SG boundary summary.
"""

import argparse
from functools import cache
from itertools import chain
import json
from pathlib import Path

from kayles_boundary_algebra import grundy, induced
from strip3_defect_algebra import RootedComb, full_head, stabilization_bound


def bridge_graph(left, right, length, phase):
    """Construct the actual graph, for finite checks and length zero."""
    if type(length) is not int or length < 0 or phase not in (0, 1):
        raise ValueError("invalid bridge parameters")
    left, right = tuple(left), tuple(right)
    left_order = len(left)
    graph = list(left) + [ns << left_order for ns in right]

    def edge(u, v):
        graph[u] |= 1 << v
        graph[v] |= 1 << u

    if not length:
        edge(0, left_order)
        return tuple(graph)
    first = len(graph)
    graph.extend([0] * length)
    edge(0, first)
    edge(first + length - 1, left_order)
    for i in range(length - 1):
        edge(first + i, first + i + 1)
    for i in range(length):
        if phase ^ (i % 2):
            graph.append(0)
            edge(first + i, len(graph) - 1)
    return tuple(graph)


class CombBridge:
    def __init__(self, left, right):
        # Both phases are needed at the right end after reversing the bridge.
        self.left = tuple(RootedComb(left, p).stabilize() for p in (0, 1))
        self.right = tuple(RootedComb(right, p).stabilize() for p in (0, 1))
        self.starts = {model: model.period_start() for model in self.left + self.right}
        self.left_start = max(self.starts[model] for model in self.left)
        self.right_start = max(self.starts[model] for model in self.right)
        self.period_start = self.left_start + self.right_start + 6
        self.order_bound = (stabilization_bound(len(left))
                            + stabilization_bound(len(right)) + 6)
        self._value = cache(self._value)

    def _tail(self, model, mask, length):
        start = self.starts[model]
        n = length if length < start else start + (length - start) % 4
        return model.values[mask][n]

    def _positions(self, length):
        if length < self.period_start:
            return range(length)
        left_end = min(length, self.left_start + 5)
        right_start = max(left_end, length - self.right_start - 2)
        return chain(range(left_end), range(right_start, length))

    def _options(self, left_mask, right_mask, length, phase):
        left = self.left[phase]
        right_phase = phase ^ ((length - 1) % 2)
        right = self.right[right_phase]
        tail_left = lambda n: self._tail(left, left_mask, n)
        tail_right = lambda n: self._tail(right, right_mask, n)
        options = []
        for side, model, mask, other_tail, first_tooth in (
                ("left", left, left_mask, tail_right, phase),
                ("right", right, right_mask, tail_left, right_phase)):
            for v, ns in enumerate(model.adjacency):
                if not mask >> v & 1:
                    continue
                rest = mask & ~model.closed[v]
                if v == 0:
                    value = model.core_value(rest) ^ other_tail(length - 1) ^ first_tooth
                elif ns & 1:
                    value = model.core_value(rest) ^ other_tail(length)
                else:
                    masks = ((rest, right_mask) if side == "left"
                             else (left_mask, rest))
                    value = self._value(*masks, length, phase)
                options.append(((side, v), value))
        for i in self._positions(length):
            tooth = phase ^ (i % 2)
            if tooth:
                options.append((("leaf", i), tail_left(i) ^ tail_right(length - i - 1)))
            if length == 1:
                value = left.core_value(left_mask & ~1) ^ right.core_value(right_mask & ~1)
            elif i == 0:
                value = left.core_value(left_mask & ~1) ^ tail_right(length - 2) ^ (1 - phase)
            elif i == length - 1:
                value = right.core_value(right_mask & ~1) ^ tail_left(length - 2) ^ (1 - right_phase)
            else:
                # Both neighbours' orphan teeth either cancel or do not exist.
                value = tail_left(i - 1) ^ tail_right(length - i - 2)
            options.append((("spine", i), value))
        return options

    def _value(self, left_mask, right_mask, length, phase):
        if not length:
            left = induced(self.left[0].adjacency, left_mask)
            right = induced(self.right[0].adjacency, right_mask)
            return grundy(bridge_graph(left, right, 0, phase))
        values = {value for _, value in self._options(left_mask, right_mask, length, phase)}
        result = 0
        while result in values:
            result += 1
        return result

    def value(self, length, phase=0):
        if type(length) is not int or length < 0 or phase not in (0, 1):
            raise ValueError("invalid bridge parameters")
        return self._value(self.left[0].full, self.right[0].full, length, phase)

    def options(self, length, phase=0):
        self.value(length, phase)  # Validate and compute all smaller core pairs.
        if not length:
            raise ValueError("move representatives require a nonempty bridge")
        return self._options(self.left[0].full, self.right[0].full, length, phase)

    def certificate(self):
        """Record generic core tables, including every root-preserving state."""
        tables = []
        for side, models in (("left", self.left), ("right", self.right)):
            for phase, model in enumerate(models):
                start = self.starts[model]
                tables.append({"side": side, "phase": phase,
                               "adjacency": model.adjacency, "start": start,
                               "states": [{"mask": mask,
                                           "values": model.values[mask][:start + 10]}
                                          for mask in model.states]})
        return {"format": "comb-bridge-period-v1", "left_start": self.left_start,
                "right_start": self.right_start, "period_start": self.period_start,
                "order_bound": self.order_bound, "one_ended_tables": tables,
                "period_values": {str(p): [self.value(self.period_start + i, p)
                                          for i in range(4)] for p in (0, 1)}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--left-head", type=int, default=2)
    parser.add_argument("--right-head", type=int, default=2)
    parser.add_argument("--length", type=int, default=10 ** 12)
    parser.add_argument("--phase", type=int, choices=(0, 1), default=0)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if min(args.left_head, args.right_head) < 1 or args.length < 1:
        parser.error("positive head widths and bridge length required")
    left = full_head(args.left_head)[1]
    right = full_head(args.right_head)[1]
    model = CombBridge(left, right)
    options = model.options(args.length, args.phase)
    winning = next((move for move, value in options if value == 0), None)
    print(json.dumps({"left_head": args.left_head, "right_head": args.right_head,
                      "length": args.length, "phase": args.phase,
                      "SG": model.value(args.length, args.phase),
                      "period_start": model.period_start, "order_bound": model.order_bound,
                      "period_values": [model.value(model.period_start + i, args.phase)
                                        for i in range(4)],
                      "representative_options": len(options), "winning_move": winning}))
    if args.output:
        args.output.write_text(json.dumps(model.certificate(), indent=2) + "\n")


if __name__ == "__main__":
    main()
