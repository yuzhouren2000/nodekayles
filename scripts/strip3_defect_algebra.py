#!/usr/bin/env python3
"""Exact interval algebra and certified period lifting for outer-row defects.

An alternating comb has a backbone path and one pendant leaf at alternating
backbone vertices. Its phase is 1 when the first vertex has a leaf. A rooted
core attaches to the first backbone vertex by ONE edge. No equality of bare
SG values is used as a replacement rule. See notes/14_outer_defect_algebra.md.
"""

import argparse
from functools import cache
import json
from pathlib import Path


def mex(values):
    result = 0
    while result in values:
        result += 1
    return result


def comb_value(phase, length):
    if phase not in (0, 1) or length < 0:
        raise ValueError("invalid alternating comb")
    return ((0, 1, 2, 1), (0, 1, 2, 3))[phase][length % 4]


class TwoRowIntervals:
    """Solve any induced two-row strip with interval and endpoint deletions.

    Column masks have bit 0 for the upper row and bit 1 for the lower row.
    The immutable interior word is retained, rather than summarized by SG.
    At most 16 states per nonempty interval; each state has O(n) options.
    """

    def __init__(self, columns):
        self.columns = tuple(columns)
        if any(type(c) is not int or not 0 <= c <= 3 for c in self.columns):
            raise ValueError("column masks must be integers in 0..3")
        self.value = cache(self._value)

    def _value(self, left, right, delete_left=0, delete_right=0):
        if left >= right:
            return 0
        options = set()
        for column in range(left, right):
            playable = self.columns[column]
            if column == left:
                playable &= ~delete_left
            if column == right - 1:
                playable &= ~delete_right
            for row in (1, 2):
                if playable & row:
                    options.add(self.value(left, column, delete_left, row)
                                ^ self.value(column + 1, right, row, delete_right))
        return mex(options)

    def full_value(self):
        return self.value(0, len(self.columns))


def grid_core(cells, root):
    cells = sorted(set(cells))
    if root not in cells:
        raise ValueError("root must be playable")
    cells.remove(root)
    cells.insert(0, root)
    adjacency = tuple(sum(1 << j for j, u in enumerate(cells)
                          if abs(v[0] - u[0]) + abs(v[1] - u[1]) == 1)
                      for v in cells)
    return cells, adjacency


def full_head(width):
    return grid_core({(r, c) for c in range(1, width + 1) for r in (1, 2, 3)}
                     - {(1, width)}, (3, width))


def validate_core(adjacency):
    if not adjacency:
        raise ValueError("a rooted core must be nonempty")
    limit = (1 << len(adjacency)) - 1
    for v, neighbours in enumerate(adjacency):
        if type(neighbours) is not int or neighbours < 0 or neighbours & ~limit:
            raise ValueError("invalid adjacency bitmask")
        if neighbours >> v & 1:
            raise ValueError("loops are forbidden")
        for u in range(len(adjacency)):
            if bool(neighbours >> u & 1) != bool(adjacency[u] >> v & 1):
                raise ValueError("asymmetric adjacency")


class RootedComb:
    """Exact move-closed recurrence for a fixed core with a comb tail.

    Root is vertex 0. Safe core moves preserve it. All safe descendants are
    retained; there is no claim that the number of core states is constant
    as the core grows. Tail length alone does not cause exponential search.
    """

    def __init__(self, adjacency, phase=0):
        validate_core(adjacency)
        if phase not in (0, 1):
            raise ValueError("invalid phase")
        self.adjacency = tuple(adjacency)
        self.phase = phase
        self.closed = tuple(ns | (1 << v) for v, ns in enumerate(adjacency))
        self.full = (1 << len(adjacency)) - 1
        states = {self.full}
        pending = [self.full]
        for mask in pending:
            for v in range(1, len(adjacency)):
                if mask >> v & 1 and not self.adjacency[v] & 1:
                    successor = mask & ~self.closed[v]
                    if successor not in states:
                        states.add(successor)
                        pending.append(successor)
        self.states = sorted(states, key=lambda mask: (mask.bit_count(), mask))
        self.values = {mask: [self.core_value(mask)] for mask in self.states}

    @cache
    def core_value(self, mask):
        if not mask:
            return 0
        reached = frontier = mask & -mask
        while frontier:
            bit = frontier & -frontier
            frontier ^= bit
            new = self.adjacency[bit.bit_length() - 1] & mask & ~reached
            reached |= new
            frontier |= new
        if reached != mask:
            return self.core_value(reached) ^ self.core_value(mask ^ reached)
        return mex({self.core_value(mask & ~self.closed[v])
                    for v in range(len(self.adjacency)) if mask >> v & 1})

    def extend(self, max_length):
        p = self.phase
        for n in range(len(self.values[self.full]), max_length + 1):
            for mask in self.states:
                options = set()
                for v in range(len(self.adjacency)):
                    if not mask >> v & 1:
                        continue
                    rest = mask & ~self.closed[v]
                    if v == 0:
                        options.add(self.core_value(rest)
                                    ^ comb_value(p ^ 1, n - 1) ^ p)
                    elif self.adjacency[v] & 1:
                        options.add(self.core_value(rest) ^ comb_value(p, n))
                    else:
                        options.add(self.values[rest][n])
                # First backbone vertex deletes the root, but not its neighbours.
                options.add(self.core_value(mask & ~1)
                            ^ comb_value(p, max(n - 2, 0))
                            ^ ((1 - p) if n >= 2 else 0))
                for i in range(n):
                    tooth = p ^ (i % 2)
                    if tooth:
                        options.add(self.values[mask][i]
                                    ^ comb_value(p ^ ((i + 1) % 2), n - i - 1))
                    if i:
                        isolated = ((int(i > 0) + int(i + 1 < n))
                                    * (1 - tooth)) % 2
                        options.add(self.values[mask][i - 1]
                                    ^ comb_value(p ^ ((i + 2) % 2), max(n - i - 2, 0))
                                    ^ isolated)
                self.values[mask].append(mex(options))
        return self

    def period_start(self):
        available = len(self.values[self.full]) - 1
        for start in range(2, available - 8):
            if all(all(values[n] == values[n + 4]
                       for n in range(start, start + 6))
                   for values in self.values.values()):
                return start
        return None

    def periodic_options(self, mask, length):
        """All option VALUES, represented by finitely many actual moves.

        Requires the six-term lifting premises for every safe core state.
        Interior tail moves repeat after four columns; the last backbone
        vertex is kept separately because it leaves a different orphan leaf.
        """
        start = self.period_start()
        if start is None or mask not in self.values or length < 0:
            raise ValueError("a certified core state and nonnegative length are required")

        def lookup(state, n):
            index = n if n < start else start + (n - start) % 4
            return self.values[state][index]

        p = self.phase
        options = []
        for v in range(len(self.adjacency)):
            if not mask >> v & 1:
                continue
            rest = mask & ~self.closed[v]
            if not length:
                value = self.core_value(rest)
            elif v == 0:
                value = self.core_value(rest) ^ comb_value(p ^ 1, length - 1) ^ p
            elif self.adjacency[v] & 1:
                value = self.core_value(rest) ^ comb_value(p, length)
            else:
                value = lookup(rest, length)
            options.append((("core", v), value))
        if length:
            options.append((("spine", 0), self.core_value(mask & ~1)
                            ^ comb_value(p, max(length - 2, 0))
                            ^ ((1 - p) if length >= 2 else 0)))
            for i in range(min(length, start + 4)):
                if p ^ (i % 2):
                    options.append((("leaf", i), lookup(mask, i)
                                    ^ comb_value(0, length - i - 1)))
            indices = set(range(1, min(length - 1, start + 5)))
            if length >= 2:
                indices.add(length - 1)
            for i in sorted(indices):
                tooth = p ^ (i % 2)
                isolated = ((1 + int(i + 1 < length)) * (1 - tooth)) % 2
                options.append((("spine", i), lookup(mask, i - 1)
                                ^ comb_value(p ^ ((i + 2) % 2), max(length - i - 2, 0))
                                ^ isolated))
        return lookup(mask, length), options

    def certificate(self, name, cells, histories, core_width):
        start = self.period_start()
        if start is None:
            raise ValueError("no period-lifting certificate in computed prefix")
        return {
            "name": name, "phase": self.phase, "core_width": core_width,
            "cells": [list(v) for v in cells], "root": list(cells[0]),
            "edges": [[v, u] for v, ns in enumerate(self.adjacency)
                      for u in range(v + 1, len(self.adjacency)) if ns >> u & 1],
            "prior_history": [list(v) for v in histories],
            "tail_history": "top row at core_width+1, core_width+3, ...",
            "period": 4, "start": start, "window": 6,
            "states": [{"mask": f"{mask:x}",
                        "values": self.values[mask][:start + 10]}
                       for mask in self.states],
            "full_core_period_block": self.values[self.full][start:start + 4],
            "full_core_zero_residues": sorted({n % 4 for n in range(start, start + 4)
                                               if self.values[self.full][n] == 0}),
        }


def build_report(max_length=24):
    examples = []
    for width in (2, 3, 4, 6):
        cells, adjacency = full_head(width)
        model = RootedComb(adjacency).extend(max_length)
        examples.append(model.certificate(f"head-{width}", cells, [], width))
    for name, mask, history in (
        ("obstruction-H", 0x7ff4d37, [(1, 3), (1, 5)]),
        ("obstruction-J", 0x7ff4df4, [(1, 1), (1, 5)]),
    ):
        cells = {(v % 3 + 1, v // 3 + 1) for v in range(27) if mask >> v & 1}
        cells, adjacency = grid_core(cells - {(1, 9)}, (3, 9))
        model = RootedComb(adjacency).extend(max_length)
        examples.append(model.certificate(name, cells, history, 9))
    return {
        "format": "nodekayles-rooted-comb-period-v1",
        "theorem": "For every recorded root-preserving core state X, "
                   "SG(X with C_p(n+4)) = SG(X with C_p(n)) for all n >= start.",
        "scope": "Single-edge attachment, alternating comb tail. This does not "
                 "establish periodicity or outcomes of empty 3xn boards.",
        "examples": examples,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--max-tail", type=int, default=24)
    parser.add_argument("--head", type=int, help="query a full three-row head, phase 0")
    parser.add_argument("--tail", type=int, help="query an arbitrarily long alternating tail")
    args = parser.parse_args()
    if args.max_tail < 14:
        parser.error("at least 14 tail vertices are needed for the example windows")
    if args.head is not None or args.tail is not None:
        if args.head is None or args.tail is None or not 1 <= args.head <= 8 or args.tail < 1:
            parser.error("query requires --head in 1..8 and a positive --tail")
        cells, adjacency = full_head(args.head)
        model = RootedComb(adjacency).extend(args.max_tail)
        value, options = model.periodic_options(model.full, args.tail)
        winning = next((move for move, g in options if g == 0), None)
        if winning:
            kind, index = winning
            cell = (cells[index] if kind == "core"
                    else (3 if kind == "spine" else 2, args.head + index + 1))
        else:
            cell = None
        print(json.dumps({"head": args.head, "tail": args.tail, "SG": value,
                          "period_start": model.period_start(),
                          "representative_options": len(options),
                          "winning_move": cell}, ensure_ascii=False))
        if not args.output:
            return
    if not args.output:
        parser.error("provide --output for certificates or --head and --tail for a query")
    report = build_report(args.max_tail)
    # Keep each finite-state row on one line for a compact, reviewable certificate.
    raw = json.dumps(report, ensure_ascii=False, indent=2)
    import re
    raw = re.sub(r'\[\n\s+(-?\d+(?:,\n\s+-?\d+)*)\n\s+\]',
                 lambda match: '[' + ', '.join(re.findall(r'-?\d+', match[1])) + ']', raw)
    args.output.write_text(raw + "\n")
    for example in report["examples"]:
        print(example["name"], "states", len(example["states"]),
              "start", example["start"], "P residues", example["full_core_zero_residues"])


if __name__ == "__main__":
    main()
