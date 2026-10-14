#!/usr/bin/env python3
"""Component-normalized one/two-root comb certificates (see notes/25).

This is an exact representation change, using disjoint-component SG xor.
Detached fragments are never silently discarded. The old raw-mask proof
format and verifier are unchanged.
"""

import argparse
from functools import cache
import gzip
import json
from pathlib import Path

from strip3_defect_algebra import comb_value, mex, stabilization_bound, validate_core
from two_comb_core import CertifiedTwoComb, DiagonalCandidate, TwoCombCore
from kayles_boundary_transfer import column_graph


class ComponentGame:
    def __init__(self, graph):
        validate_core(graph)
        self.graph = tuple(graph)
        self.full = (1 << len(graph)) - 1
        self.closed = tuple(ns | 1 << v for v, ns in enumerate(graph))
        # Per-instance caches do not retain every previously constructed graph.
        self.parts = cache(self._parts)
        self.sg = cache(self._sg)

    def _parts(self, mask):
        components = []
        remaining = mask
        while remaining:
            reached = pending = remaining & -remaining
            while pending:
                bit = pending & -pending
                pending ^= bit
                added = self.graph[bit.bit_length() - 1] & mask & ~reached
                reached |= added
                pending |= added
            components.append(reached)
            remaining ^= reached
        return tuple(components)

    def _sg(self, mask):
        components = self.parts(mask)
        if len(components) != 1:
            value = 0
            for part in components:
                value ^= self.sg(part)
            return value
        return mex({self.sg(mask & ~self.closed[v])
                    for v in range(len(self.graph)) if mask >> v & 1})

    def rooted(self, mask, root):
        """The root component and the SG of everything detached from it."""
        component = next((part for part in self.parts(mask) if part >> root & 1), 0)
        return component, self.sg(mask ^ component)


class UnaryBank:
    """Shared connected-root closure of a finite set of ambient-graph seeds."""

    def __init__(self, game, root, phase, seeds):
        self.game, self.root, self.phase = game, root, phase
        self.seeds = sorted(set(seeds))
        closure, pending = set(self.seeds), list(self.seeds)
        for mask in pending:
            if (not mask >> root & 1 or game.parts(mask) != (mask,)):
                raise ValueError("unary seeds must be connected and contain their root")
            for v, closed in enumerate(game.closed):
                if mask >> v & 1 and not closed >> root & 1:
                    rest, _ = game.rooted(mask & ~closed, root)
                    if rest not in closure:
                        closure.add(rest)
                        pending.append(rest)
        self.states = sorted(closure, key=lambda mask: (mask.bit_count(), mask))
        self.values = {mask: [game.sg(mask)] for mask in self.states}
        self.histories = {mask: [set() for _ in range(4)] for mask in self.states}
        self.start = None

    def extend(self, maximum):
        if not self.states:
            return self
        p, root, game = self.phase, self.root, self.game
        for n in range(len(self.values[self.states[0]]), maximum + 1):
            for mask in self.states:
                options = set()
                for v, closed in enumerate(game.closed):
                    if not mask >> v & 1:
                        continue
                    rest = mask & ~closed
                    if v == root:
                        options.add(game.sg(rest) ^ comb_value(p ^ 1, n - 1) ^ p)
                    elif closed >> root & 1:
                        options.add(game.sg(rest) ^ comb_value(p, n))
                    else:
                        component, offset = game.rooted(rest, root)
                        options.add(self.values[component][n] ^ offset)
                options.add(game.sg(mask & ~(1 << root))
                            ^ comb_value(p, max(n - 2, 0))
                            ^ ((1 - p) if n >= 2 else 0))
                if n < 3:
                    for i in range(n):
                        tooth = p ^ (i % 2)
                        if tooth:
                            options.add(self.values[mask][i] ^ comb_value(0, n - i - 1))
                        if i:
                            options.add(self.values[mask][i - 1]
                                        ^ comb_value(p ^ ((i + 2) % 2), max(n - i - 2, 0))
                                        ^ ((1 + int(i + 1 < n)) * (1 - tooth) % 2))
                else:
                    history = self.histories[mask]
                    history[(n - 3) % 4].add(self.values[mask][n - 3])
                    for residue, seen in enumerate(history):
                        offsets = {comb_value(p ^ ((residue + 3) % 2), (n - residue - 3) % 4)}
                        if p ^ (residue % 2):
                            offsets.add(comb_value(0, (n - residue - 1) % 4))
                        options.update(g ^ d for g in seen for d in offsets)
                    delta = p ^ (n % 2)
                    options.add(self.values[mask][n - 2] ^ delta)
                    if not delta:
                        options.add(self.values[mask][n - 1])
                self.values[mask].append(mex(options))
        return self

    def period_start(self):
        if self.start is not None:
            return self.start
        if not self.states:
            return 2
        available = len(self.values[self.states[0]]) - 1
        for start in range(2, available - 8):
            if all(all(row[n] == row[n + 4] for n in range(start, start + 6))
                   for row in self.values.values()):
                self.start = start
                return start
        return None

    def stabilize(self):
        limit = stabilization_bound(max((mask.bit_count() for mask in self.states), default=1)) + 9
        target = min(12, limit)
        while True:
            self.extend(target)
            if self.period_start() is not None:
                return self
            if target == limit:
                raise AssertionError("normalized comb stabilization bound violated")
            target = min(limit, target * 2)

    def value(self, mask, length):
        component, offset = self.game.rooted(mask, self.root)
        if not component:
            return offset ^ comb_value(self.phase, length)
        start = self.period_start()
        if start is None:
            raise ValueError("unary bank has not stabilized")
        j = length if length < start else start + (length - start) % 4
        return self.values[component][j] ^ offset

    def proof(self):
        start = self.period_start()
        if start is None:
            raise ValueError("unary bank has not stabilized")
        return {"root": self.root, "phase": self.phase, "seeds": self.seeds, "start": start,
                "states": [{"mask": mask, "values": self.values[mask][:start + 10]}
                           for mask in self.states]}


class ConnectedTwoCombCore(TwoCombCore):
    """Only retain connected components containing both active roots."""

    def __init__(self, adjacency, roots, phases=(0, 0)):
        self.game = ComponentGame(adjacency)
        self.graph, self.full, self.closed = self.game.graph, self.game.full, self.game.closed
        self.roots, self.phases = tuple(roots), tuple(phases)
        if (len(self.roots) != 2 or len(set(self.roots)) != 2
                or any(type(r) is not int or not 0 <= r < len(self.graph) for r in self.roots)
                or len(self.phases) != 2 or any(p not in (0, 1) for p in self.phases)):
            raise ValueError("two distinct valid roots and two binary phases required")
        self.root_mask = sum(1 << r for r in self.roots)
        initial = self.joint_component(self.full)
        closure, pending = ({initial}, [initial]) if initial else (set(), [])
        for mask in pending:
            for v, closed in enumerate(self.closed):
                if mask >> v & 1 and not closed & self.root_mask:
                    rest = self.joint_component(mask & ~closed)
                    if rest and rest not in closure:
                        closure.add(rest)
                        pending.append(rest)
        self.states = sorted(closure, key=lambda mask: (mask.bit_count(), mask))
        # Includes axes and every exit caused by a core/first-spine move.
        seeds = [set(), set()]
        for mask in [self.full] + self.states:
            successors = [mask] + [mask & ~closed for v, closed in enumerate(self.closed)
                                   if mask >> v & 1]
            successors += [mask & ~(1 << root) for root in self.roots]
            for rest in successors:
                for side, root in enumerate(self.roots):
                    component = next((part for part in self.game.parts(rest) if part >> root & 1), 0)
                    if component:
                        seeds[side].add(component)
        self.banks = [UnaryBank(self.game, root, p, seed).stabilize()
                      for root, p, seed in zip(self.roots, self.phases, seeds)]
        self.rows, self.maximum = None, -1
        self.core_value = self.game.sg
        self.value = cache(self._value)

    def joint_component(self, mask):
        return next((part for part in self.game.parts(mask)
                     if part & self.root_mask == self.root_mask), 0)

    def terminated(self, mask, lengths):
        """Exact xor of zero/one-root components; also handles separated roots."""
        result = 0
        for side, (root, p, n) in enumerate(zip(self.roots, self.phases, lengths)):
            if n and not mask >> root & 1:
                result ^= comb_value(p, n)
        for component in self.game.parts(mask):
            active = [side for side, root in enumerate(self.roots)
                      if lengths[side] and component >> root & 1]
            if len(active) == 2:
                raise ValueError("terminal component has two active roots")
            if active:
                side = active[0]
                result ^= self.banks[side].value(component, lengths[side])
            else:
                result ^= self.game.sg(component)
        return result

    def _sample(self, mask, lengths):
        joint = self.joint_component(mask) if all(lengths) else 0
        if not joint:
            return self.terminated(mask, lengths)
        return self.rows[joint][lengths[0]][lengths[1]] ^ self.game.sg(mask ^ joint)

    def _value(self, mask, lengths):
        if len(lengths) != 2 or any(type(n) is not int or n < 0 for n in lengths):
            raise ValueError("two nonnegative lengths required")
        joint = self.joint_component(mask) if all(lengths) else 0
        if not joint:
            return self.terminated(mask, lengths)
        if self.rows is not None and max(lengths) <= self.maximum:
            return self._sample(mask, lengths)
        return mex({g for _, g in self.literal_options(joint, lengths)}) ^ self.game.sg(mask ^ joint)

    def candidate(self, minimum, difference, boundary):
        return ConnectedCandidate(self, minimum, difference, boundary)


class ConnectedCandidate(DiagonalCandidate):
    def value(self, mask, lengths):
        joint = self.model.joint_component(mask) if all(lengths) else 0
        if not joint:
            return self.model.terminated(mask, lengths)
        n, m = self.canonical(lengths)
        return self.model.rows[joint][n][m] ^ self.model.game.sg(mask ^ joint)

    def proof(self):
        k = self.sample_bound
        return {"format": "node-kayles-connected-two-comb-v1", "adjacency": self.model.graph,
                "roots": self.model.roots, "phases": self.model.phases,
                "minimum": self.minimum, "difference": self.difference,
                "boundary": self.boundary, "sample_bound": k,
                "states": [{"mask": mask, "rows": [row[:k + 1] for row in rows[:k + 1]]}
                           for mask, rows in self.model.rows.items()],
                "unary_banks": [bank.proof() for bank in self.model.banks]}


class CertifiedConnectedTwoComb(CertifiedTwoComb):
    """Read a VERIFIED connected certificate without rebuilding its closures."""

    def __init__(self, proof):
        if proof.get("format") != "node-kayles-connected-two-comb-v1":
            raise ValueError("unknown connected-comb certificate format")
        model = ConnectedTwoCombCore.__new__(ConnectedTwoCombCore)
        model.game = ComponentGame(proof["adjacency"])
        model.graph, model.full, model.closed = model.game.graph, model.game.full, model.game.closed
        model.roots, model.phases = tuple(proof["roots"]), tuple(proof["phases"])
        model.root_mask = sum(1 << root for root in model.roots)
        model.states = [item["mask"] for item in proof["states"]]
        model.maximum = proof["sample_bound"]
        model.rows = {item["mask"]: item["rows"] for item in proof["states"]}
        model.core_value = model.game.sg
        model.banks = []
        for table in proof["unary_banks"]:
            bank = UnaryBank.__new__(UnaryBank)
            bank.game, bank.root, bank.phase = model.game, table["root"], table["phase"]
            bank.start, bank.seeds = table["start"], table["seeds"]
            bank.states = [item["mask"] for item in table["states"]]
            bank.values = {item["mask"]: item["values"] for item in table["states"]}
            model.banks.append(bank)
        self.model = model
        self.candidate = model.candidate(proof["minimum"], proof["difference"], proof["boundary"])

    def terminated(self, mask, lengths):
        return self.model.terminated(mask, lengths)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--columns", required=True)
    parser.add_argument("--sample", type=int, default=64)
    parser.add_argument("--minimum", type=int, default=8)
    parser.add_argument("--difference", type=int, default=12)
    parser.add_argument("--boundary", type=int, default=24)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if not args.columns or any(c not in "01234567" for c in args.columns):
        parser.error("columns must use digits 0..7")
    graph, cells = column_graph(tuple(map(int, args.columns)))
    roots = cells.index((2, 0)), cells.index((2, len(args.columns) - 1))
    model = ConnectedTwoCombCore(graph, roots).build(args.sample)
    candidate = model.candidate(args.minimum, args.difference, args.boundary)
    print(json.dumps({"columns": args.columns, "connected_states": len(model.states),
                      "unary_states": [len(bank.states) for bank in model.banks],
                      "unary_starts": [bank.period_start() for bank in model.banks],
                      "sample_agrees": candidate.sample_agrees()}), flush=True)
    if args.output:
        output = (json.dumps(candidate.proof(), separators=(",", ":")) + "\n").encode()
        args.output.write_bytes(gzip.compress(output, mtime=0) if args.output.suffix == ".gz" else output)


if __name__ == "__main__":
    main()
