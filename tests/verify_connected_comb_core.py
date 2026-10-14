#!/usr/bin/env python3
"""Independent verifier for component-normalized comb certificates.

Imports no production solver. Reconstructs every connected closure and checks
finite unary mex equations, followed by notes/23's complete refined spectrum.
The additional infinite lifting argument is documented in notes/25.
"""

from functools import cache
import gzip
import json
from pathlib import Path
import sys

from verify_two_comb_core import Checker, c


class ConnectedChecker(Checker):
    def __init__(self, proof):
        if proof.get("format") != "node-kayles-connected-two-comb-v1":
            raise ValueError("unknown connected-comb certificate format")
        self.graph = tuple(proof["adjacency"])
        k = len(self.graph)
        if (not k or any(type(ns) is not int or not 0 <= ns < 1 << k for ns in self.graph)
                or any(self.graph[v] >> v & 1 or bool(self.graph[v] >> u & 1) !=
                       bool(self.graph[u] >> v & 1) for v in range(k) for u in range(k))):
            raise ValueError("invalid underlying graph")
        self.roots, self.phases = tuple(proof["roots"]), tuple(proof["phases"])
        if (len(self.roots) != 2 or len(set(self.roots)) != 2
                or any(type(r) is not int or not 0 <= r < k for r in self.roots)
                or len(self.phases) != 2 or any(p not in (0, 1) for p in self.phases)):
            raise ValueError("invalid roots/phases")
        self.full, self.root_mask = (1 << k) - 1, sum(1 << root for root in self.roots)
        self.components = cache(self._components)
        self.core_sg = cache(self._core_sg)
        self.b, self.d, self.t = (proof[key] for key in ("minimum", "difference", "boundary"))
        if (any(type(x) is not int or x < 1 for x in (self.b, self.d, self.t)) or self.t < self.b):
            raise ValueError("invalid chamber thresholds")
        bound = max(self.t + 3, self.b + self.d + 6)
        if proof["sample_bound"] != bound:
            raise ValueError("invalid sample size")
        initial = self.joint(self.full)
        closure, pending = ({initial}, [initial]) if initial else (set(), [])
        for mask in pending:
            for v, ns in enumerate(self.graph):
                closed = ns | 1 << v
                if mask >> v & 1 and not closed & self.root_mask:
                    component = self.joint(mask & ~closed)
                    if component and component not in closure:
                        closure.add(component)
                        pending.append(component)
        self.rows = {}
        for item in proof["states"]:
            mask, rows = item["mask"], item["rows"]
            if (type(mask) is not int or mask not in closure or mask in self.rows
                    or len(rows) != bound + 1
                    or any(len(row) != bound + 1 or any(type(g) is not int or g < 0 for g in row)
                           for row in rows)):
                raise ValueError("invalid connected double-root table")
            self.rows[mask] = rows
        if set(self.rows) != closure:
            raise ValueError("missing connected double-root state")
        for mask, rows in self.rows.items():
            if any(g != self.candidate(mask, n, m)
                   for n, row in enumerate(rows) for m, g in enumerate(row)):
                raise ValueError("sample disagrees with its chamber formula")

        seeds = [set(), set()]
        for mask in [self.full] + list(closure):
            successors = [mask]
            successors.extend(mask & ~(ns | 1 << v) for v, ns in enumerate(self.graph)
                              if mask >> v & 1)
            successors.extend(mask & ~(1 << root) for root in self.roots)
            for rest in successors:
                for part in self.components(rest):
                    for side, root in enumerate(self.roots):
                        if part >> root & 1:
                            seeds[side].add(part)
        if len(proof["unary_banks"]) != 2:
            raise ValueError("two unary banks required")
        self.unary = {}
        self.checked_unary_values = self.unary_states = 0
        for side, table in enumerate(proof["unary_banks"]):
            self.unary[side] = self.verify_bank(table, side, seeds[side])
        self.unique_unary_tables = 2
        self.unary_bound = max(table["start"] for table in proof["unary_banks"])
        if self.t < self.unary_bound:
            raise ValueError("boundary does not cover unary periods")
        self.C = max(self.t, self.b + self.d, self.unary_bound) + 8
        self.D = self.d + 8
        self.U = self.C + self.D + 8
        self.R = self.U + 3

    def _components(self, mask):
        # Set-based graph search, independent of production bit-frontier code.
        unseen = {v for v in range(len(self.graph)) if mask >> v & 1}
        result = []
        while unseen:
            seed = min(unseen)
            reached, pending = {seed}, [seed]
            unseen.remove(seed)
            for v in pending:
                added = {u for u in unseen if self.graph[v] >> u & 1}
                unseen -= added
                reached |= added
                pending.extend(added)
            result.append(sum(1 << v for v in reached))
        return tuple(result)

    def _core_sg(self, mask):
        parts = self.components(mask)
        if len(parts) != 1:
            result = 0
            for part in parts:
                result ^= self.core_sg(part)
            return result
        options = {self.core_sg(mask & ~(ns | 1 << v))
                   for v, ns in enumerate(self.graph) if mask >> v & 1}
        result = 0
        while result in options:
            result += 1
        return result

    def joint(self, mask):
        return next((part for part in self.components(mask)
                     if part & self.root_mask == self.root_mask), 0)

    def verify_bank(self, table, side, seeds):
        root, phase = self.roots[side], self.phases[side]
        start = table["start"]
        if (table["root"] != root or table["phase"] != phase or type(start) is not int or start < 2
                or any(type(mask) is not int for mask in table["seeds"])
                or len(table["seeds"]) != len(seeds) or set(table["seeds"]) != seeds):
            raise ValueError("invalid unary seed set, root, phase or start")
        closure, pending = set(seeds), list(seeds)
        for mask in pending:
            for v, ns in enumerate(self.graph):
                closed = ns | 1 << v
                if mask >> v & 1 and not closed >> root & 1:
                    rest = mask & ~closed
                    component = next(part for part in self.components(rest) if part >> root & 1)
                    if component not in closure:
                        closure.add(component)
                        pending.append(component)
        rows = {}
        for item in table["states"]:
            mask, values = item["mask"], item["values"]
            if (type(mask) is not int or mask not in closure or mask in rows
                    or len(values) != start + 10 or any(type(g) is not int or g < 0 for g in values)):
                raise ValueError("invalid connected unary state")
            rows[mask] = values
            if any(values[n] != values[n + 4] for n in range(start, start + 6)):
                raise ValueError("unary period window fails")
        if set(rows) != closure:
            raise ValueError("missing connected unary successor")
        for mask, row in rows.items():
            if row[0] != self.core_sg(mask):
                raise ValueError("unary zero-length SG fails")
            self.checked_unary_values += 1
            for n in range(1, start + 10):
                options = set()
                for v, ns in enumerate(self.graph):
                    if not mask >> v & 1:
                        continue
                    rest = mask & ~(ns | 1 << v)
                    if v == root:
                        value = self.core_sg(rest) ^ c(phase ^ 1, n - 1) ^ phase
                    elif ns >> root & 1:
                        value = self.core_sg(rest) ^ c(phase, n)
                    else:
                        component = next(part for part in self.components(rest) if part >> root & 1)
                        value = rows[component][n] ^ self.core_sg(rest ^ component)
                    options.add(value)
                options.add(self.core_sg(mask & ~(1 << root)) ^ c(phase, max(n - 2, 0))
                            ^ ((1 - phase) if n >= 2 else 0))
                # Expand ALL tail moves here; no production residue histories.
                for i in range(n):
                    tooth = phase ^ (i % 2)
                    if tooth:
                        options.add(row[i] ^ c(0, n - i - 1))
                    if i:
                        options.add(row[i - 1] ^ c(phase ^ ((i + 2) % 2), max(n - i - 2, 0))
                                    ^ ((1 + int(i + 1 < n)) * (1 - tooth) % 2))
                value = 0
                while value in options:
                    value += 1
                if value != row[n]:
                    raise ValueError(f"unary mex fails at mask {mask}, length {n}")
                self.checked_unary_values += 1
        self.unary_states += len(rows)
        return start, rows

    def terminal(self, mask, lengths):
        result = 0
        for root, phase, n in zip(self.roots, self.phases, lengths):
            if n and not mask >> root & 1:
                result ^= c(phase, n)
        for part in self.components(mask):
            active = [side for side, root in enumerate(self.roots)
                      if lengths[side] and part >> root & 1]
            if len(active) > 1:
                raise ValueError("two active roots in a terminal component")
            if active:
                side = active[0]
                start, rows = self.unary[side]
                n = lengths[side]
                result ^= rows[part][n if n < start else start + (n - start) % 4]
            else:
                result ^= self.core_sg(part)
        return result

    def at(self, mask, n, m):
        joint = self.joint(mask) if n and m else 0
        if joint:
            return self.candidate(joint, n, m) ^ self.core_sg(mask ^ joint)
        return self.terminal(mask, (n, m))

    def verify(self):
        result = super().verify()
        # These are checked finite equations, not raw whole-graph unary searches.
        result["finite_unary_mex_identities"] = result.pop("original_unary_values")
        result["connected_unary_states"] = self.unary_states
        result["format"] = "node-kayles-connected-two-comb-audit-v1"
        return result


if __name__ == "__main__":
    path = Path(sys.argv[1])
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as source:
        print(json.dumps(ConnectedChecker(json.load(source)).verify()))
