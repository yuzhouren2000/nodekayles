#!/usr/bin/env python3
"""Independent finite premises for a two-comb semilinear SG certificate.

No production solver is imported. All unary tables are checked by the original
named-vertex game. The infinite step is the bounded option-spectrum lemma in
notes/23: all coordinate and diagonal chambers have a bounded representative.
"""

from functools import cache
import json
import gzip
from pathlib import Path
import sys

from verify_comb_bridge import original_game


def c(phase, n):
    return ((0, 1, 2, 1), (0, 1, 2, 3))[phase][n % 4]


def remap(graph, mask, root=None):
    vertices = ([] if root is None else [root]) + [v for v in range(len(graph))
                 if mask >> v & 1 and v != root]
    return tuple(sum(1 << j for j, u in enumerate(vertices) if graph[v] >> u & 1)
                 for v in vertices)


def verify_unary_table(table):
    """Original vertex-subset game; share its cache across recorded subsets."""
    graph = table["adjacency"]
    k = len(graph)
    start, phase = table["start"], table["phase"]
    if (not k or type(start) is not int or start < 2 or phase not in (0, 1)
            or any(type(ns) is not int or not 0 <= ns < 1 << k for ns in graph)
            or any(graph[v] >> v & 1 or bool(graph[v] >> u & 1) != bool(graph[u] >> v & 1)
                   for v in range(k) for u in range(k))):
        raise ValueError("invalid unary graph, phase or start")
    limit = start + 9
    edges = {("core", v): set() for v in range(k)}
    def edge(v, u):
        edges.setdefault(v, set()).add(u)
        edges.setdefault(u, set()).add(v)
    for v in range(k):
        for u in range(v):
            if graph[v] >> u & 1:
                edge(("core", v), ("core", u))
    for i in range(limit):
        edge(("spine", i), ("core", 0) if not i else ("spine", i - 1))
        if phase ^ (i % 2):
            edge(("spine", i), ("leaf", i))
    @cache
    def sg(position):
        if not position:
            return 0
        seen = {next(iter(position))}
        pending = list(seen)
        while pending:
            v = pending.pop()
            new = (edges[v] & position) - seen
            seen.update(new)
            pending.extend(new)
        if len(seen) != len(position):
            return sg(frozenset(seen)) ^ sg(position - seen)
        options = {sg(position - edges[v] - {v}) for v in position}
        value = 0
        while value in options:
            value += 1
        return value
    full = frozenset(("core", v) for v in range(k))
    closure, pending = {full}, [full]
    for position in pending:
        for v in position:
            if v == ("core", 0) or ("core", 0) in edges[v]:
                continue
            rest = frozenset(position - edges[v] - {v})
            if rest not in closure:
                closure.add(rest)
                pending.append(rest)
    tails = [frozenset()]
    for i in range(limit):
        added = {("spine", i)}
        if phase ^ (i % 2):
            added.add(("leaf", i))
        tails.append(tails[-1] | added)
    checked, recorded = 0, set()
    for row in table["states"]:
        mask, values = row["mask"], row["values"]
        position = frozenset(("core", v) for v in range(k) if mask >> v & 1)
        if (type(mask) is not int or not 0 < mask < 1 << k or position not in closure
                or position in recorded or len(values) != start + 10
                or any(type(value) is not int or value < 0 for value in values)):
            raise ValueError("invalid unary state or value row")
        recorded.add(position)
        if any(values[n] != values[n + 4] for n in range(start, start + 6)):
            raise ValueError("unary period window fails")
        for n, expected in enumerate(values):
            if sg(position | tails[n]) != expected:
                raise ValueError("unary original vertex-subset SG fails")
            checked += 1
    if recorded != closure:
        raise ValueError("missing unary root-preserving state")
    return checked


class Checker:
    def __init__(self, proof, check_unary=True):
        if proof.get("format") != "node-kayles-two-comb-semilinear-v1":
            raise ValueError("unknown two-comb certificate format")
        self.graph = tuple(proof["adjacency"])
        k = len(self.graph)
        if (not k or any(type(ns) is not int or not 0 <= ns < 1 << k for ns in self.graph)
                or any(self.graph[v] >> v & 1 or bool(self.graph[v] >> u & 1) !=
                       bool(self.graph[u] >> v & 1) for v in range(k) for u in range(k))):
            raise ValueError("invalid core adjacency")
        self.roots, self.phases = tuple(proof["roots"]), tuple(proof["phases"])
        if (len(self.roots) != 2 or len(set(self.roots)) != 2
                or any(type(r) is not int or not 0 <= r < k for r in self.roots)
                or len(self.phases) != 2 or any(p not in (0, 1) for p in self.phases)):
            raise ValueError("invalid roots or phases")
        self.full = (1 << k) - 1
        self.root_mask = sum(1 << r for r in self.roots)
        self.b, self.d, self.t = (proof[key] for key in ("minimum", "difference", "boundary"))
        if any(type(x) is not int or x < 1 for x in (self.b, self.d, self.t)) or self.t < self.b:
            raise ValueError("invalid chamber thresholds")
        bound = max(self.t + 3, self.b + self.d + 6)
        if proof["sample_bound"] != bound:
            raise ValueError("wrong sample dimensions")
        self.rows = {}
        closure, pending = {self.full}, [self.full]
        for mask in pending:
            for v, ns in enumerate(self.graph):
                closed = ns | 1 << v
                if mask >> v & 1 and not closed & self.root_mask:
                    rest = mask & ~closed
                    if rest not in closure:
                        closure.add(rest)
                        pending.append(rest)
        for item in proof["states"]:
            mask, rows = item["mask"], item["rows"]
            if (mask not in closure or mask in self.rows or len(rows) != bound + 1
                    or any(len(row) != bound + 1 or any(type(v) is not int or v < 0 for v in row)
                           for row in rows)):
                raise ValueError("incorrect, missing, extra or duplicate core table")
            self.rows[mask] = rows
        if set(self.rows) != closure:
            raise ValueError("missing a root-preserving core state")
        for mask, rows in self.rows.items():
            if any(value != self.candidate(mask, n, m)
                   for n, row in enumerate(rows) for m, value in enumerate(row)):
                raise ValueError("stored sample disagrees with its chamber formula")
        self.unary = {}
        self.checked_unary_values = 0
        checked_tables = {}
        for table in proof["unary"]:
            key = table["core_mask"], table["side"]
            if key in self.unary or key[1] not in (0, 1):
                raise ValueError("duplicate or invalid unary table")
            mask, side = key
            if not 0 < mask <= self.full or not mask >> self.roots[side] & 1:
                raise ValueError("unary root has been deleted")
            if (tuple(table["adjacency"]) != remap(self.graph, mask, self.roots[side])
                    or table["phase"] != self.phases[side]):
                raise ValueError("unary table has the wrong underlying core")
            if check_unary:
                table_key = (tuple(table["adjacency"]), table["phase"], table["start"],
                             tuple((row["mask"], tuple(row["values"])) for row in table["states"]))
                if table_key not in checked_tables:
                    checked_tables[table_key] = verify_unary_table(table)
                self.checked_unary_values += checked_tables[table_key]
            full = (1 << len(table["adjacency"])) - 1
            values = next(item["values"] for item in table["states"] if item["mask"] == full)
            self.unary[key] = table["start"], values
        self.unary_bound = max((start for start, _ in self.unary.values()), default=2)
        if self.t < self.unary_bound:
            raise ValueError("short-axis threshold does not cover unary periods")
        # Eight exceeds the period (4), the largest prefix offset (3),
        # and their sum. Thus each unbounded residue interval is saturated.
        self.C = max(self.t, self.b + self.d, self.unary_bound) + 8
        self.D = self.d + 8
        self.U = self.C + self.D + 8
        self.R = self.U + 3
        self.unique_unary_tables = len(checked_tables)

    def canonical(self, n, m):
        if min(n, m) < self.b:
            return (n if n < self.t else self.t + (n - self.t) % 4,
                    m if m < self.t else self.t + (m - self.t) % 4)
        if abs(n - m) < self.d:
            shift = 4 * ((min(n, m) - self.b) // 4)
            return n - shift, m - shift
        if n <= m:
            a = self.b + (n - self.b) % 4
            return a, a + self.d + (m - n - self.d) % 4
        a = self.b + (m - self.b) % 4
        return a + self.d + (n - m - self.d) % 4, a

    def candidate(self, mask, n, m):
        a, b = self.canonical(n, m)
        return self.rows[mask][a][b]

    @cache
    def core_sg(self, mask):
        return original_game(remap(self.graph, mask), None, 0, 0)[0]

    def terminal(self, mask, lengths):
        offset, active = 0, []
        for side, (root, p, length) in enumerate(zip(self.roots, self.phases, lengths)):
            if length and not mask >> root & 1:
                offset ^= c(p, length)
            elif length:
                active.append(side)
        if not active:
            return offset ^ self.core_sg(mask)
        if len(active) != 1:
            raise ValueError("terminal query has two active arms")
        side = active[0]
        if (mask, side) not in self.unary:
            raise ValueError("missing required unary table")
        start, values = self.unary[mask, side]
        n = lengths[side]
        return offset ^ values[n if n < start else start + (n - start) % 4]

    def at(self, mask, n, m):
        if not n or not m or mask & self.root_mask != self.root_mask:
            return self.terminal(mask, (n, m))
        return self.candidate(mask, n, m)

    def prefix_positions(self, upper, other):
        """All candidate values with residue: finite pieces of an integer interval."""
        if upper < 0:
            return []
        # The short cells and the narrow diagonal strip have exact coordinates.
        # Far intervals depend only on residue, so at most four representatives.
        intervals = []
        exact = set()
        if other < self.b:
            exact.update(range(min(upper + 1, self.t)))
            intervals.append((self.t, upper))
        else:
            exact.update(range(min(upper + 1, self.b)))
            low, high = max(self.b, other - self.d + 1), min(upper, other + self.d - 1)
            if low <= high:
                exact.update(range(low, high + 1))
            intervals.extend(((self.b, min(upper, other - self.d)),
                              (max(self.b, other + self.d), upper)))
        for low, high in intervals:
            if low <= high:
                exact.update(range(low, min(high + 1, low + 4)))
        return sorted(exact)

    def options(self, mask, n, m):
        lengths = n, m
        options = set()
        for v, ns in enumerate(self.graph):
            if not mask >> v & 1:
                continue
            rest, new, offset = mask & ~(ns | 1 << v), [n, m], 0
            for side, (root, p, length) in enumerate(zip(self.roots, self.phases, lengths)):
                if length and v == root:
                    new[side] = 0
                    offset ^= c(p ^ 1, length - 1) ^ p
            options.add(self.at(rest, *new) ^ offset)
        for side, (root, p, length, other) in enumerate(zip(self.roots, self.phases, lengths,
                                                           lengths[::-1])):
            if not length:
                continue
            new = [n, m]
            new[side] = 0
            options.add(self.terminal(mask & ~(1 << root), tuple(new))
                        ^ c(p, max(length - 2, 0)) ^ ((1 - p) if length >= 2 else 0))
            # For old prefix j <= length-3 use the proven coefficient identity.
            if length >= 3:
                for j in self.prefix_positions(length - 3, other):
                    a, b = (j, other) if side == 0 else (other, j)
                    prefix = self.at(mask, a, b)
                    options.add(prefix ^ c(p ^ ((j + 3) % 2), length - j - 3))
                    if p ^ (j % 2):
                        options.add(prefix ^ c(0, length - j - 1))
                a, b = (length - 2, other) if side == 0 else (other, length - 2)
                delta = p ^ (length % 2)
                options.add(self.at(mask, a, b) ^ delta)
                if not delta:
                    a, b = (length - 1, other) if side == 0 else (other, length - 1)
                    options.add(self.at(mask, a, b))
            else:
                for i in range(length):
                    tooth = p ^ (i % 2)
                    if tooth:
                        a, b = (i, other) if side == 0 else (other, i)
                        options.add(self.at(mask, a, b) ^ c(0, length - i - 1))
                    if i:
                        a, b = (i - 1, other) if side == 0 else (other, i - 1)
                        options.add(self.at(mask, a, b) ^ c(p ^ ((i + 2) % 2),
                                                           max(length - i - 2, 0))
                                    ^ ((1 + int(i + 1 < length)) * (1 - tooth) % 2))
        return options

    def refined(self, n, m):
        if min(n, m) < self.C:
            return (n if n < self.U else self.U + (n - self.U) % 4,
                    m if m < self.U else self.U + (m - self.U) % 4)
        if abs(n - m) < self.D:
            shift = 4 * ((min(n, m) - self.C) // 4)
            return n - shift, m - shift
        if n <= m:
            a = self.C + (n - self.C) % 4
            return a, a + self.D + (m - n - self.D) % 4
        a = self.C + (m - self.C) % 4
        return a + self.D + (n - m - self.D) % 4, a

    def verify(self):
        positions = {(n, m) for n in range(self.R + 1) for m in range(self.R + 1)
                     if self.refined(n, m) == (n, m)}
        checked = 0
        for mask in self.rows:
            for n, m in positions:
                expected = self.candidate(mask, n, m)
                if not n or not m:
                    actual = self.terminal(mask, (n, m))
                else:
                    options = self.options(mask, n, m)
                    actual = 0
                    while actual in options:
                        actual += 1
                if expected != actual:
                    raise ValueError(f"mex identity fails: core mask {mask}, lengths {n},{m}, "
                                     f"claimed {expected}, actual {actual}")
                checked += 1
        return {"core_states": len(self.rows), "unary_tables": len(self.unary),
                "distinct_unary_tables_checked": self.unique_unary_tables,
                "original_unary_values": self.checked_unary_values,
                "finite_mex_identities": checked, "representative_pairs": len(positions),
                "proof_rectangle": self.R, "thresholds": [self.b, self.d, self.t]}


def verify(proof):
    return Checker(proof).verify()


if __name__ == "__main__":
    path = Path(sys.argv[1])
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as source:
        print(json.dumps(verify(json.load(source))))
