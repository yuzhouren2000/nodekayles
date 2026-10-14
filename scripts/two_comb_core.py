#!/usr/bin/env python3
"""Exact two-arm comb recurrence and finite diagonal-period candidates.

Two distinct roots of an arbitrary core each receive one alternating comb.
The central core may stay connected after its own moves. No independent
period reduction is assumed for the two lengths. See notes/23 for the finite
semilinear certificate criterion.
"""

import argparse
from functools import cache
import gzip
import json
from pathlib import Path

from kayles_boundary_algebra import grundy, induced
from kayles_boundary_transfer import column_graph
from strip3_defect_algebra import RootedComb, comb_value, validate_core


def tail_graph(core, roots, phases, lengths):
    graph = list(core)
    for root, phase, length in zip(roots, phases, lengths):
        if type(length) is not int or length < 0:
            raise ValueError("arm lengths must be nonnegative integers")
        if not length:
            continue
        first = len(graph)
        graph += [0] * length
        def edge(v, u):
            graph[v] |= 1 << u
            graph[u] |= 1 << v
        edge(root, first)
        for i in range(length):
            if i:
                edge(first + i - 1, first + i)
            if phase ^ (i % 2):
                graph.append(0)
                edge(first + i, len(graph) - 1)
    return tuple(graph)


def reorder_root(graph, mask, root):
    vertices = [root] + [v for v in range(len(graph)) if mask >> v & 1 and v != root]
    return tuple(sum(1 << j for j, u in enumerate(vertices) if graph[v] >> u & 1)
                 for v in vertices)


class TwoCombCore:
    def __init__(self, adjacency, roots, phases=(1, 1)):
        validate_core(adjacency)
        self.graph = tuple(adjacency)
        self.roots, self.phases = tuple(roots), tuple(phases)
        if (len(self.roots) != 2 or len(set(self.roots)) != 2
                or any(type(r) is not int or not 0 <= r < len(self.graph) for r in self.roots)
                or len(self.phases) != 2 or any(p not in (0, 1) for p in self.phases)):
            raise ValueError("two distinct valid roots and two binary phases required")
        self.full = (1 << len(self.graph)) - 1
        self.closed = tuple(ns | 1 << v for v, ns in enumerate(self.graph))
        self.root_mask = sum(1 << r for r in self.roots)
        closure, pending = {self.full}, [self.full]
        for mask in pending:
            for v in range(len(self.graph)):
                if mask >> v & 1 and not self.closed[v] & self.root_mask:
                    rest = mask & ~self.closed[v]
                    if rest not in closure:
                        closure.add(rest)
                        pending.append(rest)
        self.states = sorted(closure, key=lambda mask: (mask.bit_count(), mask))
        self.unary = {}
        self.rows = None
        self.maximum = -1
        self.core_value = cache(self.core_value)
        self.value = cache(self._value)

    def core_value(self, mask):
        return grundy(induced(self.graph, mask))

    def unary_model(self, mask, side):
        key = mask, side
        if key not in self.unary:
            self.unary[key] = RootedComb(reorder_root(self.graph, mask, self.roots[side]),
                                        self.phases[side]).stabilize()
        return self.unary[key]

    def terminated(self, mask, lengths):
        """Handle at most one surviving nonempty attached arm exactly."""
        offset = 0
        active = []
        for side, (root, phase, n) in enumerate(zip(self.roots, self.phases, lengths)):
            if n and not mask >> root & 1:
                offset ^= comb_value(phase, n)
            elif n:
                active.append(side)
        if not active:
            return self.core_value(mask) ^ offset
        if len(active) > 1:
            raise ValueError("both arms are active")
        side = active[0]
        model = self.unary_model(mask, side)
        start = model.period_start()
        n = lengths[side]
        j = n if n < start else start + (n - start) % 4
        return model.values[model.full][j] ^ offset

    def direct_options(self, mask, lengths, lookup=None, terminated=None):
        """Core moves and the first main-chain move of either arm."""
        lookup = self.value if lookup is None else lookup
        terminated = self.terminated if terminated is None else terminated
        options = []
        for v in range(len(self.graph)):
            if not mask >> v & 1:
                continue
            rest, new, offset = mask & ~self.closed[v], list(lengths), 0
            for side, (root, p, n) in enumerate(zip(self.roots, self.phases, lengths)):
                if n and v == root:
                    new[side] = 0
                    offset ^= comb_value(p ^ 1, n - 1) ^ p
            value = (lookup(rest, tuple(new)) if rest & self.root_mask == self.root_mask
                     and all(new) else terminated(rest, tuple(new)))
            options.append((("core", v), value ^ offset))
        for side, (root, p, n) in enumerate(zip(self.roots, self.phases, lengths)):
            if not n:
                continue
            new = list(lengths)
            new[side] = 0
            value = terminated(mask & ~(1 << root), tuple(new))
            value ^= comb_value(p, max(n - 2, 0)) ^ ((1 - p) if n >= 2 else 0)
            options.append((("spine", side, 0), value))
        return options

    def literal_options(self, mask, lengths, lookup=None):
        lookup = self.value if lookup is None else lookup
        options = self.direct_options(mask, lengths, lookup)
        for side, (p, n) in enumerate(zip(self.phases, lengths)):
            for i in range(n):
                tooth = p ^ (i % 2)
                if tooth:
                    new = list(lengths)
                    new[side] = i
                    options.append((("leaf", side, i), lookup(mask, tuple(new))
                                    ^ comb_value(0, n - i - 1)))
                if i:
                    new = list(lengths)
                    new[side] = i - 1
                    isolated = (1 + int(i + 1 < n)) * (1 - tooth) % 2
                    options.append((("spine", side, i), lookup(mask, tuple(new))
                                    ^ comb_value(p ^ ((i + 2) % 2), max(n - i - 2, 0))
                                    ^ isolated))
        return options

    def _value(self, mask, lengths):
        if (len(lengths) != 2 or any(type(n) is not int or n < 0 for n in lengths)):
            raise ValueError("two nonnegative lengths required")
        if mask & self.root_mask != self.root_mask or not all(lengths):
            return self.terminated(mask, lengths)
        if self.rows is not None and max(lengths) <= self.maximum and mask in self.rows:
            return self.rows[mask][lengths[0]][lengths[1]]
        values = {g for _, g in self.literal_options(mask, lengths)}
        result = 0
        while result in values:
            result += 1
        return result

    def build(self, maximum):
        """Finite exact sample, with four residue-history sets per direction."""
        if type(maximum) is not int or maximum < 0:
            raise ValueError("nonnegative sample bound required")
        self.rows = {mask: [[0] * (maximum + 1) for _ in range(maximum + 1)]
                     for mask in self.states}
        self.maximum = maximum
        for mask in self.states:
            rows = self.rows[mask]
            for n in range(maximum + 1):
                rows[n][0] = self.terminated(mask, (n, 0))
                rows[0][n] = self.terminated(mask, (0, n))
            column_history = [[set() for _ in range(4)] for _ in range(maximum + 1)]
            for n in range(1, maximum + 1):
                row_history = [set() for _ in range(4)]
                for m in range(1, maximum + 1):
                    values = {g for _, g in self.direct_options(mask, (n, m), self._sample)}
                    for side, (p, length, other, histories) in enumerate((
                            (self.phases[0], n, m, column_history[m]),
                            (self.phases[1], m, n, row_history))):
                        if length < 3:
                            for i in range(length):
                                tooth = p ^ (i % 2)
                                if tooth:
                                    prefix = rows[i][other] if side == 0 else rows[other][i]
                                    values.add(prefix ^ comb_value(0, length - i - 1))
                                if i:
                                    prefix = rows[i - 1][other] if side == 0 else rows[other][i - 1]
                                    values.add(prefix ^ comb_value(p ^ ((i + 2) % 2),
                                                                  max(length - i - 2, 0))
                                               ^ ((1 + int(i + 1 < length)) * (1 - tooth) % 2))
                        else:
                            prefix = rows[length - 3][other] if side == 0 else rows[other][length - 3]
                            histories[(length - 3) % 4].add(prefix)
                            for residue, seen in enumerate(histories):
                                offsets = {comb_value(p ^ ((residue + 3) % 2),
                                                      (length - residue - 3) % 4)}
                                if p ^ (residue % 2):
                                    offsets.add(comb_value(0, (length - residue - 1) % 4))
                                values.update(g ^ d for g in seen for d in offsets)
                            prefix = rows[length - 2][other] if side == 0 else rows[other][length - 2]
                            delta = p ^ (length % 2)
                            values.add(prefix ^ delta)
                            if not delta:
                                prefix = rows[length - 1][other] if side == 0 else rows[other][length - 1]
                                values.add(prefix)
                    result = 0
                    while result in values:
                        result += 1
                    rows[n][m] = result
        self.value.cache_clear()
        return self

    def _sample(self, mask, lengths):
        return self.rows[mask][lengths[0]][lengths[1]]

    def candidate(self, minimum, difference, boundary):
        return DiagonalCandidate(self, minimum, difference, boundary)


class DiagonalCandidate:
    def __init__(self, model, minimum, difference, boundary):
        if min(minimum, difference, boundary) < 1 or boundary < minimum:
            raise ValueError("positive thresholds and boundary >= minimum required")
        self.model = model
        self.minimum, self.difference, self.boundary = minimum, difference, boundary
        self.sample_bound = max(boundary + 3, minimum + difference + 6)
        if self.sample_bound > model.maximum:
            raise ValueError("exact sample is too short")

    def canonical(self, lengths):
        n, m = lengths
        b, d, t = self.minimum, self.difference, self.boundary
        if min(n, m) < b:
            return (n if n < t else t + (n - t) % 4,
                    m if m < t else t + (m - t) % 4)
        if abs(m - n) < d:
            shift = ((min(n, m) - b) // 4) * 4
            return n - shift, m - shift
        if m >= n:
            a = b + (n - b) % 4
            return a, a + d + (m - n - d) % 4
        a = b + (m - b) % 4
        return a + d + (n - m - d) % 4, a

    def value(self, mask, lengths):
        n, m = self.canonical(lengths)
        return self.model.rows[mask][n][m]

    def sample_agrees(self):
        return all(g == self.value(mask, (n, m))
                   for mask, rows in self.model.rows.items()
                   for n, row in enumerate(rows) for m, g in enumerate(row))

    def prefix_positions(self, upper, other):
        """Actual prefix lengths covering each value and residue in 0..upper."""
        if upper < 0:
            return []
        b, d, t = self.minimum, self.difference, self.boundary
        positions = set()
        if other < b:
            positions.update(range(min(upper + 1, t)))
            intervals = ((t, upper),)
        else:
            positions.update(range(min(upper + 1, b)))
            low, high = max(b, other - d + 1), min(upper, other + d - 1)
            if low <= high:
                positions.update(range(low, high + 1))
            intervals = ((b, min(upper, other - d)),
                         (max(b, other + d), upper))
        for low, high in intervals:
            if low <= high:
                positions.update(range(low, min(high + 1, low + 4)))
        return sorted(positions)

    def proof(self):
        unary = []
        for (mask, side), model in sorted(self.model.unary.items()):
            start = model.period_start()
            unary.append({"core_mask": mask, "side": side, "start": start,
                          "adjacency": model.adjacency, "phase": model.phase,
                          "states": [{"mask": submask,
                                      "values": model.values[submask][:start + 10]}
                                     for submask in model.states]})
        k = self.sample_bound
        return {"format": "node-kayles-two-comb-semilinear-v1", "adjacency": self.model.graph,
                "roots": self.model.roots, "phases": self.model.phases,
                "minimum": self.minimum, "difference": self.difference,
                "boundary": self.boundary, "sample_bound": k,
                "states": [{"mask": mask, "rows": [row[:k + 1] for row in rows[:k + 1]]}
                           for mask, rows in self.model.rows.items()],
                "unary": unary}


class CertifiedTwoComb:
    """SG and actual representative moves from a previously verified proof.

    This is a query reader, not an independent certificate verifier. Its
    working graph is the finite core; neither long arm is materialized.
    """

    def __init__(self, proof):
        if proof.get("format") != "node-kayles-two-comb-semilinear-v1":
            raise ValueError("unknown two-comb certificate format")
        self.model = TwoCombCore(proof["adjacency"], proof["roots"], proof["phases"])
        self.model.maximum = proof["sample_bound"]
        self.model.rows = {row["mask"]: row["rows"] for row in proof["states"]}
        self.candidate = self.model.candidate(proof["minimum"], proof["difference"],
                                              proof["boundary"])
        self.unary = {}
        for table in proof["unary"]:
            full = (1 << len(table["adjacency"])) - 1
            values = next(row["values"] for row in table["states"] if row["mask"] == full)
            self.unary[table["core_mask"], table["side"]] = table["start"], values

    def terminated(self, mask, lengths):
        offset, active = 0, []
        for side, (root, phase, n) in enumerate(zip(self.model.roots, self.model.phases, lengths)):
            if n and not mask >> root & 1:
                offset ^= comb_value(phase, n)
            elif n:
                active.append(side)
        if not active:
            return self.model.core_value(mask) ^ offset
        if len(active) != 1:
            raise ValueError("terminal query has two active arms")
        side = active[0]
        start, values = self.unary[mask, side]
        n = lengths[side]
        return offset ^ values[n if n < start else start + (n - start) % 4]

    def value(self, lengths, mask=None):
        if len(lengths) != 2 or any(type(n) is not int or n < 0 for n in lengths):
            raise ValueError("two nonnegative lengths required")
        mask = self.model.full if mask is None else mask
        if not all(lengths) or mask & self.model.root_mask != self.model.root_mask:
            return self.terminated(mask, lengths)
        return self.candidate.value(mask, lengths)

    def representative_options(self, lengths, mask=None):
        """Named original moves covering the entire option-value set."""
        mask = self.model.full if mask is None else mask
        self.value(lengths, mask)  # Validate lengths before constructing moves.
        if mask & self.model.root_mask != self.model.root_mask:
            raise ValueError("representative moves require both core roots")
        lookup = lambda rest, pair: self.value(pair, rest)
        # The same exact core/first-spine recurrence uses stored unary tables.
        # No discovery solver is needed while querying the certificate.
        options = self.model.direct_options(mask, lengths, lookup, self.terminated)
        for side, (p, length, other) in enumerate(zip(self.model.phases, lengths, lengths[::-1])):
            def prefix(j):
                pair = (j, other) if side == 0 else (other, j)
                return self.value(pair, mask)
            if length < 3:
                for i in range(length):
                    tooth = p ^ (i % 2)
                    if tooth:
                        options.append((("leaf", side, i), prefix(i)
                                        ^ comb_value(0, length - i - 1)))
                    if i:
                        options.append((("spine", side, i), prefix(i - 1)
                                        ^ comb_value(p ^ ((i + 2) % 2), max(length - i - 2, 0))
                                        ^ ((1 + int(i + 1 < length)) * (1 - tooth) % 2)))
                continue
            for j in self.candidate.prefix_positions(length - 3, other):
                value = prefix(j)
                options.append((("spine", side, j + 1), value
                                ^ comb_value(p ^ ((j + 3) % 2), length - j - 3)))
                if p ^ (j % 2):
                    options.append((("leaf", side, j), value ^ comb_value(0, length - j - 1)))
            delta = p ^ (length % 2)
            options.append((("spine", side, length - 1), prefix(length - 2) ^ delta))
            if not delta:
                options.append((("leaf", side, length - 1), prefix(length - 1)))
        return options

    def winning_move(self, lengths):
        if self.value(lengths) == 0:
            return None
        for move, value in self.representative_options(lengths):
            if value == 0:
                return move
        raise ValueError("certificate is inconsistent: nonzero SG has no zero-valued option")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--certificate", type=Path, help="query a previously verified table")
    parser.add_argument("--lengths", nargs=2, type=int, metavar=("LEFT", "RIGHT"))
    parser.add_argument("--columns", default="67776")
    parser.add_argument("--roots", nargs=2, type=int)
    parser.add_argument("--phases", nargs=2, type=int, choices=(0, 1), default=(0, 0))
    parser.add_argument("--sample", type=int, default=64)
    parser.add_argument("--minimum", type=int, default=8)
    parser.add_argument("--difference", type=int, default=12)
    parser.add_argument("--boundary", type=int, default=24)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.certificate:
        if args.lengths is None or min(args.lengths) < 0:
            parser.error("--certificate requires two nonnegative --lengths")
        opener = gzip.open if args.certificate.suffix == ".gz" else open
        with opener(args.certificate, "rt") as source:
            proof = json.load(source)
        # This reads the proved table; the independent verifier is a separate
        # command so arbitrary-length queries do not repeat its finite proof.
        query = CertifiedTwoComb(proof)
        lengths = tuple(args.lengths)
        print(json.dumps({"lengths": lengths,
                          "canonical_lengths": query.candidate.canonical(lengths),
                          "SG": query.value(lengths), "winning_move": query.winning_move(lengths),
                          "move_indexing": "zero-based; arm side 0=left, 1=right; index from core",
                          "thresholds": [proof[k] for k in ("minimum", "difference", "boundary")]}))
        return
    if not args.columns or any(c not in "01234567" for c in args.columns):
        parser.error("column word must use digits 0..7")
    graph, cells = column_graph(tuple(map(int, args.columns)))
    default_roots = (cells.index((2, 0)) if int(args.columns[0]) & 4 else 0,
                     cells.index((2, len(args.columns) - 1)) if int(args.columns[-1]) & 4
                     else len(graph) - 1)
    model = TwoCombCore(graph, args.roots or default_roots, args.phases).build(args.sample)
    candidate = model.candidate(args.minimum, args.difference, args.boundary)
    agrees = candidate.sample_agrees()
    print(json.dumps({"core": args.columns, "states": len(model.states),
                      "sample": args.sample, "unary_tables": len(model.unary),
                      "sample_agrees": agrees,
                      "thresholds": [args.minimum, args.difference, args.boundary]}))
    if args.output:
        output = (json.dumps(candidate.proof(), separators=(",", ":")) + "\n").encode()
        args.output.write_bytes(gzip.compress(output, mtime=0) if args.output.suffix == ".gz"
                                else output)


if __name__ == "__main__":
    main()
