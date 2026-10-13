#!/usr/bin/env python3
"""Check bridge certificates by named-vertex Node Kayles, without solvers."""

from functools import cache
import json
from pathlib import Path
import sys


def original_game(left, right, length, phase, left_mask=None, right_mask=None):
    """Build edges and every genuine move with no interface reduction."""
    edges = {}

    def vertex(v):
        edges.setdefault(v, set())

    def edge(v, u):
        vertex(v)
        vertex(u)
        edges[v].add(u)
        edges[u].add(v)

    for side, core, mask in (("left", left, left_mask), ("right", right, right_mask)):
        if core is None:
            continue
        if mask is None:
            mask = (1 << len(core)) - 1
        for v, ns in enumerate(core):
            if not mask >> v & 1:
                continue
            vertex((side, v))
            for u in range(v):
                if mask >> u & 1 and ns >> u & 1:
                    edge((side, v), (side, u))
    if length:
        edge(("left", 0), ("spine", 0))
        if right is not None:
            edge(("right", 0), ("spine", length - 1))
        for i in range(length):
            vertex(("spine", i))
            if i:
                edge(("spine", i - 1), ("spine", i))
            if phase ^ (i % 2):
                edge(("spine", i), ("leaf", i))
    elif right is not None:
        edge(("left", 0), ("right", 0))

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
        outcomes = {sg(position - edges[v] - {v}) for v in position}
        result = 0
        while result in outcomes:
            result += 1
        return result

    initial = frozenset(edges)
    options = {v: sg(initial - edges[v] - {v}) for v in initial}
    return sg(initial), options


def verify_table(table):
    graph = table["adjacency"]
    k = len(graph)
    if not k or any(type(ns) is not int or ns < 0 or ns >= 1 << k for ns in graph):
        raise ValueError("invalid adjacency")
    if any(graph[v] >> v & 1 or bool(graph[v] >> u & 1) != bool(graph[u] >> v & 1)
           for v in range(k) for u in range(k)):
        raise ValueError("loops or asymmetric adjacency")
    phase, start = table["phase"], table["start"]
    if phase not in (0, 1) or type(start) is not int or start < 2:
        raise ValueError("invalid phase or period start")
    full = (1 << k) - 1
    closure, pending = {full}, [full]
    for mask in pending:
        for v in range(1, k):
            if mask >> v & 1 and not graph[v] & 1:
                rest = mask & ~(graph[v] | 1 << v)
                if rest not in closure:
                    closure.add(rest)
                    pending.append(rest)
    rows = {}
    checked = 0
    for row in table["states"]:
        mask, values = row["mask"], row["values"]
        if mask not in closure or mask in rows or len(values) != start + 10:
            raise ValueError("missing, extra or duplicate state or wrong table length")
        if any(type(value) is not int or value < 0 for value in values):
            raise ValueError("invalid SG entry")
        for n, expected in enumerate(values):
            actual, _ = original_game(graph, None, n, phase, mask)
            if actual != expected:
                raise ValueError("incorrect original-game SG entry")
            checked += 1
        if any(values[n] != values[n + 4] for n in range(start, start + 6)):
            raise ValueError("one-ended lifting window fails")
        rows[mask] = values
    if set(rows) != closure:
        raise ValueError("incomplete root-preserving closure")
    return checked


def verify(proof):
    if proof.get("format") != "comb-bridge-period-v1":
        raise ValueError("unknown certificate format")
    tables = proof["one_ended_tables"]
    if len(tables) != 4 or {(t["side"], t["phase"]) for t in tables} != {
            (s, p) for s in ("left", "right") for p in (0, 1)}:
        raise ValueError("both cores and both phases are required")
    checked = sum(verify_table(table) for table in tables)
    by_side = {s: [t for t in tables if t["side"] == s] for s in ("left", "right")}
    if any(pair[0]["adjacency"] != pair[1]["adjacency"] for pair in by_side.values()):
        raise ValueError("phase tables use different cores")
    starts = {s: max(t["start"] for t in pair) for s, pair in by_side.items()}
    onset = starts["left"] + starts["right"] + 6
    if (proof["left_start"] != starts["left"] or proof["right_start"] != starts["right"]
            or proof["period_start"] != onset):
        raise ValueError("wrong bridge threshold")
    left = by_side["left"][0]["adjacency"]
    right = by_side["right"][0]["adjacency"]
    bound = sum(192 * len(core) ** 2 + 984 * len(core) + 3 for core in (left, right)) + 6
    if proof["order_bound"] != bound or onset > bound:
        raise ValueError("wrong order bound")
    for phase in (0, 1):
        values = proof["period_values"][str(phase)]
        if len(values) != 4:
            raise ValueError("wrong bridge period block")
        for i, value in enumerate(values):
            actual, _ = original_game(left, right, onset + i, phase)
            if actual != value:
                raise ValueError("incorrect bridge period value")
            checked += 1
    return {"one_ended_tables": 4, "original_game_values": checked,
            "bridge_period_start": onset, "period_values": proof["period_values"]}


if __name__ == "__main__":
    print(json.dumps(verify(json.loads(Path(sys.argv[1]).read_text()))))
