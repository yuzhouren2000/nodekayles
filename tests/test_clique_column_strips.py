#!/usr/bin/env python3
"""Graph-level checks of interval closure, period use, and the ladder formula."""

import copy
from functools import cache
from itertools import combinations, product
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from periodic_clique_strip import (CliqueColumnIntervals, DefectiveThreeRowIntervals, PeriodicStrip,
                                  RootedPeriodicStrip, ladder_endpoint_value)
from strip3_defect_algebra import full_head
from verify_clique_strip_period import check, cells, after_history
from verify_strip3_two_move import recompute


def raw_solver(adjacency):
    @cache
    def sg(mask):
        if not mask:
            return 0
        reached = frontier = mask & -mask
        while frontier:
            bit = frontier & -frontier
            frontier ^= bit
            new = adjacency[bit.bit_length() - 1] & mask & ~reached
            reached |= new
            frontier |= new
        if reached != mask:
            return sg(reached) ^ sg(mask ^ reached)
        forbidden = {sg(mask & ~(neighbours | (1 << i)))
                     for i, neighbours in enumerate(adjacency) if mask >> i & 1}
        g = 0
        while g in forbidden:
            g += 1
        return g
    return sg


def graph(position):
    vertices = sorted(position)
    adjacency = tuple(sum(1 << j for j, u in enumerate(vertices)
                          if abs(v[0] - u[0]) + abs(v[1] - u[1]) == 1)
                      for v in vertices)
    return vertices, adjacency


def test_arbitrary_words():
    checked = 0
    for width in range(5):
        vertices, adjacency = graph({(r, c) for r in (1, 2, 3) for c in range(1, width + 1)})
        raw = raw_solver(adjacency)
        for word in product((0, 1, 2, 3, 4, 6), repeat=width):
            intervals = CliqueColumnIntervals(word)
            for left, right in product(range(8), repeat=2):
                position = cells(word, width, left=left, right=right) if width else set()
                mask = sum(1 << i for i, v in enumerate(vertices) if v in position)
                assert intervals.value(0, width, left, right) == raw(mask)
                checked += 1
    print(f"arbitrary three-row clique-column words: {checked} independent SG checks")
    for invalid in ((), (5,), (7,), (-1,), (True,), (8,)):
        try:
            PeriodicStrip(invalid)
        except ValueError:
            pass
        else:
            raise AssertionError("non-clique column accepted")


def test_ladder():
    # A proof-oriented complete option-spectrum check, including endpoints.
    checked = 0
    for n in range(1, 501):
        options = {ladder_endpoint_value(i, 0, bit) ^
                   ladder_endpoint_value(n - i - 1, bit, 1)
                   for i in range(n) for bit in (1, 2)
                   if not (i == n - 1 and bit == 1)}
        assert options == set(range(n)), (n, options)
        assert ladder_endpoint_value(n, 0, 1) == n
        checked += 1
    for n in range(13):
        vertices, adjacency = graph({(r, c) for r in (1, 2) for c in range(1, n + 1)})
        raw = raw_solver(adjacency)
        for left, right in product(range(4), repeat=2):
            position = cells((3,), n, left=left, right=right)
            mask = sum(1 << i for i, v in enumerate(vertices) if v in position)
            assert ladder_endpoint_value(n, left, right) == raw(mask)
            checked += 1
    print(f"unbounded ladder defect: {checked} SG or complete spectrum checks")
    attached = RootedPeriodicStrip((0,), PeriodicStrip((3,)), attachment=1).extend(100)
    assert all(attached.value(1, n) == n + 1 for n in range(101))
    assert attached.certificate(2) is None


def test_snake_history():
    for k in range(31):
        n = 3 * k + 1
        position = {(r, c) for r in (1, 2, 3) for c in range(1, n + 1)}
        for j in range(k + 1):
            move = (1 if j % 2 == 0 else 3, 3 * j + 1)
            assert move in position
            position = {v for v in position if abs(v[0] - move[0]) + abs(v[1] - move[1]) > 1}
        _, adjacency = graph(position)
        assert len(adjacency) == 5 * k + 1
        assert all(a.bit_count() <= 2 for a in adjacency)
        assert sum(a.bit_count() for a in adjacency) == 2 * (len(adjacency) - 1)
        reached = frontier = 1
        while frontier:
            bit = frontier & -frontier
            frontier ^= bit
            new = adjacency[bit.bit_length() - 1] & ~reached
            reached |= new
            frontier |= new
        assert reached == (1 << len(adjacency)) - 1
    print("31 alternating outer-row histories produce connected snake paths")


def test_defective_words():
    checked = 0
    for width in range(6):
        vertices, adjacency = graph({(r, c) for r in (1, 2, 3) for c in range(1, width + 1)})
        raw = raw_solver(adjacency)
        for word in product(range(8), repeat=width):
            model = DefectiveThreeRowIntervals(word)
            position = cells(word, width) if width else set()
            mask = sum(1 << i for i, v in enumerate(vertices) if v in position)
            assert model.full_value() == raw(mask), word
            checked += 1
    # Non-clique mask 5 can arise at an endpoint by deleting row 2 in an
    # originally full column. Explicitly check arbitrary endpoint deletion
    # states and all tracked special choices, not just the initial graph.
    width = 4
    vertices, adjacency = graph({(r, c) for r in (1, 2, 3) for c in range(1, width + 1)})
    raw = raw_solver(adjacency)
    for word in ((7, 3, 6, 7), (5, 7, 5, 7), (7, 7, 7, 7)):
        model = DefectiveThreeRowIntervals(word)
        for choices in range(3 ** len(model.special)):
            transformed = model.transformed(choices)
            for left, right in product(range(8), repeat=2):
                position = cells(transformed, width, left=left, right=right)
                mask = sum(1 << i for i, v in enumerate(vertices) if v in position)
                assert model.value(choices, 0, width, left, right) == raw(mask)
                checked += 1
    print(f"all induced 3x5 graphs and tracked defect/endpoint states: {checked} independent SG checks")


def test_periodic_options():
    checked = 0
    report = json.loads((ROOT / "data/proofs/clique_strip_period.json").read_text())
    for item in report["strips"]:
        strip = PeriodicStrip(item["word"]).extend(130)
        start, period = item["start"], item["period"]
        assert strip.certificate(period) == start
        for n in (0, 1, 2, 11, 35, 70, 119, 130):
            for phase in range(strip.degree):
                for left, right in product(range(8), repeat=2):
                    actual = set(g for _, g in strip.options(n, phase, left, right))
                    representatives = strip.options(n, phase, left, right, (start, period))
                    assert actual == {g for _, g in representatives}
                    assert strip.value(n, phase, left, right) == strip.lookup(
                        n, phase, left, right, start, period)
                    checked += 1
    for item in report["cores"]:
        head = item["head"]
        vertices, adjacency = full_head(head)
        strip = PeriodicStrip(item["word"])
        core = RootedPeriodicStrip(adjacency, strip).extend(130)
        start, period = item["start"], item["period"]
        assert core.certificate(period) == start
        for mask in core.core.states:
            for n in (0, 1, 2, 11, 35, 70, 119, 130):
                for right in range(8):
                    actual = {g for _, g in core.options(mask, n, right)}
                    representatives = core.options(mask, n, right, (start, period))
                    assert actual == {g for _, g in representatives}
                    assert core.value(mask, n, right) == core.lookup(mask, n, right, start, period)
                    checked += 1
        # Compare claimed real moves against direct graph deletion, not just
        # against the recurrence used to select representatives.
        for n in (1, 2, 7, 11, 12, 17):
            position = set(vertices) | cells(item["word"], n, offset=head)
            ordered, adj = graph(position)
            raw = raw_solver(adj)
            full = (1 << len(ordered)) - 1
            assert raw(full) == core.value(core.core.full, n)
            for move, value in core.options(core.core.full, n, certificate=(start, period)):
                vertex = vertices[move[1]] if move[0] == "core" else (move[1], head + move[2])
                i = ordered.index(vertex)
                assert raw(full & ~(adj[i] | (1 << i))) == value
                checked += 1
        for n in range(1, 101):
            if n % 10 not in (0, 5):
                assert after_history(head, n) == set(vertices) | cells(item["word"], n, offset=head)
    print(f"period lookup, complete representative spectra, actual moves: {checked} checks")


def test_small_arbitrary_cores():
    checked = 0
    for order in range(1, 4):
        possible = list(combinations(range(order), 2))
        for edge_mask in range(1 << len(possible)):
            adjacency = [0] * order
            for i, (v, u) in enumerate(possible):
                if edge_mask >> i & 1:
                    adjacency[v] |= 1 << u
                    adjacency[u] |= 1 << v
            for word in ((4, 6, 3, 1), (4, 6, 4, 6, 3, 1, 3, 1, 3, 6), (3,)):
                core = RootedPeriodicStrip(adjacency, PeriodicStrip(word)).extend(20)
                for n in range(21):
                    tail_vertices, tail_adj = graph(cells(word, n))
                    attached = list(adjacency) + [a << order for a in tail_adj]
                    if n:
                        i = tail_vertices.index((core.first.bit_length(), 1))
                        attached[0] |= 1 << (order + i)
                        attached[order + i] |= 1
                    raw = raw_solver(attached)
                    for state in core.core.states:
                        for right in range(8):
                            mask = state | sum(1 << (order + i) for i, (r, c) in enumerate(tail_vertices)
                                               if c != n or not (right >> (r - 1) & 1))
                            assert raw(mask) == core.value(state, n, right)
                            actual = {raw(mask & ~(a | (1 << i))) for i, a in enumerate(attached)
                                      if mask >> i & 1}
                            assert actual == {g for _, g in core.options(state, n, right)}
                            checked += 1
    print(f"all simple cores of order <=3, three tails: {checked} raw graph/spectrum checks")


def test_bad_certificates():
    report = json.loads((ROOT / "data/proofs/clique_strip_period.json").read_text())
    mutations = [
        lambda x: x["strips"][0]["word"].__setitem__(0, 7),
        lambda x: x["strips"][0].__setitem__("period", 3),
        lambda x: x["strips"][0]["rows"].pop(),
        lambda x: x["strips"][0]["rows"][3].pop(),
        lambda x: x["cores"][0].__setitem__("start", 1),
        lambda x: x["cores"][1]["states"].pop(),
        lambda x: x["cores"][1]["adjacency"].__setitem__(0, 0),
        lambda x: x["cores"][1].__setitem__("reachable_length_residues", list(range(10))),
        lambda x: x["cores"][1]["full_period_block"].__setitem__(0, 99),
    ]
    for mutate in mutations:
        changed = copy.deepcopy(report)
        mutate(changed)
        try:
            check(changed)
        except (ValueError, KeyError):
            pass
        else:
            raise AssertionError("altered certificate accepted")
    premises, _ = check(report)
    false = dict(next(iter(premises.values())))
    false["value"] += 1
    try:
        recompute([false], 20)
    except ValueError:
        pass
    else:
        raise AssertionError("false finite SG premise accepted")
    print("10 invalid certificates or SG premises rejected")


if __name__ == "__main__":
    test_arbitrary_words()
    test_ladder()
    test_snake_history()
    test_defective_words()
    test_periodic_options()
    test_small_arbitrary_cores()
    test_bad_certificates()
