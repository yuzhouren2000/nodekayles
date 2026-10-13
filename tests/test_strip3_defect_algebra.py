#!/usr/bin/env python3
"""Independent graph checks for interval closure, comb values and period use."""

import copy
from functools import cache
from itertools import combinations, product
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from strip3_defect_algebra import RootedComb, TwoRowIntervals, comb_value, stabilization_bound
from verify_strip3_defect_period import check
from verify_strip3_two_move import recompute


def reference(adjacency):
    """Only raw graph deletion and disjoint components; no strip formulas."""
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
        options = set()
        for v, neighbours in enumerate(adjacency):
            if mask >> v & 1:
                options.add(sg(mask & ~((1 << v) | neighbours)))
        g = 0
        while g in options:
            g += 1
        return g
    return sg


def ladder(width):
    cells = [(r, c) for c in range(width) for r in (0, 1)]
    return tuple(sum(1 << j for j, u in enumerate(cells)
                     if abs(u[0] - v[0]) + abs(u[1] - v[1]) == 1)
                 for v in cells)


def attach(adjacency, length, phase):
    k = len(adjacency)
    result = list(adjacency) + [0] * (length + sum(phase ^ (i % 2) for i in range(length)))

    def edge(v, u):
        result[v] |= 1 << u
        result[u] |= 1 << v

    if length and k:
        edge(0, k)
    leaf = k + length
    masks = [(1 << k) - 1]
    for i in range(length):
        if i:
            edge(k + i - 1, k + i)
        mask = masks[-1] | (1 << (k + i))
        if phase ^ (i % 2):
            edge(k + i, leaf)
            mask |= 1 << leaf
            leaf += 1
        masks.append(mask)
    return tuple(result), masks


def test_intervals():
    checked = 0
    for width in range(6):
        sg = reference(ladder(width))
        for word in product(range(4), repeat=width):
            model = TwoRowIntervals(word)
            for left_deleted, right_deleted in product(range(4), repeat=2):
                mask = sum(col << (2 * i) for i, col in enumerate(word))
                if width:
                    mask &= ~left_deleted
                    mask &= ~(right_deleted << (2 * (width - 1)))
                assert model.value(0, width, left_deleted, right_deleted) == sg(mask)
                checked += 1
    print(f"arbitrary two-row words and endpoint deletions: {checked} independent SG checks")


def test_combs_and_contexts():
    for phase in (0, 1):
        adjacency, masks = attach((), 24, phase)
        sg = reference(adjacency)
        for length, mask in enumerate(masks):
            assert sg(mask) == comb_value(phase, length)
    checked = 0
    contexts = 0
    for order in range(1, 5):
        possible = list(combinations(range(order), 2))
        for edge_mask in range(1 << len(possible)):
            adjacency = [0] * order
            for i, (v, u) in enumerate(possible):
                if edge_mask >> i & 1:
                    adjacency[v] |= 1 << u
                    adjacency[u] |= 1 << v
            for phase in (0, 1):
                # Incremental history updates and unconditional discovery are
                # checked against raw graph deletion, including full spectra.
                model = RootedComb(adjacency, phase).extend(5).extend(8).stabilize().extend(24)
                graph, masks = attach(adjacency, 24, phase)
                sg = reference(graph)
                assert model.period_start() is not None
                assert model.period_start() <= stabilization_bound(order)
                for state in model.states:
                    assert max(model.values[state]) < 4 * (state.bit_count() + 2)
                    for length, full in enumerate(masks):
                        position = full & ~((1 << order) - 1) | state
                        assert model.values[state][length] == sg(position)
                        value, options = model.periodic_options(state, length)
                        actual_options = {sg(position & ~((1 << v) | neighbours))
                                          for v, neighbours in enumerate(graph) if position >> v & 1}
                        assert value == sg(position)
                        assert {g for _, g in options} == actual_options
                        # Each representative is an actual legal move with its
                        # claimed value, including the separately kept endpoint.
                        for (kind, i), expected in options:
                            if kind == "core":
                                vertex = i
                            elif kind == "spine":
                                vertex = order + i
                            else:
                                vertex = order + 24 + sum(phase ^ (j % 2) for j in range(i))
                            assert position >> vertex & 1
                            assert sg(position & ~((1 << vertex) | graph[vertex])) == expected
                        checked += 1
                contexts += 1
    print(f"all rooted simple graphs of order <=4, both phases: {contexts} contexts, "
          f"{checked} independent SG and complete option-spectrum checks")


def test_mutations():
    report = json.loads((ROOT / "data/proofs/strip3_defect_period.json").read_text())
    small = {**report, "examples": [report["examples"][0]]}
    mutations = [
        lambda x: x["examples"][0].__setitem__("phase", 1),
        lambda x: x["examples"][0].__setitem__("window", 5),
        lambda x: x["examples"][0].__setitem__("period", 8),
        lambda x: x["examples"][0]["states"].pop(),
        lambda x: x["examples"][0]["states"][0]["values"].__setitem__(2, 99),
        lambda x: x["examples"][0].__setitem__("root", [2, 2]),
        lambda x: x["examples"][0]["edges"].pop(),
        lambda x: x["examples"][0].__setitem__("prior_history", [[1, 1], [1, 1]]),
        lambda x: x["examples"][0].__setitem__("full_core_zero_residues", [0]),
    ]
    for mutate in mutations:
        changed = copy.deepcopy(small)
        mutate(changed)
        try:
            check(changed)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid certificate accepted")
    # A changed early SG value can leave the period window intact. The separate
    # recomputation layer must reject it rather than trusting consistent rows.
    changed = copy.deepcopy(small)
    changed["examples"][0]["states"][0]["values"][0] += 1
    premises, _ = check(changed)
    try:
        recompute(premises, 60)
    except ValueError:
        pass
    else:
        raise AssertionError("false finite SG premise accepted")
    print("10 altered certificates rejected, including a false SG outside the period window")


if __name__ == "__main__":
    test_intervals()
    test_combs_and_contexts()
    test_mutations()
