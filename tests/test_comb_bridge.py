#!/usr/bin/env python3
"""Original-game checks for two arbitrary cores and a comb bridge."""

from copy import deepcopy
from itertools import combinations
import json
from pathlib import Path
import random
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from comb_bridge import CombBridge
from strip3_defect_algebra import full_head
from verify_comb_bridge import original_game, verify


def graphs(n):
    possible = list(combinations(range(n), 2))
    for code in range(1 << len(possible)):
        graph = [0] * n
        for i, (u, v) in enumerate(possible):
            if code >> i & 1:
                graph[u] |= 1 << v
                graph[v] |= 1 << u
        yield tuple(graph)


def check(model, left, right, length, phase):
    expected, all_options = original_game(left, right, length, phase)
    assert model.value(length, phase) == expected
    if length:
        representatives = model.options(length, phase)
        assert {value for _, value in representatives} == set(all_options.values())
        for move, value in representatives:
            assert all_options[move] == value


def main():
    cores = [graph for n in range(1, 4) for graph in graphs(n)]
    exhaustive = random_cases = periodic_cases = 0
    for left in cores:
        for right in cores:
            model = CombBridge(left, right)
            for phase in (0, 1):
                for length in range(10):
                    check(model, left, right, length, phase)
                    exhaustive += 1
                for length in range(model.period_start, model.period_start + 8):
                    check(model, left, right, length, phase)
                    assert model.value(length, phase) == model.value(length + 4, phase)
                    periodic_cases += 1
    rng = random.Random(221026)
    for _ in range(50):
        pair = []
        for side in range(2):
            n = rng.randrange(1, 7)
            graph = [0] * n
            for u, v in combinations(range(n), 2):
                if rng.randrange(3) == 0:
                    graph[u] |= 1 << v
                    graph[v] |= 1 << u
            pair.append(tuple(graph))
        model = CombBridge(*pair)
        for phase in (0, 1):
            for length in (0, 1, 2, 7, model.period_start + rng.randrange(8)):
                check(model, *pair, length, phase)
                random_cases += 1
        huge = 10 ** 18 + rng.randrange(16)
        for phase in (0, 1):
            reduced = model.period_start + (huge - model.period_start) % 4
            assert model.value(huge, phase) == model.value(reduced, phase)
            assert len(model.options(huge, phase)) <= sum(map(len, pair)) + 2 * (
                model.left_start + model.right_start + 7)

    # Every odd bridge length is an actual history in an empty three-row board.
    reachable_cases = []
    for left_width, right_width in ((1, 1), (1, 4), (2, 4), (3, 4)):
        left, right = full_head(left_width)[1], full_head(right_width)[1]
        model = CombBridge(left, right)
        for n in range(1, 20, 2):
            width = left_width + n + right_width
            available = {(r, c) for c in range(width) for r in range(3)}
            history = [(0, left_width + i) for i in range(0, n, 2)]
            for r, c in history:
                assert (r, c) in available
                available -= {(r, c), (r - 1, c), (r + 1, c), (r, c - 1), (r, c + 1)}
            left_cells = [(r - 1, c - 1) for r, c in full_head(left_width)[0]]
            right_cells = [(r - 1, width - c) for r, c in full_head(right_width)[0]]
            actual_cells = set(left_cells + right_cells)
            actual_cells |= {(2, left_width + i) for i in range(n)}
            actual_cells |= {(1, left_width + i) for i in range(1, n, 2)}
            assert actual_cells == available
            # In the actual grid each core meets the bridge by just one edge.
            bridge_cells = actual_cells - set(left_cells + right_cells)
            cross = {(v, u) for v in left_cells + right_cells for u in bridge_cells
                     if abs(v[0] - u[0]) + abs(v[1] - u[1]) == 1}
            assert cross == {(left_cells[0], (2, left_width)),
                             (right_cells[0], (2, left_width + n - 1))}
        reachable_cases.append({"left_head": left_width, "right_head": right_width,
                                "period_start": model.period_start,
                                "odd_residue_values": {str(r): model.value(
                                    model.period_start + (r - model.period_start) % 4, 0)
                                    for r in (1, 3)}})

    example = CombBridge(full_head(2)[1], full_head(4)[1])
    proof = example.certificate()
    result = verify(proof)
    tampered = []
    bad = deepcopy(proof); bad["period_start"] -= 1; tampered.append(bad)
    bad = deepcopy(proof); bad["one_ended_tables"][0]["states"].pop(); tampered.append(bad)
    bad = deepcopy(proof); bad["one_ended_tables"][0]["states"][0]["values"][0] ^= 1; tampered.append(bad)
    bad = deepcopy(proof); bad["period_values"]["0"][0] ^= 1; tampered.append(bad)
    for bad in tampered:
        try:
            verify(bad)
        except ValueError:
            continue
        raise AssertionError("tampered bridge certificate accepted")
    (ROOT / "data/proofs/comb_bridge_2_4.json").write_text(json.dumps(proof, indent=2) + "\n")
    audit = {"format": "comb-bridge-audit-v1", "exhaustive_original_games": exhaustive,
             "periodic_original_games": periodic_cases, "random_original_games": random_cases,
             "reachable_families": reachable_cases, "certificate_checks": result,
             "tamper_rejections": len(tampered)}
    (ROOT / "data/proofs/comb_bridge_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps(audit))


if __name__ == "__main__":
    main()
