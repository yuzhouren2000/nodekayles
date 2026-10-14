#!/usr/bin/env python3
"""Real histories, raw SG comparisons and infinite-bridge reader checks."""

from functools import cache
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from parametric_comb_bridge import (read_proof, reduce_bridge, thick_chain_columns,
                                    thick_chain_history, uniform_onset)
from strip3_ports import run_masks
from two_comb_core import TwoCombCore
from verify_parametric_comb_bridge import verify_uniform
from verify_two_comb_core import verify as verify_base


def raw_sg(columns):
    cells = frozenset((r, c) for c, col in enumerate(columns) for r in range(3) if col >> r & 1)
    edges = {v: {u for u in cells if abs(v[0] - u[0]) + abs(v[1] - u[1]) == 1}
             for v in cells}

    @cache
    def sg(position):
        if not position:
            return 0
        seen, pending = {next(iter(position))}, [next(iter(position))]
        while pending:
            v = pending.pop()
            added = (edges[v] & position) - seen
            seen |= added
            pending.extend(added)
        if len(seen) != len(position):
            return sg(frozenset(seen)) ^ sg(position - seen)
        outcomes = {sg(position - edges[v] - {v}) for v in position}
        value = 0
        while value in outcomes:
            value += 1
        return value
    return sg(cells)


def main():
    proofs = {w: read_proof(ROOT / f"data/proofs/two_comb_strip_middle{w}.json.gz")
              for w in (2, 3, 4, 5)}
    onsets = {w: uniform_onset(p) for w, p in proofs.items()}
    assert onsets == {2: 8, 3: 9, 4: 8, 5: 9}
    for width, proof in proofs.items():
        assert verify_uniform(proof, check_base=False)["uniform_onset"] == onsets[width]
    histories = 0
    for w, u in ((2, 2), (2, 3), (3, 5), (5, 4)):
        for n in range(1, 8):
            for m in range(1, 8):
                for ell in (1, 3, 5, 23, 25):
                    cols = thick_chain_columns(w, u, n, m, ell)
                    size, history = thick_chain_history(w, u, n, m, ell)
                    available = {(r, c) for c in range(size) for r in range(3)}
                    for r, c in history:
                        assert (r, c) in available
                        available -= {(r, c), (r - 1, c), (r + 1, c), (r, c - 1), (r, c + 1)}
                    actual = tuple(sum(1 << r for r in range(3) if (r, c) in available)
                                   for c in range(size))
                    assert cols == actual
                    histories += 1

    positions, pairings, raw_inputs = [], [], []
    def add(cols):
        index = len(positions)
        positions.append((len(cols), sum(x << (3 * i) for i, x in enumerate(cols))))
        return index
    for w, u in ((2, 2), (2, 3), (2, 4), (2, 5), (3, 3), (3, 4), (4, 4)):
        start = onsets[w] + onsets[u] + 6
        odd_start = start if start % 2 else start + 1
        for n in range(1, 6):
            for m in range(1, 6):
                for ell in range(odd_start, odd_start + 5, 2):
                    a = thick_chain_columns(w, u, n, m, ell)
                    b = thick_chain_columns(w, u, n, m, ell + 4)
                    if len(b) <= 42:
                        pairings.append((add(a), add(b)))
    for n in (1, 2):
        for m in (1, 2):
            for ell in (1, 3):
                cols = thick_chain_columns(2, 2, n, m, ell)
                raw_inputs.append((add(cols), raw_sg(cols)))
    # A real reachable core shows why a small onset cannot be assumed for
    # every thick geometry. These are finite counterexamples, not a new
    # arbitrary-width period assertion.
    counter_core = tuple(map(int, "277646313776"))
    counter_indices = []
    for length in (12, 16):
        cols = counter_core + tuple(4 if i % 2 == 0 else 6 for i in range(length))
        counter_indices.append(add(cols))
        from strip3_history_language import decode
        code = decode(sum(x << (3 * i) for i, x in enumerate(cols)))
        size, history, target = code.witness()
        available = {(r, c) for c in range(size) for r in range(3)}
        for r, c in history:
            assert (r, c) in available
            available -= {(r, c), (r - 1, c), (r + 1, c), (r, c - 1), (r, c + 1)}
        desired = {(i % 3, i // 3) for i in range(target.bit_length()) if target >> i & 1}
        assert desired <= available
        assert not any(abs(r - s) + abs(c - d) == 1
                       for r, c in desired for s, d in available - desired)
    values, provenance = run_masks(positions, ht_bits=25, timeout=60)
    assert provenance["returncode"] == 0 and None not in values, provenance
    for i, j in pairings:
        assert values[i] == values[j], (positions[i], positions[j], values[i], values[j])
    for index, expected in raw_inputs:
        assert values[index] == expected
    assert [values[i] for i in counter_indices] == [9, 13]

    huge = 10 ** 100 + 1
    result = reduce_bridge(proofs[2], proofs[3], huge)
    assert result["bridge_onset"] == 23 and result["reduced_length"] == 25
    # An independently certified phase-zero counterexample keeps a diagonal
    # band forever, so this extra module premise must reject it.
    absorbed = TwoCombCore((26, 5, 74, 5, 33, 16, 132, 64), (4, 6), (0, 0)).build(32)
    diagonal_proof = absorbed.candidate(8, 8, 16).proof()
    verify_base(diagonal_proof)  # Reject the extra premise, not a faulty base certificate.
    rejections = 0
    for check in (uniform_onset, verify_uniform):
        try:
            check(diagonal_proof)
        except ValueError:
            rejections += 1
        else:
            raise AssertionError("a diagonal-only module was accepted as uniformly periodic")
    assert rejections == 2
    for length in (0, -1, 2, 10):
        try:
            reduce_bridge(proofs[2], proofs[3], length)
        except ValueError:
            continue
        raise AssertionError("invalid reachable bridge length accepted")
    audit = {"format": "node-kayles-parametric-bridge-test-audit-v1",
             "uniform_onsets": onsets, "real_histories": histories,
             "original_C_period_pairs": len(pairings), "raw_subset_SG_checks": len(raw_inputs),
             "C_provenance": provenance, "huge_bridge_reduction": result,
             "diagonal_module_rejections": rejections,
             "real_onset_counterexample": {"core": "277646313776", "lengths": [12, 16],
                                             "SG_values": [9, 13]}}
    (ROOT / "data/proofs/parametric_bridge_test_audit.json").write_text(
        json.dumps(audit, indent=2) + "\n")
    print(json.dumps(audit), flush=True)


if __name__ == "__main__":
    main()
