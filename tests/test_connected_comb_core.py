#!/usr/bin/env python3
"""Raw named-vertex, old-representation, and certificate mutation checks."""

from copy import deepcopy
from itertools import combinations
import json
from pathlib import Path
import random
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from connected_comb_core import CertifiedConnectedTwoComb, ConnectedTwoCombCore
from two_comb_core import TwoCombCore
from test_two_comb_core import original, literal_candidate_options
from verify_connected_comb_core import ConnectedChecker


def main():
    rng = random.Random(251014)
    raw_checks = old_checks = certified_checks = spectra = 0
    graphs = []
    for k in (2, 3):
        possible = list(combinations(range(k), 2))
        for code in range(1 << len(possible)):
            graph = [0] * k
            for i, (v, u) in enumerate(possible):
                if code >> i & 1:
                    graph[v] |= 1 << u
                    graph[u] |= 1 << v
            graphs.append((graph, (0, k - 1)))
    for _ in range(24):
        k = rng.randrange(2, 7)
        graph = [0] * k
        for v, u in combinations(range(k), 2):
            if rng.randrange(3) == 0:
                graph[v] |= 1 << u
                graph[u] |= 1 << v
        graphs.append((graph, tuple(rng.sample(range(k), 2))))
    square = [10, 5, 10, 5]
    graphs.append((square, (0, 2)))
    for graph, roots in graphs:
        for phases in ((0, 0), (0, 1), (1, 0), (1, 1)):
            model = ConnectedTwoCombCore(graph, roots, phases).build(5)
            old = TwoCombCore(graph, roots, phases).build(5)
            for n in range(6):
                for m in range(6):
                    lengths = n, m
                    expected, options = original(graph, roots, phases, lengths)
                    actual = model.literal_options(model.full, lengths)
                    assert model.value(model.full, lengths) == expected
                    assert {g for _, g in actual} == set(options.values())
                    assert all(options[move] == g for move, g in actual)
                    raw_checks += 1
                    # Include all raw descendants, with their detached fragments.
                    for mask in old.states:
                        assert old.value(mask, lengths) == model.value(mask, lengths)
                        old_checks += 1

    # Includes connected and disconnected initial roots, and every phase pair.
    certs = []
    chain = [sum(1 << u for u in range(6) if abs(v - u) == 1) for v in range(6)]
    for graph, roots in graphs[:10] + [(chain, (0, 1)), (square, (0, 2))]:
        for phases in ((0, 0), (0, 1), (1, 0), (1, 1)):
            if graph == square and phases in ((0,1),(1,0)):
                continue  # These mixed-phase samples do not fit this partition.
            model = ConnectedTwoCombCore(graph, roots, phases).build(64)
            candidate = model.candidate(12, 12, 24)
            assert candidate.sample_agrees()
            proof = candidate.proof()
            checker = ConnectedChecker(proof)
            checker.verify()
            query = CertifiedConnectedTwoComb(proof)
            for n in range(6):
                for m in range(6):
                    expected, options = original(graph, roots, phases, (n, m))
                    actual = query.representative_options((n, m))
                    assert query.value((n, m)) == expected
                    assert {g for _, g in actual} == set(options.values())
                    assert all(options[move] == g for move, g in actual)
                    move = query.winning_move((n, m))
                    assert (move is None) == (expected == 0)
                    if move is not None:
                        assert options[move] == 0
                    certified_checks += 1
            for _ in range(8):
                n, m = rng.randrange(10**18), rng.randrange(10**18)
                assert query.value((n, m)) == checker.at(checker.full, n, m)
                assert {g for _, g in query.representative_options((n, m))} == checker.options(checker.full, n, m)
                a, b = checker.refined(n, m)
                assert checker.options(checker.full, n, m) == checker.options(checker.full, a, b)
                n, m = rng.randrange(180), rng.randrange(180)
                assert checker.options(checker.full, n, m) == literal_candidate_options(checker, checker.full, n, m)
                spectra += 1
            certs.append(proof)
    mixed_samples = []
    for phases in ((0,1),(1,0)):
        model = ConnectedTwoCombCore(square, (0,2), phases).build(64)
        old = TwoCombCore(square, (0,2), phases).build(64)
        assert not model.candidate(12,12,24).sample_agrees()
        values = [model.value(model.full,(n,n)) for n in (16,32,64)]
        assert values == [old.value(old.full,(n,n)) for n in (16,32,64)]
        mixed_samples.append({"phases":phases,"diagonal_lengths":[16,32,64],"SG":values,
                              "status":"finite exact samples; no unboundedness theorem asserted"})
    # Normalization must retain the known PERMANENT diagonal obstruction,
    # including its phase-zero realization; it is not an independent-period rule.
    absorbed = ConnectedTwoCombCore((26,5,74,5,33,16,132,64), (4,6), (0,0)).build(64)
    candidate = absorbed.candidate(12,12,24)
    assert candidate.sample_agrees()
    proof = candidate.proof()
    ConnectedChecker(proof).verify()
    query = CertifiedConnectedTwoComb(proof)
    n = 10**100
    assert query.value((n,n)) == 0
    assert query.value((n+4,n)) == 2
    assert query.value((n+2,n)) == 5
    certs.append(proof)
    proof = next(p for p in reversed(certs) if len(p["states"]) > 1)
    mutations = []
    bad = deepcopy(proof); bad["states"].pop(); mutations.append(bad)
    bad = deepcopy(proof); bad["states"][0]["rows"][12][12] ^= 1; mutations.append(bad)
    bad = deepcopy(proof); bad["unary_banks"][0]["states"].pop(); mutations.append(bad)
    bad = deepcopy(proof); bad["unary_banks"][0]["seeds"].pop(); mutations.append(bad)
    bad = deepcopy(proof); bad["unary_banks"][0]["states"][0]["values"][0] ^= 1; mutations.append(bad)
    for bad in mutations:
        try:
            ConnectedChecker(bad).verify()
        except ValueError:
            continue
        raise AssertionError("altered connected certificate accepted")
    audit = {"raw_named_graph_checks": raw_checks, "raw_descendant_cross_checks": old_checks,
             "verified_certificates": len(certs), "query_original_checks": certified_checks,
             "literal_and_refined_spectra": spectra, "tamper_rejections": len(mutations),
             "mixed_phase_square_samples": mixed_samples,
             "permanent_phase_zero_difference": {"n":n,"equal":0,"difference4":2,"difference2":5}}
    (ROOT / "data/proofs/connected_comb_core_test_audit.json").write_text(json.dumps(audit,indent=2)+"\n")
    print(json.dumps(audit))


if __name__ == "__main__":
    main()
