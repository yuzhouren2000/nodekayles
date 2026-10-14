#!/usr/bin/env python3
"""Whole-strip SG, legal histories, large-integer and exceptional cases."""

import json
from pathlib import Path
import random
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from parametric_comb_bridge import read_proof, thick_chain_columns, thick_chain_history
from strip3_ports import run_masks
from three_comb_chain import CertifiedThreeComb
from test_parametric_comb_bridge import raw_sg
from verify_connected_comb_core import ConnectedChecker


def main():
    bundle = read_proof(ROOT / "data/proofs/three_comb_chain_2_3.json.gz")
    query = CertifiedThreeComb(bundle)
    histories = raw_checks = 0
    positions, expected, parameters = [], [], []
    def add(n, m, ell):
        cols = thick_chain_columns(2, 3, n, m, ell)
        assert len(cols) <= 42
        positions.append((len(cols), sum(bits << (3*i) for i,bits in enumerate(cols))))
        expected.append(query.query(n, m, ell)["SG"])
        parameters.append((n, m, ell))
    # Every reduced bridge, axes as abstract graphs, and actual nonzero arms.
    for ell in range(1, 26, 2):
        for n in (1, 2, 5, 7, 8, 11, 12):
            for m in (1, 3, 5, 7, 8, 11, 12):
                cols = thick_chain_columns(2, 3, n, m, ell)
                size, history = thick_chain_history(2, 3, n, m, ell)
                available = {(r,c) for r in range(3) for c in range(size)}
                for r,c in history:
                    assert (r,c) in available
                    available -= {(r,c), (r-1,c), (r+1,c), (r,c-1), (r,c+1)}
                actual = tuple(sum(1 << r for r in range(3) if (r,c) in available) for c in range(size))
                assert actual == cols
                histories += 1
                if size <= 42:
                    add(n, m, ell)
    for n in (1, 2):
        for m in (1, 2):
            for ell in (1, 3):
                cols = thick_chain_columns(2, 3, n, m, ell)
                assert query.query(n,m,ell)["SG"] == raw_sg(cols)
                raw_checks += 1
    # Direct comparison with the original C solver, including bridge >25.
    for ell in (27, 29, 31):
        add(1,1,ell)
    values, provenance = run_masks(positions, ht_bits=25, timeout=60)
    assert provenance["returncode"] == 0 and None not in values, provenance
    for actual, wanted, p in zip(values, expected, parameters):
        assert actual == wanted, (p, actual, wanted)

    rng = random.Random(253214)
    r = (0, 3, 2, 1)
    huge_checks = 0
    for _ in range(1000):
        n, m = rng.randrange(5,10**100), rng.randrange(5,10**100)
        ell = 2*rng.randrange(5,10**100)+1
        result = query.query(n,m,ell)
        assert result["SG"] == r[n%4] ^ r[m%4] ^ ((ell+1)%4)
        assert (result["SG"] == 0) == ((n-m-ell-1)%4 == 0)
        huge_checks += 1
    # The smaller P onset does not imply the full SG formula at ell=9.
    assert query.query(7,7,9)["SG"] == 4
    assert (r[7%4] ^ r[7%4] ^ 2) == 2
    assert query.query(7,7,3)["SG"] == 4
    assert query.query(7,7,5)["SG"] == 0
    assert query.query(4,4,11)["SG"] == 10
    for n in range(5,21):
        for m in range(5,21):
            for ell in (7,9,11,13,23,25,10**100+1):
                sg = query.query(n,m,ell)["SG"]
                assert (sg == 0) == ((n-m-ell-1)%4 == 0)
                if ell >= 11:
                    assert sg == r[n%4] ^ r[m%4] ^ ((ell+1)%4)
    # Every stored representative reads the same abstract zero axes as its checker.
    axes = 0
    for item in bundle["connected_proofs"]:
        checker = ConnectedChecker(item["proof"])
        for n,m in ((0,0),(0,1),(1,0),(0,10**100),(10**100,0)):
            assert query.query(n,m,item["bridge"])["SG"] == checker.at(checker.full,n,m)
            axes += 1
    invalids = [(1,1,0),(1,1,2),(-1,1,3),(1,-1,3),(1.5,1,3)]
    for p in invalids:
        try:
            query.query(*p)
        except ValueError:
            continue
        raise AssertionError("invalid lengths accepted")
    audit = {"format": "node-kayles-thick-three-comb-test-audit-v1",
             "real_histories": histories, "original_C_SG_checks": len(positions),
             "named_vertex_SG_checks": raw_checks, "huge_formula_checks": huge_checks,
             "abstract_zero_axis_checks": axes, "invalid_inputs_rejected": len(invalids),
             "short_bridge_SG_counterexample": {"left":7,"right":7,"bridge":9,"SG":4},
             "C_provenance": provenance}
    (ROOT / "data/proofs/three_comb_chain_2_3_test_audit.json").write_text(json.dumps(audit,indent=2)+"\n")
    print(json.dumps(audit))


if __name__ == "__main__":
    main()
