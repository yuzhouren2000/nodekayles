#!/usr/bin/env python3
"""Verify an entire three-length result, including all short bridge cases."""

import gzip
import hashlib
import json
from pathlib import Path
import sys

from verify_connected_comb_core import ConnectedChecker
from verify_parametric_comb_bridge import verify_uniform

ROOT = Path(__file__).resolve().parents[1]


def coordinate_graph(columns):
    cells = [(r, col) for col, bits in enumerate(columns) for r in range(3) if bits >> r & 1]
    graph = tuple(sum(1 << j for j, u in enumerate(cells)
                      if abs(v[0] - u[0]) + abs(v[1] - u[1]) == 1) for v in cells)
    roots = cells.index((2, 0)), cells.index((2, len(columns) - 1))
    return graph, roots


def check_geometry(proof, columns):
    graph, roots = coordinate_graph(columns)
    if (tuple(proof["adjacency"]) != graph or tuple(proof["roots"]) != roots
            or tuple(proof["phases"]) != (0, 0)):
        raise ValueError("certificate has the wrong thick-chain graph or roots")


def verify(bundle):
    if bundle.get("format") != "node-kayles-thick-three-comb-v1":
        raise ValueError("unknown three-comb bundle format")
    widths = bundle["widths"]
    if len(widths) != 2 or any(type(w) is not int or w not in (2, 3, 4, 5) for w in widths):
        raise ValueError("unsupported thicknesses")
    if len(bundle["module_proofs"]) != 2:
        raise ValueError("two module premises required")
    modules = []
    for width, proof in zip(widths, bundle["module_proofs"]):
        check_geometry(proof, (6,) + (7,) * width + (6,))
        modules.append(verify_uniform(proof))
    start = sum(module["uniform_onset"] for module in modules) + 6
    if bundle["bridge_onset"] != start:
        raise ValueError("wrong uniform bridge onset")
    required = set(range(1, start + 4, 2))
    bridges = [item["bridge"] for item in bundle["connected_proofs"]]
    if (any(type(ell) is not int for ell in bridges) or len(bridges) != len(required)
            or set(bridges) != required):
        raise ValueError("missing, extra or duplicate representative bridge")
    records = []
    for item in bundle["connected_proofs"]:
        ell, proof = item["bridge"], item["proof"]
        columns = ((6,) + (7,) * widths[0] + (6,)
                   + tuple(4 if i % 2 == 0 else 6 for i in range(ell))
                   + (6,) + (7,) * widths[1] + (6,))
        check_geometry(proof, columns)
        checker = ConnectedChecker(proof)
        audit = checker.verify()
        bound = proof["sample_bound"]
        # Check ALL canonical cells, not a sample of remote coordinates.
        minimum = None
        for b in range(1, checker.b + 1):
            matrix = [[checker.at(checker.full, b + i, b + j) for j in range(4)] for i in range(4)]
            if all(checker.at(checker.full, n, m) == matrix[(n-b) % 4][(m-b) % 4]
                   for n in range(b, bound + 1) for m in range(b, bound + 1)):
                minimum = b
                break
        if minimum is not None:
            matrix_mod4 = [[checker.at(checker.full, minimum + (i-minimum) % 4,
                                       minimum + (j-minimum) % 4) for j in range(4)] for i in range(4)]
            audit["full_core_independent_onset"] = minimum
            audit["full_core_matrix_mod4"] = matrix_mod4
            audit["P_pairs_mod4"] = [[i,j] for i in range(4) for j in range(4) if not matrix_mod4[i][j]]
        # Descendants determine whether this entire cluster can act as a module.
        audit["all_connected_descendants_independent"] = all(
            checker.at(mask, n, m) == checker.at(mask, checker.b + (n-checker.b) % 4,
                                                checker.b + (m-checker.b) % 4)
            for mask in checker.rows for n in range(checker.b, bound + 1)
            for m in range(checker.b, bound + 1))
        audit["bridge"] = ell
        records.append(audit)
        print(json.dumps(audit), flush=True)
    result = {"format": "node-kayles-thick-three-comb-audit-v1", "widths": widths,
            "bridge_onset": start, "module_premises": modules, "bridges": records,
            "total_finite_mex_identities": sum(r["finite_mex_identities"] for r in records),
            "total_finite_unary_mex_identities": sum(r["finite_unary_mex_identities"] for r in records)}
    if widths == [2, 3]:
        r = (0, 3, 2, 1)
        for record in records:
            ell = record["bridge"]
            if ell >= 7:
                if record.get("full_core_independent_onset", 100) > 5:
                    raise ValueError("claimed five-column outer onset fails")
                expected_P = [[i,j] for i in range(4) for j in range(4)
                              if (i-j-ell-1) % 4 == 0]
                if record["P_pairs_mod4"] != expected_P:
                    raise ValueError("claimed three-length P formula fails")
            if ell >= 11:
                expected = [[r[i] ^ r[j] ^ ((ell+1) % 4) for j in range(4)] for i in range(4)]
                if record["full_core_matrix_mod4"] != expected:
                    raise ValueError("claimed three-length SG formula fails")
        result["closed_SG_formula"] = {"outer_minimum": 5, "odd_bridge_minimum": 11,
                                       "r_mod4": r, "SG": "r(n) xor r(m) xor ((ell+1) mod4)"}
        result["closed_P_formula"] = {"outer_minimum": 5, "odd_bridge_minimum": 7,
                                      "condition": "n-m == ell+1 (mod4)"}
    return result


if __name__ == "__main__":
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "data/proofs/three_comb_chain_2_3.json.gz"
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as source:
        result = verify(json.load(source))
    result["certificate_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    output = ROOT / "data/proofs/three_comb_chain_2_3_audit.json"
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key:value for key,value in result.items() if key not in ("bridges", "module_premises")}))
