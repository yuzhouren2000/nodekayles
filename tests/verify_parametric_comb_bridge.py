#!/usr/bin/env python3
"""Independent extra premises for uniform parameterized module bridges.

The base proofs are independently checked by verify_two_comb_core. This file
does not import the production bridge reducer or candidate reader.
"""

import gzip
import hashlib
import json
from pathlib import Path
import sys

from verify_two_comb_core import Checker

ROOT = Path(__file__).resolve().parents[1]


def verify_uniform(proof, check_base=True):
    checker = Checker(proof, check_unary=check_base)
    base = checker.verify() if check_base else None
    if checker.phases != (0, 0):
        raise ValueError("uniform phase-zero module premises required")
    long_checks = short_checks = 0
    for mask in checker.rows:
        matrix = [[checker.at(mask, checker.b + i, checker.b + j) for j in range(4)]
                  for i in range(4)]
        for n in range(checker.b, proof["sample_bound"] + 1):
            for m in range(checker.b, proof["sample_bound"] + 1):
                if checker.at(mask, n, m) != matrix[(n - checker.b) % 4][(m - checker.b) % 4]:
                    raise ValueError("a long-long canonical chamber depends on exact difference")
                long_checks += 1
    onset = max(checker.b, checker.unary_bound)
    for mask in checker.rows:
        for short in range(checker.b):
            for length in range(checker.b, max(checker.t, checker.unary_bound) + 4):
                for side in (0, 1):
                    a = (short, length) if side == 0 else (length, short)
                    b = (short, length + 4) if side == 0 else (length + 4, short)
                    if checker.at(mask, *a) != checker.at(mask, *b):
                        onset = max(onset, length + 1)
                    short_checks += 1
    return {"uniform_onset": onset, "core_states": len(checker.rows),
            "long_chamber_checks": long_checks, "short_axis_checks": short_checks,
            "base_certificate": base}


def main():
    records = []
    for width in (2, 3, 4, 5):
        path = ROOT / f"data/proofs/two_comb_strip_middle{width}.json.gz"
        with gzip.open(path, "rt") as source:
            result = verify_uniform(json.load(source))
        result.update({"complete_middle_columns": width,
                       "certificate_sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
        records.append(result)
        print(json.dumps(result), flush=True)
    output = {"format": "node-kayles-parametric-bridge-uniform-audit-v1",
              "modules": records,
              "bridge_onsets": [[a["uniform_onset"] + b["uniform_onset"] + 6
                                  for b in records] for a in records]}
    (ROOT / "data/proofs/parametric_bridge_uniform_audit.json").write_text(
        json.dumps(output, indent=2) + "\n")


if __name__ == "__main__":
    main()
