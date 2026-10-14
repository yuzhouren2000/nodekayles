#!/usr/bin/env python3
"""Uniform bridge periods between certified two-arm parameterized modules.

The certificate reader checks the extra uniform premises, not the base SG
certificate. Run tests/verify_parametric_comb_bridge.py for that independent
verification. No long graph or remote computation is needed for reduction.
"""

import argparse
import gzip
import json
from pathlib import Path

from two_comb_core import CertifiedTwoComb


def read_proof(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as source:
        return json.load(source)


def uniform_onset(proof):
    """Prove a sufficient onset from the already certified finite partition.

    Both phases must be zero. Every stored two-root descendant must have an
    independent period once both coordinates are >= minimum. Short outer
    arms and all recorded one-root descendants are checked separately.
    """
    query = CertifiedTwoComb(proof)
    if query.model.phases != (0, 0):
        raise ValueError("this theorem requires two phase-zero module arms")
    b, t = proof["minimum"], proof["boundary"]
    limit = proof["sample_bound"]
    for mask in query.model.rows:
        for n in range(b, limit + 1):
            for m in range(b, limit + 1):
                if query.value((n, m), mask) != query.value((b + (n - b) % 4,
                                                            b + (m - b) % 4), mask):
                    raise ValueError("two long arms still depend on their length difference")
    unary_bound = max((table["start"] for table in proof["unary"]), default=2)
    onset = max(b, unary_bound)
    # The independent long-long chamber is covered above. In a short-long
    # chamber the candidate is periodic after t; the zero axis is unary.
    for mask in query.model.rows:
        for short in range(b):
            for length in range(b, max(t, unary_bound) + 4):
                for side in (0, 1):
                    a = (short, length) if side == 0 else (length, short)
                    c = (short, length + 4) if side == 0 else (length + 4, short)
                    if query.value(a, mask) != query.value(c, mask):
                        onset = max(onset, length + 1)
    return onset


def reduce_bridge(left_proof, right_proof, length):
    if type(length) is not int or length < 1 or length % 2 == 0:
        raise ValueError("the reachable phase-zero bridge length must be positive and odd")
    left_onset, right_onset = uniform_onset(left_proof), uniform_onset(right_proof)
    start = left_onset + right_onset + 6
    reduced = length if length < start else start + (length - start) % 4
    return {"left_uniform_onset": left_onset, "right_uniform_onset": right_onset,
            "bridge_period": 4, "bridge_onset": start, "length": length,
            "reduced_length": reduced,
            "scope": "same SG for every outer-arm pair, including all parameterized root-safe descendants"}


def thick_chain_columns(left_width, right_width, left_arm, right_arm, bridge):
    if (min(left_width, right_width, left_arm, right_arm, bridge) < 1
            or bridge % 2 == 0):
        raise ValueError("positive widths/arms and an odd bridge required")
    arm = lambda n: tuple(4 if i % 2 == 0 else 6 for i in range(n))
    return (arm(left_arm)[::-1] + (6,) + (7,) * left_width + (6,)
            + arm(bridge) + (6,) + (7,) * right_width + (6,) + arm(right_arm))


def thick_chain_history(left_width, right_width, left_arm, right_arm, bridge):
    """Actual zero-based top-row history on an empty three-row board."""
    n, m, w, u, ell = left_arm, right_arm, left_width, right_width, bridge
    thick_chain_columns(w, u, n, m, ell)  # Validate the geometry.
    total = n + m + ell + w + u + 4
    columns = list(range(n - 1, -1, -2))
    columns += list(range(n + w + 2, n + w + 2 + ell, 2))
    columns += list(range(n + w + ell + u + 4, total, 2))
    return total, tuple((0, c) for c in columns)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--left-certificate", type=Path, required=True)
    parser.add_argument("--right-certificate", type=Path, required=True)
    parser.add_argument("--bridge", type=int, required=True)
    args = parser.parse_args()
    try:
        result = reduce_bridge(read_proof(args.left_certificate),
                               read_proof(args.right_certificate), args.bridge)
    except ValueError as error:
        parser.error(str(error))
    print(json.dumps(result))


if __name__ == "__main__":
    main()
