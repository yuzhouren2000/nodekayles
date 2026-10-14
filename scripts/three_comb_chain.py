#!/usr/bin/env python3
"""Certified SG of two thick blocks with three arbitrary comb lengths.

The bundle contains both uniform module premises and ALL reduced odd bridges.
No arm or long bridge is materialized during an SG query.
"""

import argparse
import gc
import gzip
import json
from pathlib import Path
import time

from connected_comb_core import CertifiedConnectedTwoComb, ConnectedTwoCombCore
from kayles_boundary_transfer import column_graph
from parametric_comb_bridge import read_proof, reduce_bridge

ROOT = Path(__file__).resolve().parents[1]
DEFAULT = ROOT / "data/proofs/three_comb_chain_2_3.json.gz"


def chain_core(left_width, right_width, bridge):
    columns = ((6,) + (7,) * left_width + (6,)
               + tuple(4 if i % 2 == 0 else 6 for i in range(bridge))
               + (6,) + (7,) * right_width + (6,))
    graph, cells = column_graph(columns)
    roots = cells.index((2, 0)), cells.index((2, len(columns) - 1))
    return graph, roots


def build_bundle(left_width=2, right_width=3, sample=64):
    modules = [read_proof(ROOT / f"data/proofs/two_comb_strip_middle{w}.json.gz")
               for w in (left_width, right_width)]
    start = reduce_bridge(*modules, 1)["bridge_onset"]
    # Each odd length below start+4 is needed, including ALL short bridges.
    proofs = []
    for bridge in range(1, start + 4, 2):
        began = time.monotonic()
        graph, roots = chain_core(left_width, right_width, bridge)
        model = ConnectedTwoCombCore(graph, roots).build(sample)
        for minimum in (12, 16, 20, 24, 32):
            candidate = model.candidate(minimum, 4, max(24, minimum))
            if candidate.sample_agrees():
                break
        else:
            raise ValueError(f"no sampled partition for bridge {bridge}")
        proofs.append({"bridge": bridge, "proof": candidate.proof()})
        print(json.dumps({"bridge": bridge, "connected_states": len(model.states),
                          "unary_states": [len(bank.states) for bank in model.banks],
                          "unary_starts": [bank.period_start() for bank in model.banks],
                          "thresholds": [candidate.minimum, candidate.difference, candidate.boundary],
                          "seconds": round(time.monotonic() - began, 2)}), flush=True)
        del model, candidate
        gc.collect()
    return {"format": "node-kayles-thick-three-comb-v1", "widths": [left_width, right_width],
            "module_proofs": modules, "bridge_onset": start, "connected_proofs": proofs}


class CertifiedThreeComb:
    """Finite query reader, requiring prior independent bundle verification."""

    def __init__(self, bundle):
        if bundle.get("format") != "node-kayles-thick-three-comb-v1":
            raise ValueError("unknown three-comb bundle")
        self.bundle = bundle
        self.start = reduce_bridge(*bundle["module_proofs"], 1)["bridge_onset"]
        if self.start != bundle["bridge_onset"]:
            raise ValueError("bundle disagrees with its module bridge onset")
        self.queries = {item["bridge"]: CertifiedConnectedTwoComb(item["proof"])
                        for item in bundle["connected_proofs"]}

    def query(self, left, right, bridge):
        if any(type(n) is not int or n < 0 for n in (left, right)):
            raise ValueError("outer lengths must be nonnegative integers")
        if type(bridge) is not int or bridge < 1 or not bridge % 2:
            raise ValueError("bridge must have positive odd length")
        representative = bridge if bridge < self.start else self.start + (bridge-self.start) % 4
        query = self.queries[representative]
        return {"lengths": {"left": left, "right": right, "bridge": bridge},
                "reduced_bridge": representative,
                "reduced_outer_lengths": query.candidate.canonical((left, right)),
                "SG": query.value((left, right)),
                "bridge_onset": self.start}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", action="store_true", help="generate finite premises; then run verifier")
    parser.add_argument("--certificate", type=Path, default=DEFAULT)
    parser.add_argument("--widths", nargs=2, type=int, choices=(2, 3, 4, 5), default=(2, 3))
    parser.add_argument("--sample", type=int, default=64)
    parser.add_argument("--lengths", nargs=3, type=int, metavar=("LEFT", "RIGHT", "BRIDGE"))
    args = parser.parse_args()
    if args.build:
        bundle = build_bundle(*args.widths, args.sample)
        output = (json.dumps(bundle, separators=(",", ":")) + "\n").encode()
        args.certificate.write_bytes(gzip.compress(output, mtime=0) if args.certificate.suffix == ".gz" else output)
    elif args.lengths:
        try:
            result = CertifiedThreeComb(read_proof(args.certificate)).query(*args.lengths)
        except ValueError as error:
            parser.error(str(error))
        print(json.dumps(result))
    else:
        parser.error("--build or --lengths required")


if __name__ == "__main__":
    main()
