#!/usr/bin/env python3
"""Verify finite premises of the rooted-comb period-lifting theorem.

No import of the constructor, interval recurrence or three-row solver. Core
closure and real histories use coordinate sets; EVERY recorded SG entry is
recomputed by the original multiword C implementation. The infinite step is
the proved option correspondence in notes/14_outer_defect_algebra.md.
"""

import argparse
import hashlib
import json
from pathlib import Path

from verify_strip3_two_move import recompute


def neighbourhood(v):
    r, c = v
    return {(r, c), (r - 1, c), (r + 1, c), (r, c - 1), (r, c + 1)}


def board(width):
    return {(r, c) for c in range(1, width + 1) for r in (1, 2, 3)}


def tail(width, length):
    return ({(3, width + j) for j in range(1, length + 1)}
            | {(2, width + j) for j in range(2, length + 1, 2)})


def encode(cells):
    return format(sum(1 << (3 * (c - 1) + r - 1) for r, c in cells), "x")


def check(report):
    if report.get("format") != "nodekayles-rooted-comb-period-v1":
        raise ValueError("unknown certificate format")
    premises = {}
    summaries = []
    if not report.get("examples"):
        raise ValueError("empty certificate")
    names = set()
    for example in report["examples"]:
        if example["name"] in names:
            raise ValueError("duplicate example")
        names.add(example["name"])
        width, start = example["core_width"], example["start"]
        if (type(width) is not int or type(start) is not int
                or width < 1 or start < 2 or width + start + 9 > 128):
            raise ValueError("invalid width or start")
        if example["phase"] != 0 or example["period"] != 4 or example["window"] != 6:
            raise ValueError("these reachable grid certificates require phase 0 and window 6")
        if example["tail_history"] != "top row at core_width+1, core_width+3, ...":
            raise ValueError("unexpected tail history")
        cells = [tuple(v) for v in example["cells"]]
        core = frozenset(cells)
        root = tuple(example["root"])
        if len(core) != len(cells) or not core <= board(width):
            raise ValueError("invalid core coordinates")
        if root != (3, width) or cells[0] != root:
            raise ValueError("incorrect attachment root")
        edges = [[i, j] for i, v in enumerate(cells)
                 for j, u in enumerate(cells) if i < j and u in neighbourhood(v)]
        if example["edges"] != edges:
            raise ValueError("edges disagree with the grid geometry")
        initial = board(width)
        for move in example["prior_history"]:
            move = tuple(move)
            if move not in initial:
                raise ValueError("illegal prior history")
            initial -= neighbourhood(move)
        if initial - {(1, width)} != set(core):
            raise ValueError("prior history does not give the claimed core")

        # Build the complete set of root-preserving descendants independently.
        closure = {core}
        pending = [core]
        for position in pending:
            for move in position:
                if root not in neighbourhood(move):
                    rest = frozenset(position - neighbourhood(move))
                    if rest not in closure:
                        closure.add(rest)
                        pending.append(rest)
        states = {}
        for row in example["states"]:
            mask = int(row["mask"], 16)
            if mask <= 0 or mask >= 1 << len(cells) or not mask & 1 or mask in states:
                raise ValueError("invalid or duplicate core state")
            position = frozenset(v for i, v in enumerate(cells) if mask >> i & 1)
            if position not in closure:
                raise ValueError("state is not a root-preserving descendant")
            values = row["values"]
            if len(values) != start + 10 or any(type(v) is not int or v < 0 for v in values):
                raise ValueError("incorrect finite value table")
            if any(values[n] != values[n + 4] for n in range(start, start + 6)):
                raise ValueError("period-lifting window fails")
            states[mask] = values
            for n, expected in enumerate(values):
                attached = position | tail(width, n)
                links = {(v, u) for v in position for u in tail(width, n)
                         if u in neighbourhood(v)}
                if n and links != {(root, (3, width + 1))}:
                    raise ValueError("attachment is not a single edge")
                key = width + n, encode(attached)
                if key in premises and premises[key]["value"] != expected:
                    raise ValueError("conflicting value entries")
                premises[key] = {"family": "rooted-comb", "n": key[0],
                                 "mask": key[1], "value": expected}
        if len(states) != len(closure):
            raise ValueError("missing root-preserving descendant")
        full_mask = (1 << len(cells)) - 1
        values = states[full_mask]
        if example["full_core_period_block"] != values[start:start + 4]:
            raise ValueError("incorrect advertised period block")
        residues = sorted(n % 4 for n in range(start, start + 4) if not values[n])
        if example["full_core_zero_residues"] != residues:
            raise ValueError("incorrect advertised P residues")

        # All positive tail lengths are genuine positions reached from an empty
        # three-row board. Auxiliary core states need not have the same history.
        for n in range(1, start + 10):
            reachable = board(width + n)
            history = [tuple(v) for v in example["prior_history"]]
            history += [(1, width + j) for j in range(1, n + 1, 2)]
            for move in history:
                if move not in reachable:
                    raise ValueError("extended history contains an illegal move")
                reachable -= neighbourhood(move)
            if reachable != core | tail(width, n):
                raise ValueError("history does not create the advertised tail")
        summaries.append({"name": example["name"], "states": len(states),
                          "start": start, "P_residues": residues})
    return list(premises.values()), summaries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--audit", type=Path)
    parser.add_argument("--timeout", type=float, default=180)
    args = parser.parse_args()
    report = json.loads(args.certificate.read_text())
    premises, summaries = check(report)
    independent = recompute(premises, args.timeout)
    for summary in summaries:
        print(summary)
    print(f"original multiword C: all {len(premises)} distinct SG premises agree")
    print("closed core states and six-term windows verified; period lifting applies")
    if args.audit:
        audit = {
            "format": "nodekayles-rooted-comb-audit-v1",
            "certificate_sha256": hashlib.sha256(args.certificate.read_bytes()).hexdigest(),
            "verifier_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "original_c": independent, "examples": summaries,
            "conclusion": "finite premises of the period-4 lifting theorem independently checked",
        }
        args.audit.write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n")


if __name__ == "__main__":
    main()
