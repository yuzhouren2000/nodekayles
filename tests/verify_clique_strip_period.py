#!/usr/bin/env python3
"""Independently verify simultaneous period premises and reachable histories.

No imports from the constructor or the interval/comb solvers. All recorded
graph values are recomputed by the original multiword C graph solver.
"""

import argparse
import hashlib
import json
from pathlib import Path

from verify_strip3_two_move import recompute


def integer(value, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError("invalid nonnegative integer")
    return value


def cells(word, n, phase=0, left=0, right=0, offset=0):
    result = set()
    for i in range(n):
        mask = word[(phase + i) % len(word)]
        if i == 0:
            mask &= ~left
        if i == n - 1:
            mask &= ~right
        result.update((r + 1, offset + i + 1) for r in range(3) if mask >> r & 1)
    return result


def neighbours(v, u):
    return abs(v[0] - u[0]) + abs(v[1] - u[1]) <= 1


def insert(premises, position, width, value):
    # Empty strips use a one-column bounding box in the independent C input.
    width = max(width, 1)
    encoded = format(sum(1 << (3 * (c - 1) + r - 1) for r, c in position), "x")
    key = width, encoded
    if key in premises and premises[key]["value"] != value:
        raise ValueError("inconsistent finite SG premises")
    premises[key] = {"family": "clique-column", "n": width, "mask": encoded, "value": value}


def check_window(rows, start, period, size):
    integer(start, 1)
    integer(period, 1)
    end = 2 * start + 2 * period - 1
    if len(rows) != end + 1:
        raise ValueError("incorrect finite period-window length")
    for row in rows:
        if len(row) != size:
            raise ValueError("missing simultaneous interface states")
        for value in row:
            integer(value)
    if any(rows[n] != rows[n + period] for n in range(start, 2 * start + period)):
        raise ValueError("simultaneous period comparison fails")
    return end


def history(head, n):
    # One-based offsets. All same-row distances are >=2. Opposite rows do
    # not attack one another; the entire finite list is a legal history.
    return [(r, head + c) for c in range(1, n + 1)
            for r, residue in ((1, 1), (1, 3), (3, 6), (3, 8)) if c % 10 == residue]


def after_history(head, n):
    result = {(r, c) for r in (1, 2, 3) for c in range(1, head + n + 1)}
    for move in history(head, n):
        if move not in result:
            raise ValueError("illegal purported history")
        result = {v for v in result if not neighbours(v, move)}
    return result


def check(report):
    if report.get("format") != "nodekayles-clique-strip-period-suite-v1":
        raise ValueError("unknown period certificate format")
    premises, strips = {}, {}
    for item in report["strips"]:
        word = tuple(item["word"])
        if not word or any(type(c) is not int or c not in (0, 1, 2, 3, 4, 6) for c in word):
            raise ValueError("a column is not a clique")
        if word in strips or item["period"] % len(word):
            raise ValueError("duplicate word or phase-changing period")
        end = check_window(item["rows"], item["start"], item["period"], len(word) * 64)
        if item["maximum_length"] != end:
            raise ValueError("incorrect maximum length")
        digest = hashlib.sha256(json.dumps(item["rows"], separators=(",", ":")).encode()).hexdigest()
        if item["rows_sha256"] != digest:
            raise ValueError("table digest differs")
        blocks = [[item["rows"][n][phase * 64] for n in
                   range(item["start"], item["start"] + item["period"])]
                  for phase in range(len(word))]
        if item["full_period_blocks"] != blocks:
            raise ValueError("wrong advertised full period block")
        strips[word] = item
        for n, row in enumerate(item["rows"]):
            for phase in range(len(word)):
                for left in range(8):
                    for right in range(8):
                        insert(premises, cells(word, n, phase, left, right), n,
                               row[phase * 64 + left * 8 + right])
    if not strips:
        raise ValueError("no pure-strip certificate")
    cores = []
    for item in report["cores"]:
        word = tuple(item["word"])
        if word not in strips or word != tuple(map(int, "4646313136")):
            raise ValueError("unsupported reachable-history tail")
        if item["attachment"] != 4:
            raise ValueError("wrong attachment row")
        pure = strips[word]
        period, start = item["period"], item["start"]
        if period != pure["period"] or start < pure["start"]:
            raise ValueError("core and pure tail need a common window")
        head = integer(item["head"], 1)
        vertices = [tuple(v) for v in item["cells"]]
        expected = {(r, c) for r in (1, 2, 3) for c in range(1, head + 1)} - {(1, head)}
        root = (3, head)
        if len(vertices) != len(expected) or set(vertices) != expected or vertices[0] != root:
            raise ValueError("wrong rooted grid head")
        adjacency = [sum(1 << j for j, u in enumerate(vertices) if u != v and neighbours(v, u))
                     for v in vertices]
        if item["adjacency"] != adjacency:
            raise ValueError("adjacency disagrees with coordinates")
        # Independently close under every move whose closed neighbourhood
        # does not delete the attachment vertex.
        closure, pending = {frozenset(expected)}, [frozenset(expected)]
        for position in pending:
            for move in position:
                if not neighbours(root, move):
                    rest = frozenset(v for v in position if not neighbours(v, move))
                    if rest not in closure:
                        closure.add(rest)
                        pending.append(rest)
        states = {}
        for state in item["states"]:
            mask = int(state["mask"], 16)
            if not 0 < mask < 1 << len(vertices) or mask in states:
                raise ValueError("invalid or duplicate core mask")
            position = frozenset(v for i, v in enumerate(vertices) if mask >> i & 1)
            if position not in closure:
                raise ValueError("not a root-preserving core descendant")
            end = check_window(state["rows"], start, period, 8)
            if item["maximum_length"] != end:
                raise ValueError("wrong core window end")
            states[mask] = state["rows"]
            for n, row in enumerate(state["rows"]):
                for right in range(8):
                    tail = cells(word, n, right=right, offset=head)
                    edges = {(v, u) for v in position for u in tail if neighbours(v, u)}
                    desired = {(root, (3, head + 1))} if (3, head + 1) in tail else set()
                    if edges != desired:
                        raise ValueError("attachment is not a single edge")
                    insert(premises, position | tail, head + n, row[right])
        if len(states) != len(closure):
            raise ValueError("missing safe core descendants")
        full = (1 << len(vertices)) - 1
        block = [states[full][n][0] for n in range(start, start + period)]
        if item["full_period_block"] != block:
            raise ValueError("wrong full-core period block")
        if item["reachable_length_residues"] != [1, 2, 3, 4, 6, 7, 8, 9]:
            raise ValueError("right-boundary history condition lost")
        for n in range(1, item["maximum_length"] + 1):
            if n % 10 not in (0, 5):
                if after_history(head, n) != expected | cells(word, n, offset=head):
                    raise ValueError("history does not produce the claimed graph")
        cores.append({"head": head, "core_states": len(states), "start": start,
                      "period": period,
                      "zero_residues": sorted(n % period for n in range(start, start + period)
                                              if states[full][n][0] == 0)})
    return premises, {"pure_words": len(strips), "core_examples": cores}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--timeout", type=int, default=120)
    args = parser.parse_args()
    premises, summary = check(json.loads(args.certificate.read_text()))
    metadata = recompute(list(premises.values()), args.timeout)
    report = {"certificate": str(args.certificate), "certificate_sha256":
              hashlib.sha256(args.certificate.read_bytes()).hexdigest(),
              **summary, "different_sg_premises": len(premises), "recomputation": metadata}
    if args.output:
        args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
