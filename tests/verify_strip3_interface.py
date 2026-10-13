#!/usr/bin/env python3
"""Independently recompute the boundary-summary obstruction.

Uses coordinate-set recursion with no symmetry or component reductions.
Also checks all distinct masks with the original multiword C implementation.
Imports neither the search script nor the specialised three-row solver.
"""

import argparse
from functools import cache
import hashlib
import json
from pathlib import Path

from verify_strip3_two_move import recompute


def full(width):
    return frozenset((row, column) for row in (1, 2, 3)
                     for column in range(1, width + 1))


def neighbourhood(vertex):
    row, column = vertex
    return frozenset(((row, column), (row - 1, column), (row + 1, column),
                      (row, column - 1), (row, column + 1)))


def encode(cells):
    return format(sum(1 << (3 * (column - 1) + row - 1)
                      for row, column in cells), "x")


@cache
def value(cells):
    options = {value(cells - neighbourhood(vertex)) for vertex in cells}
    g = 0
    while g in options:
        g += 1
    return g


def check(report):
    if report.get("format") != "nodekayles-strip3-interface-obstruction-v1":
        raise ValueError("unknown certificate")
    width = report["n"]
    if type(width) is not int or not 3 <= width <= 11:
        raise ValueError("unsupported certificate width")
    if report["boundary"] != [[row, width] for row in (1, 2, 3)]:
        raise ValueError("incorrect boundary")
    if report["identical_full_suffix_columns"] != 3:
        raise ValueError("incorrect suffix specification")
    table = report["shared_signature"]
    if len(table) != 8 or any(len(row) != 8 for row in table):
        raise ValueError("expected an 8 by 8 table")
    if any(type(g) is not int or g < 0 for row in table for g in row):
        raise ValueError("invalid SG value")
    if len(report["positions"]) != 2:
        raise ValueError("expected two positions")

    premises = {}

    def verify(cells, expected, board_width=width):
        if type(expected) is not int or expected < 0 or value(cells) != expected:
            raise ValueError("incorrect SG value")
        key = board_width, encode(cells)
        premises[key] = {"family": "interface", "n": board_width,
                         "mask": key[1], "value": expected}

    for position in report["positions"]:
        cells = full(width)
        for move in position["history"]:
            vertex = tuple(move)
            if vertex not in cells:
                raise ValueError("history contains an illegal move")
            cells -= neighbourhood(vertex)
        if encode(cells) != position["mask"]:
            raise ValueError("history does not yield the recorded mask")
        seen, frontier = set(), {min(cells)}
        while frontier:
            seen |= frontier
            frontier = set().union(*(neighbourhood(v) for v in frontier)) & cells - seen
        if seen != cells:
            raise ValueError("position is disconnected")
        if not frozenset((r, c) for r in (1, 2, 3)
                         for c in range(width - 2, width + 1)) <= cells:
            raise ValueError("the last three columns are not full")
        verify(cells, position["value"])

        for selected in range(8):
            blocked = frozenset().union(*(
                neighbourhood((row, width)) for row in (1, 2, 3)
                if selected >> (row - 1) & 1))
            for deleted in range(8):
                boundary_deleted = frozenset((row, width) for row in (1, 2, 3)
                                            if deleted >> (row - 1) & 1)
                verify(cells - blocked - boundary_deleted, table[selected][deleted])

        extension = cells | (full(width + 1) - full(width))
        if encode(extension) != position["extension_mask"]:
            raise ValueError("incorrect extension mask")
        verify(extension, position["extension_value"], width + 1)
        options = position["extension_options"]
        moves = [tuple(option["move"]) for option in options]
        if len(moves) != len(set(moves)) or set(moves) != extension:
            raise ValueError("option records do not cover every legal move")
        for option in options:
            verify(extension - neighbourhood(tuple(option["move"])),
                   option["value"], width + 1)
        g = 0
        while g in {option["value"] for option in options}:
            g += 1
        if g != position["extension_value"]:
            raise ValueError("incorrect mex")

        if "cut_at_2_6" in position:
            if width != 9 or (2, 6) not in cells:
                raise ValueError("this diagnostic is specific to width nine")
            cut = extension - neighbourhood((2, 6))
            left = frozenset(v for v in cut if v[1] < 6)
            right = frozenset(v for v in cut if v[1] > 6)
            if left | right != cut or any(neighbourhood(v) & right for v in left):
                raise ValueError("diagnostic is not a disjoint split")
            verify(left, position["cut_at_2_6"]["left_value"], width + 1)
            verify(right, position["cut_at_2_6"]["right_value"], width + 1)
            base_cut = cells - neighbourhood((2, 6))
            verify(base_cut, position["cut_at_2_6"]["base_option_value"])
            verify(frozenset(v for v in base_cut if v[1] > 6),
                   position["cut_at_2_6"]["base_right_value"])

    outcomes = [position["extension_value"] == 0 for position in report["positions"]]
    if outcomes[0] == outcomes[1]:
        raise ValueError("the two extended outcomes do not differ")
    return list(premises.values())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--audit", type=Path)
    args = parser.parse_args()
    report = json.loads(args.certificate.read_text())
    premises = check(report)
    independent = recompute(premises, timeout=60)
    print(f"legal histories, connected masks, common suffix, "
          f"128 table entries, and all extension moves checked")
    print(f"coordinate recursion: {value.cache_info().currsize} subsets; "
          f"original multiword C: {len(premises)} distinct positions agree")
    print(f"same declared summary, different extension outcomes: "
          f"{[p['extension_value'] for p in report['positions']]}")
    if args.audit:
        audit = {
            "format": "nodekayles-strip3-interface-audit-v1",
            "certificate_sha256": hashlib.sha256(args.certificate.read_bytes()).hexdigest(),
            "verifier_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "coordinate_recursion_subsets": value.cache_info().currsize,
            "original_c": independent,
            "conclusion": "the declared boundary SG table does not preserve column extension",
        }
        args.audit.write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n")


if __name__ == "__main__":
    main()
