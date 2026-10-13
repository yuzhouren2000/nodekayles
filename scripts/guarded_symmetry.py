#!/usr/bin/env python3
"""Exact SG certificates from a guarded fixed set of an involution.

Every nonfixed vertex must be nonadjacent to its partner. Each fixed centre
has an odd number of fixed pendant leaves. Arbitrary edges between centres
and to the paired vertices are allowed. The SG is the parity of the number
of centres plus fixed isolated vertices. See notes/20_guarded_symmetry.md.
"""

import argparse
from dataclasses import asdict, dataclass
import json

from strip3_defect_algebra import validate_core


@dataclass(frozen=True)
class GuardedCertificate:
    centres: tuple
    leaves: tuple
    isolates: tuple

    @property
    def grundy(self):
        return (len(self.centres) + len(self.isolates)) % 2


def guarded_certificate(adjacency, involution):
    """Check a supplied symmetry; return None when the theorem does not apply.

    Invalid graphs, permutations or automorphisms raise ValueError. This
    checks the structural hypotheses only; it never evaluates game options.
    Empty fixed sets include the usual fixed-point-free pairing theorem.
    """
    graph, sigma = tuple(adjacency), tuple(involution)
    n = len(graph)
    if graph:
        validate_core(graph)
    if (len(sigma) != n or any(type(v) is not int for v in sigma)
            or set(sigma) != set(range(n))
            or any(sigma[sigma[v]] != v for v in range(n))):
        raise ValueError("an involutive permutation is required")
    for v, neighbours in enumerate(graph):
        mapped = sum(1 << sigma[u] for u in range(n) if neighbours >> u & 1)
        if mapped != graph[sigma[v]]:
            raise ValueError("the involution is not a graph automorphism")
    if any(sigma[v] != v and graph[v] >> sigma[v] & 1 for v in range(n)):
        return None

    fixed = {v for v in range(n) if sigma[v] == v}
    isolates = {v for v in fixed if not graph[v]}
    leaves = {v for v in fixed if graph[v].bit_count() == 1}
    # An isolated fixed K2 permits either endpoint as its centre. Choose the
    # lower label; all other pendant leaves have an unambiguous centre.
    for v in sorted(tuple(leaves)):
        u = graph[v].bit_length() - 1
        if u in leaves and v < u:
            leaves.remove(v)
    centres = fixed - isolates - leaves
    groups = tuple(tuple(sorted(v for v in leaves if graph[v] == 1 << u))
                   for u in sorted(centres))
    if any(not group or len(group) % 2 == 0 for group in groups):
        return None
    if set().union(set(), *map(set, groups)) != leaves:
        return None
    return GuardedCertificate(tuple(sorted(centres)), groups, tuple(sorted(isolates)))


def central_pair_history(width):
    """One-based (row,column) history on an odd empty three-row board; SG=1."""
    if type(width) is not int or width < 3 or width % 2 == 0:
        raise ValueError("odd width at least 3 required")
    centre = (width + 1) // 2
    return ((1, centre - 1), (1, centre + 1))


def four_move_p_history(width):
    """One-based legal four-move history on every width >=8; SG=0.

    The left three columns are the fixed SG=1 core 775. The two middle-row
    moves separate it from a symmetric SG=1 remainder. No claim is made
    that either player can force the opponent to follow this history.
    """
    if type(width) is not int or width < 8:
        raise ValueError("width at least 8 required")
    even = width - width % 2
    return ((2, 4), (1, even // 2 + 1), (1, even // 2 + 3), (2, even))


def remaining_columns(width, history):
    """Apply the history by coordinates, rejecting any illegal move."""
    if type(width) is not int or width < 1:
        raise ValueError("positive board width required")
    columns = [7] * width
    for row, column in history:
        if row not in (1, 2, 3) or not 1 <= column <= width:
            raise ValueError("move outside the board")
        bit, i = 1 << (row - 1), column - 1
        if not columns[i] & bit:
            raise ValueError("illegal move")
        columns[i] &= ~(bit | bit << 1 | bit >> 1)
        if i:
            columns[i - 1] &= ~bit
        if i + 1 < width:
            columns[i + 1] &= ~bit
    return tuple(columns)


def reflection_certificate(columns):
    """Recognise horizontal guarded symmetry of an actual three-row mask."""
    from kayles_boundary_transfer import column_graph
    columns = tuple(columns)
    graph, cells = column_graph(columns)
    indices = {cell: i for i, cell in enumerate(cells)}
    reflected = [(row, len(columns) - column - 1) for row, column in cells]
    if any(cell not in indices for cell in reflected):
        return None
    return guarded_certificate(graph, tuple(indices[cell] for cell in reflected))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("family", choices=("central", "four-move-p"))
    parser.add_argument("width", type=int)
    parser.add_argument("--check-geometry", action="store_true",
                        help="materialise the board and check its symmetry")
    args = parser.parse_args()
    if args.family == "central":
        history = central_pair_history(args.width)
        centre = (args.width + 1) // 2
        result = {"width": args.width, "history": history, "grundy": 1,
                  "winning_moves": ((2, centre), (3, centre)),
                  "playable_cells": 3 * args.width - (5 if args.width == 3 else 7)}
    else:
        history = four_move_p_history(args.width)
        result = {"width": args.width, "history": history, "grundy": 0,
                  "region_values": (1, 1, 0) if args.width % 2 else (1, 1),
                  "playable_cells": 3 * args.width - (12 if args.width < 10 else 16) - args.width % 2}
    if args.check_geometry:
        columns = remaining_columns(args.width, history)
        even = args.width - args.width % 2
        certificate = reflection_certificate(columns if args.family == "central" else columns[4:even - 1])
        if certificate is None or certificate.grundy != 1:
            raise AssertionError("family failed the structural certificate")
        result["symmetric_part"] = asdict(certificate)
    print(json.dumps(result, ensure_ascii=False))


if __name__ == "__main__":
    main()
