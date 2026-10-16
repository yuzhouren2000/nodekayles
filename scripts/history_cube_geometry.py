#!/usr/bin/env python3
"""Exact combinatorics of independent histories; no Node Kayles search.

An equality cube cannot contain two independent histories. Maximal-history
parity gives a separate sufficient SG=0/1 rule, computable by a strip transfer.
"""
import argparse
import json


def neighbourhoods(columns):
    vertices = range(3*columns)
    return [sum(1 << u for u in vertices
                if abs(v//3-u//3) + abs(v%3-u%3) == 1) for v in vertices]


def cube_witness(adjacency, yes, no):
    """Return two one-move-related members, or None for a singleton cube.

    The cube must contain an independent history. Its minimal member is yes;
    an unconstrained vertex not adjacent to yes gives the witnessed extension.
    """
    full = (1 << len(adjacency)) - 1
    if min(yes, no) < 0 or (yes | no) & ~full or yes & no:
        raise ValueError("invalid cube")
    blocked = yes
    for v, neighbours in enumerate(adjacency):
        if yes >> v & 1:
            if yes & neighbours:
                raise ValueError("cube has no independent histories")
            blocked |= neighbours
    free = full & ~(blocked | no)
    return (yes, yes | (free & -free)) if free else None


def independent_history_count(n):
    """Count histories using column-independent masks, with exact integers."""
    states = {0: 1}
    letters = (0, 1, 2, 4, 5)
    for _ in range(n):
        states = {b: sum(count for a, count in states.items() if not a & b) for b in letters}
    return sum(states.values())


def maximal_history_parities(mask):
    """Even/odd support of maximal independent sets in a three-row mask.

    The DP keeps the last two chosen columns. Once the next column is chosen,
    independence and domination of the middle column can both be finalized.
    Counts are Boolean parity support (bit 0 even, bit 1 odd), not signed sums.
    """
    if mask < 0:
        raise ValueError("negative mask")
    n = (mask.bit_length()+2)//3
    if not n:
        return 1
    active = [(mask >> (3*c)) & 7 for c in range(n)] + [0]
    choices = [[a for a in (0, 1, 2, 4, 5) if not a & ~col] for col in active]
    states = {(0, b): 1 << (b.bit_count() % 2) for b in choices[0]}
    for c in range(n):
        following = {}
        for (a, b), parities in states.items():
            for d in choices[c+1]:
                if b & d:
                    continue
                dominated = a | b | (b << 1) | (b >> 1) | d
                if active[c] & ~dominated:
                    continue
                shifted = ((parities & 1) << 1) | (parities >> 1) if d.bit_count() % 2 else parities
                key = b, d
                following[key] = following.get(key, 0) | shifted
        states = following
    result = 0
    for support in states.values():
        result |= support
    assert result in (1, 2, 3)
    return result


def parity_sg(mask):
    """Uniform maximal-history parity fixes SG exactly; mixed means unknown."""
    support = maximal_history_parities(mask)
    return {1: 0, 2: 1, 3: None}[support]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-n", type=int, default=100)
    args = parser.parse_args()
    if not 1 <= args.max_n <= 100:
        parser.error("require 1..100")
    for n in range(1, args.max_n+1):
        mask = (1 << (3*n))-1
        print(json.dumps(dict(n=n, histories=independent_history_count(n),
                              maximal_history_parities=maximal_history_parities(mask),
                              parity_SG=parity_sg(mask))))


if __name__ == "__main__":
    main()
