#!/usr/bin/env python3
"""A reachable witness family that defeats bounded-thickness component oracles.

This bounds complete memoized SG expansion with the specified reductions,
not the complexity of Node Kayles or of all possible solvers. See notes/27.
No SG search is run by this script.
"""

import argparse
import json

from strip3_history_language import outer_mask


def seed_count(length):
    """Regular outer words of this length whose first letter is top=1."""
    if type(length) is not int or length < 1:
        raise ValueError("positive seed length required")
    u, v, w = 0, 1, 0
    for _ in range(1, length):
        u, v, w = u + w, 2 * u + w, v
    return u + v + w


def seed_words(length):
    if type(length) is not int or length < 1:
        raise ValueError("positive seed length required")

    def visit(prefix, side, age):
        if len(prefix) == length:
            yield prefix
            return
        yield from visit(prefix + (0,), side, min(age + 1, 3))
        for bit in (1, 4):
            if age >= (2 if bit == side else 3):
                yield from visit(prefix + (bit,), bit, 1)
    yield from visit((1,), 1, 1)


def protected_word(seed, cutoff):
    if type(cutoff) is not int or cutoff < 2:
        raise ValueError("thickness cutoff must be at least 2")
    if not seed or seed[0] != 1:
        raise ValueError("seed must start with the top-left move")
    word = tuple(seed) + (0,) * (cutoff + 2)
    outer_mask(word)  # Validate the complete history.
    return word


def family(width, cutoff=5):
    if type(width) is not int or width < cutoff + 3:
        raise ValueError("width must be at least cutoff+3")
    for seed in seed_words(width - cutoff - 2):
        word = protected_word(seed, cutoff)
        yield word, outer_mask(word)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--width", type=int, default=24)
    parser.add_argument("--cutoff", type=int, default=5)
    args = parser.parse_args()
    if args.cutoff < 2 or args.width < args.cutoff + 3 or args.width > 10000:
        parser.error("2<=cutoff, cutoff+3<=width<=10000")
    length = args.width - args.cutoff - 2
    count = seed_count(length)
    print(json.dumps({"width": args.width, "oracle_max_full_run": args.cutoff,
                      "protected_full_columns": args.cutoff + 1,
                      "free_history_columns": length, "mandatory_distinct_states": count,
                      "mandatory_enumerated_moves_lower_bound": args.width * count,
                      "recurrence": "b1=1,b2=1,b3=2; bt=bt-1+bt-2+bt-3",
                      "scope": "complete memoized SG expansion; paths/cycles, rectangle symmetries, equal-component cancellation, and whole-component oracles with bounded full-column runs",
                      "not_claimed": "lower bound for all algorithms, unbounded-width rules, or P/N-only search"}))


if __name__ == "__main__":
    main()
