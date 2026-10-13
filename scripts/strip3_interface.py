#!/usr/bin/env python3
"""Find an outcome obstruction to a boundary-deletion SG summary.

Enumerates reachable, connected three-row positions with three full final
columns. Uses an exact subset recurrence, without the production solver.
This refutes the declared summary, not arbitrary finite interface methods.
"""

import argparse
from functools import cache
import json
from pathlib import Path


def evaluator(width):
    closed = []
    for v in range(3 * width):
        c, r = divmod(v, 3)
        closed.append(sum(1 << (3 * cc + rr)
                          for cc, rr in ((c, r), (c - 1, r), (c + 1, r),
                                         (c, r - 1), (c, r + 1))
                          if 0 <= cc < width and 0 <= rr < 3))

    @cache
    def value(mask):
        if not mask:
            return 0
        shift = ((mask & -mask).bit_length() - 1) // 3 * 3
        if shift:
            return value(mask >> shift)
        options = set()
        bits = mask
        while bits:
            bit = bits & -bits
            bits ^= bit
            options.add(value(mask & ~closed[bit.bit_length() - 1]))
        g = 0
        while g in options:
            g += 1
        return g

    return closed, value


def connected(mask, closed):
    reached = frontier = mask & -mask
    while frontier:
        neighbours = 0
        while frontier:
            bit = frontier & -frontier
            frontier ^= bit
            neighbours |= closed[bit.bit_length() - 1]
        frontier = neighbours & mask & ~reached
        reached |= frontier
    return bool(mask) and reached == mask


def reachable(width, closed):
    # An independent set is exactly a legal history, in any ordering.
    # A full last column excludes any occupied point in that column.
    def visit(column, previous, occupied, deleted):
        if column == width - 1:
            yield ((1 << (3 * width)) - 1) & ~deleted, occupied
            return
        for rows in (0, 1, 2, 4, 5):
            if rows & previous:
                continue
            bits = rows << (3 * column)
            blocked = 0
            while bits:
                bit = bits & -bits
                bits ^= bit
                blocked |= closed[bit.bit_length() - 1]
            yield from visit(column + 1, rows,
                             occupied | (rows << (3 * column)), deleted | blocked)

    yield from visit(0, 0, 0, 0)


def signature(mask, width, closed, value):
    boundary_shift = 3 * (width - 1)
    table = []
    for selected in range(8):
        blocked = 0
        for row in range(3):
            if selected >> row & 1:
                blocked |= closed[boundary_shift + row]
        table.append([value(mask & ~blocked & ~(deleted << boundary_shift))
                      for deleted in range(8)])
    return table


def coordinates(mask):
    return sorted([v % 3 + 1, v // 3 + 1]
                  for v in range(mask.bit_length()) if mask >> v & 1)


def certificate(width, pair, closed, value):
    positions = []
    for mask, history in pair:
        extension = mask | (7 << (3 * width))
        cut = (extension & ~closed[3 * 5 + 1]
               if width == 9 and mask >> (3 * 5 + 1) & 1 else None)
        record = {
            "history": coordinates(history), "mask": f"{mask:x}",
            "value": value(mask), "extension_mask": f"{extension:x}",
            "extension_value": value(extension),
            "extension_options": [
                {"move": move, "value": value(extension & ~closed[
                    3 * (move[1] - 1) + move[0] - 1])}
                for move in coordinates(extension)],
        }
        if cut is not None:
            base_cut = mask & ~closed[3 * 5 + 1]
            record["cut_at_2_6"] = {
                "left_value": value(cut & ((1 << 15) - 1)),
                "right_value": value(cut >> 18),
                "base_option_value": value(base_cut),
                "base_right_value": value(base_cut >> 18),
            }
        positions.append(record)
    return {
        "format": "nodekayles-strip3-interface-obstruction-v1", "n": width,
        "boundary": [[row, width] for row in (1, 2, 3)],
        "identical_full_suffix_columns": 3,
        "signature_definition": "entry[S][D] = SG(G minus (union of closed "
                                "neighbourhoods of S) minus D); row 1 is bit 0",
        "shared_signature": signature(pair[0][0], width, closed, value),
        "positions": positions,
        "scope": "Refutes a summary consisting of width, the last three columns, "
                 "and these 64 SG values. Does not refute all finite interface "
                 "summaries or periodicity of full strips.",
    }


def search(max_width):
    closed, value = evaluator(max_width + 1)
    for width in range(3, max_width + 1):
        seen, representatives = set(), {}
        suffix = ((1 << 9) - 1) << (3 * (width - 3))
        for mask, history in reachable(width, closed):
            if mask in seen or mask & suffix != suffix or not connected(mask, closed):
                continue
            seen.add(mask)
            table = signature(mask, width, closed, value)
            key = tuple(tuple(row) for row in table)
            extension_value = value(mask | (7 << (3 * width)))
            if key in representatives:
                other, other_history, other_value = representatives[key]
                if (other_value == 0) != (extension_value == 0):
                    return certificate(width, [(other, other_history), (mask, history)],
                                       closed, value)
            else:
                representatives[key] = mask, history, extension_value
        print(f"width={width}: {len(seen)} positions, "
              f"{len(representatives)} signatures", flush=True)
    raise ValueError("no outcome obstruction within the requested search")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-n", type=int, default=9)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not 3 <= args.max_n <= 11:
        parser.error("this small exact search supports 3 <= max-n <= 11")
    report = search(args.max_n)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    print(f"wrote {args.output}: extension values "
          f"{[p['extension_value'] for p in report['positions']]}")


if __name__ == "__main__":
    main()
