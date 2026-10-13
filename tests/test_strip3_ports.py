#!/usr/bin/env python3
"""Independent coordinate-level checks for the B/E/H port families.

Everything here is generated from one-based coordinates: cell sets, adjacency,
closed neighbourhoods, components and the reference SG.  The cutting formulas
under test are re-expressed as raw coordinate sets, so a mistake in the shared
mask construction in scripts/strip3_ports.py cannot hide here.

Checks:
  A  L1 geometry: for every legal middle move of B/E/H with n <= max_n the
     successor splits into the two declared pieces (left/right of the cut
     column), the pieces match the declared template sets exactly (including
     the reflection for the i>j left piece of B), and no edge crosses the cut.
  B  L1 values: with the independent recursion, for n <= --max-sg-n every
     middle move's successor value equals rec(left) xor rec(right), and the
     position value equals the mex over all legal moves.
  C  L2 audit: premises d(t)=0 / a(t)=0, witness-move legality, index algebra
     and the witness option value for the four extension relations.
  D  Reflection identities: B and H are horizontally symmetric as cell sets.

Run: python3 tests/test_strip3_ports.py [--max-n 12] [--max-sg-n 8]
"""

import argparse
from functools import lru_cache

# ------------------------------------------------- coordinates (one-based)

def full_cells(n):
    return frozenset((r, c) for r in (1, 2, 3) for c in range(1, n + 1))


def hole_set(n, j):
    return frozenset(
        (r, c) for r, c in ((1, j - 1), (1, j), (1, j + 1), (2, j))
        if 1 <= c <= n)


def cells_of(name, n, j=None):
    """Raw coordinate definitions, deliberately re-typed here."""
    base = full_cells(n)
    if name == "F":
        return base
    if name == "A":
        return base - {(2, n)} if n >= 1 else frozenset()
    if name == "D":
        return base - {(2, 1), (2, n)} if n >= 1 else frozenset()
    if name == "B":
        return base - hole_set(n, j)
    if name == "E":
        return base - hole_set(n, j) - {(2, 1)}
    if name == "H":
        return base - hole_set(n, j) - {(2, 1)} - {(2, n)}
    raise ValueError(name)


def move_targets(cell):
    r, c = cell
    return frozenset((rr, cc) for rr, cc in
                     ((r, c), (r - 1, c), (r + 1, c), (r, c - 1), (r, c + 1))
                     if 1 <= rr <= 3 and cc >= 1)


def reflect(cells, width):
    return frozenset((r, width + 1 - c) for r, c in cells)


def shift(cells, dc):
    return frozenset((r, c + dc) for r, c in cells)


def edge_between(left, right):
    for (r1, c1) in left:
        for (r2, c2) in right:
            if abs(r1 - r2) + abs(c1 - c2) == 1:
                return True
    return False


# --------------------------------------------------- independent recursion

def closed_neighbourhoods(cells):
    return {v: frozenset(u for u in move_targets(v) if u in cells)
            for v in cells}


def component_split(cells):
    remaining = set(cells)
    comps = []
    while remaining:
        start = min(remaining, key=lambda x: (x[1], x[0]))
        comp = {start}
        frontier = {start}
        while frontier:
            nxt = set()
            for v in frontier:
                for u in move_targets(v):
                    if u in remaining and u not in comp:
                        nxt.add(u)
            frontier = nxt
            comp |= frontier
        remaining -= comp
        comps.append(frozenset(comp))
    return comps


@lru_cache(maxsize=None)
def rec_value(cells):
    if not cells:
        return 0
    comps = component_split(cells)
    if len(comps) > 1:
        g = 0
        for comp in comps:
            g ^= rec_value(comp)
        return g
    closed = closed_neighbourhoods(cells)
    seen = {rec_value(frozenset(cells - closed[v])) for v in cells}
    g = 0
    while g in seen:
        g += 1
    return g


def mex(values):
    g = 0
    while g in values:
        g += 1
    return g


# ------------------------------------------------------------- check A / B

def middle_legal(name, n, j, i):
    """Declared legality rule for a middle move at column i."""
    if i < 1 or i > n:
        return False
    if name == "B":
        return i != j
    if name == "E":
        return i not in (1, j)
    return i not in (1, n, j)


def expected_pieces(name, n, j, i):
    """Declared L1 template for the two pieces, as coordinate sets."""
    if name == "B":
        if i < j:
            left = cells_of("A", i - 1)
            right = shift(cells_of("E", n - i, j - i), i) if n - i >= 1 else frozenset()
        else:
            left = reflect(cells_of("E", i - 1, i - j), i - 1) if i - 1 >= 1 else frozenset()
            right = shift(reflect(cells_of("A", n - i), n - i), i) if n - i >= 1 else frozenset()
    elif name == "E":
        if i < j:
            left = cells_of("D", i - 1)
            right = shift(cells_of("E", n - i, j - i), i) if n - i >= 1 else frozenset()
        else:
            left = cells_of("H", i - 1, j) if i - 1 >= 1 else frozenset()
            right = shift(reflect(cells_of("A", n - i), n - i), i) if n - i >= 1 else frozenset()
    else:
        if i < j:
            left = cells_of("D", i - 1)
            right = shift(cells_of("H", n - i, j - i), i) if n - i >= 1 else frozenset()
        else:
            left = cells_of("H", i - 1, j) if i - 1 >= 1 else frozenset()
            right = shift(cells_of("D", n - i), i) if n - i >= 1 else frozenset()
    return left, right


def check_geometry(max_n):
    moves = 0
    for n in range(1, max_n + 1):
        for j in range(1, n + 1):
            for name in ("B", "E", "H"):
                base = cells_of(name, n, j)
                for i in range(1, n + 1):
                    present = (2, i) in base
                    if present != middle_legal(name, n, j, i):
                        raise AssertionError(f"legality mismatch {name}({n},{j}) i={i}")
                    if not present:
                        continue
                    moves += 1
                    succ = base - move_targets((2, i))
                    left = frozenset(c for c in succ if c[1] < i)
                    right = frozenset(c for c in succ if c[1] > i)
                    if left | right != succ:
                        raise AssertionError(f"cut column not empty {name}({n},{j}) i={i}")
                    if edge_between(left, right):
                        raise AssertionError(f"cross edge {name}({n},{j}) i={i}")
                    exp_left, exp_right = expected_pieces(name, n, j, i)
                    if left != exp_left:
                        raise AssertionError(
                            f"left piece {name}({n},{j}) i={i}: "
                            f"diff {sorted(left ^ exp_left)}")
                    if right != exp_right:
                        raise AssertionError(
                            f"right piece {name}({n},{j}) i={i}: "
                            f"diff {sorted(right ^ exp_right)}")
    return moves


def check_values(max_n):
    for n in range(1, max_n + 1):
        for j in range(1, n + 1):
            for name in ("B", "E", "H"):
                base = cells_of(name, n, j)
                options = {rec_value(frozenset(base - move_targets(v))) for v in base}
                value = rec_value(base)
                if value != mex(options):
                    raise AssertionError(f"mex mismatch {name}({n},{j})")
                for i in range(1, n + 1):
                    if not middle_legal(name, n, j, i):
                        continue
                    succ = frozenset(base - move_targets((2, i)))
                    left, right = expected_pieces(name, n, j, i)
                    if rec_value(succ) != (rec_value(left) ^ rec_value(right)):
                        raise AssertionError(
                            f"L1 value mismatch {name}({n},{j}) i={i}")
    return max_n


# --------------------------------------------------------------- check C

def check_l2(max_n, zero_d, zero_a):
    """Audit the four extension relations: premises, witness legality, indices.

    zero_d / zero_a are the computed zero tails (from A/D sequences).
    """
    checks = 0
    for t in zero_d:
        q = t + 1
        for n in range(1, max_n + 1):
            for j in range(1, n + 1):
                # e(n+q, j+q) != e(n, j) via middle move at column q
                big = cells_of("E", n + q, j + q)
                if (2, q) not in big or q == 1 or q == j + q:
                    raise AssertionError(f"L2 e-witness illegal at t={t},n={n},j={j}")
                succ = frozenset(big - move_targets((2, q)))
                left = frozenset(c for c in succ if c[1] < q)
                right = frozenset(c for c in succ if c[1] > q)
                if left != cells_of("D", q - 1):
                    raise AssertionError(f"L2 e-witness left t={t},n={n},j={j}")
                if right != shift(cells_of("E", n, j), q):
                    raise AssertionError(f"L2 e-witness right t={t},n={n},j={j}")
                if rec_value(left) != 0 or rec_value(big) == rec_value(right):
                    raise AssertionError(f"L2 e-relation fails t={t},n={n},j={j}")
                # h(n+q, j+q) != h(n, j) via the same cut
                big_h = cells_of("H", n + q, j + q)
                if (2, q) not in big_h or q in (1, n + q, j + q):
                    raise AssertionError(f"L2 h-witness illegal at t={t},n={n},j={j}")
                succ_h = frozenset(big_h - move_targets((2, q)))
                left_h = frozenset(c for c in succ_h if c[1] < q)
                right_h = frozenset(c for c in succ_h if c[1] > q)
                if left_h != cells_of("D", q - 1):
                    raise AssertionError(f"L2 h-witness left t={t},n={n},j={j}")
                if right_h != shift(cells_of("H", n, j), q):
                    raise AssertionError(f"L2 h-witness right t={t},n={n},j={j}")
                if rec_value(left_h) != 0 or rec_value(big_h) == rec_value(right_h):
                    raise AssertionError(f"L2 h-relation fails t={t},n={n},j={j}")
                # h(n+q, j) != h(n, j) via middle move at column n+1 (right cut)
                big_r = cells_of("H", n + q, j)
                if (2, n + 1) not in big_r or n + 1 in (1, n + q, j):
                    raise AssertionError(f"L2 h-right witness illegal t={t},n={n},j={j}")
                succ_r = frozenset(big_r - move_targets((2, n + 1)))
                left_r = frozenset(c for c in succ_r if c[1] < n + 1)
                right_r = frozenset(c for c in succ_r if c[1] > n + 1)
                if left_r != cells_of("H", n, j):
                    raise AssertionError(f"L2 h-right left t={t},n={n},j={j}")
                if right_r != shift(cells_of("D", q - 1), n + 1):
                    raise AssertionError(f"L2 h-right right t={t},n={n},j={j}")
                if rec_value(right_r) != 0 or rec_value(big_r) == rec_value(left_r):
                    raise AssertionError(f"L2 h-right relation fails t={t},n={n},j={j}")
                checks += 3
    for t in zero_a:
        for n in range(1, max_n + 1):
            for j in range(1, n + 1):
                big = cells_of("E", n + t + 1, j)
                if (2, n + 1) not in big or n + 1 in (1, j):
                    raise AssertionError(f"L2 eh-witness illegal t={t},n={n},j={j}")
                succ = frozenset(big - move_targets((2, n + 1)))
                left = frozenset(c for c in succ if c[1] < n + 1)
                right = frozenset(c for c in succ if c[1] > n + 1)
                if left != cells_of("H", n, j):
                    raise AssertionError(f"L2 eh left t={t},n={n},j={j}")
                if right != shift(reflect(cells_of("A", t), t), n + 1):
                    raise AssertionError(f"L2 eh right t={t},n={n},j={j}")
                if rec_value(right) != 0 or rec_value(big) == rec_value(left):
                    raise AssertionError(f"L2 eh relation fails t={t},n={n},j={j}")
                checks += 1
    return checks


# --------------------------------------------------------------- check D

def check_reflection(max_n):
    for n in range(1, max_n + 1):
        for j in range(1, n + 1):
            if cells_of("B", n, j) != reflect(cells_of("B", n, n + 1 - j), n):
                raise AssertionError(f"B reflection fails n={n} j={j}")
            if cells_of("H", n, j) != reflect(cells_of("H", n, n + 1 - j), n):
                raise AssertionError(f"H reflection fails n={n} j={j}")
            # E maps to the right-gap family under reflection; check that too.
            right_gap = full_cells(n) - hole_set(n, n + 1 - j) - {(2, n)}
            if reflect(cells_of("E", n, j), n) != right_gap:
                raise AssertionError(f"E right-gap image fails n={n} j={j}")
    return max_n


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--max-n", type=int, default=12,
                        help="parameter range for geometry checks")
    parser.add_argument("--max-sg-n", type=int, default=8,
                        help="parameter range for the slow independent SG checks")
    args = parser.parse_args()

    moves = check_geometry(args.max_n)
    print(f"geometry: {moves} legal middle moves over B/E/H with n<={args.max_n} "
          f"split exactly into the declared pieces")
    check_values(args.max_sg_n)
    print(f"values: independent recursion agrees with mex and with the piece xor "
          f"for all B/E/H with n<={args.max_sg_n}")
    # zero tails from the independently verified sequences
    seq_d = [0, 0, 0, 1, 1, 2, 0, 3, 1, 1, 0, 5, 3, 4, 2, 2, 0, 1]
    seq_a = [0, 0, 3, 1, 2, 0, 0, 1, 5, 2, 6, 3, 3, 2, 0]
    zero_d = [t for t in range(1, len(seq_d)) if seq_d[t] == 0]
    zero_a = [t for t in range(len(seq_a)) if seq_a[t] == 0]
    # cross-check the tails lazily for the range used in the L2 audit
    for t in zero_d:
        if t + 1 <= args.max_sg_n + 2:
            if rec_value(cells_of("D", t)) != 0:
                raise AssertionError(f"D({t}) should be zero")
    for t in zero_a:
        if t <= args.max_sg_n + 2:
            if rec_value(cells_of("A", t)) != 0:
                raise AssertionError(f"A({t}) should be zero")
    checks = check_l2(args.max_sg_n, [t for t in zero_d if t <= 6],
                      [t for t in zero_a if t <= 5])
    print(f"L2: {checks} extension audits passed (witness cells, piece shape, "
          f"zero tail, inequality) for the small zero tails")
    check_reflection(args.max_n)
    print(f"reflection: B and H are horizontally symmetric as cell sets, "
          f"E maps to the right-gap family, n<={args.max_n}")
    print("all port-family coordinate checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
