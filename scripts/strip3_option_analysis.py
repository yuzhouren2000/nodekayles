#!/usr/bin/env python3
"""Independent structural analysis of the full-board option spectrum of 3xn
Node Kayles.

Deliberately independent of every C solver in ``src/cfast``: a plain bitmask +
connected-component memoised SG recursion over induced subgraphs.  It is used
only to test structural hypotheses about the option values

    opt_F(n, c, r) = SG( F_n minus the closed neighbourhood of cell (r, c) ),

with rows r in {top, mid, bot} and columns c = 1..n.  ``A(k)`` is the single
right-end gap ``F_k minus {(2, k)}`` and ``a(k) = SG(A(k))``.

Findings (documented in notes/40_option_spectrum_structure.md):

* middle row is separable over GF(2):
      opt_F(n, c, mid) = a(c-1) xor a(n-c).
  This is the known F/A/D middle-split theorem (notes/08, notes/10 L1); here
  it is reconfirmed with zero violations for all n <= N.
* edge row is NOT separable: there are no sequences p, q with
      opt_F(n, c, top) = p(c-1) xor q(n-c).
  A function on the interior grid is of that form iff every 2x2 "rectangle"
  has even parity; the script checks every admissible rectangle and many fail.
  So the option spectrum does not collapse to universal port sequences, and
  the middle columns alone do not determine G(F_n).

Run:  python3 scripts/strip3_option_analysis.py --sizes 16
"""

import argparse
from functools import lru_cache


def solve(n):
    """Return (G, top, mid, bot, A) for the full 3xn board."""

    def idx(r, c):          # r = 0,1,2 ; c = 1..n ; column-major bit index
        return 3 * (c - 1) + r

    closed = []
    for c in range(1, n + 1):
        for r in range(3):
            b = 0
            for dr, dc in ((0, 0), (1, 0), (-1, 0), (0, 1), (0, -1)):
                rr, cc = r + dr, c + dc
                if 0 <= rr < 3 and 1 <= cc <= n:
                    b |= 1 << idx(rr, cc)
            closed.append(b)

    memo = {0: 0}

    def parts(mask):
        out, rem = [], mask
        while rem:
            comp, frontier = 0, rem & -rem
            while frontier:
                bit = frontier & -frontier
                frontier ^= bit
                v = bit.bit_length() - 1
                comp |= bit
                frontier |= (closed[v] & rem) & ~comp
            out.append(comp)
            rem &= ~comp
        return out

    def sg(mask):
        if mask in memo:
            return memo[mask]
        cs = parts(mask)
        if len(cs) > 1:
            x = 0
            for cc in cs:
                x ^= sg(cc)
            memo[mask] = x
            return x
        seen, mm = set(), mask
        while mm:
            bit = mm & -mm
            mm ^= bit
            seen.add(sg(mask & ~closed[bit.bit_length() - 1]))
        g = 0
        while g in seen:
            g += 1
        memo[mask] = g
        return g

    full = (1 << (3 * n)) - 1
    top, mid, bot = [], [], []
    for c in range(1, n + 1):
        for r in range(3):
            value = sg(full & ~closed[idx(r, c)])
            (top if r == 0 else mid if r == 1 else bot).append(value)
    g = sg(full)
    gap = sg(full & ~(1 << idx(1, n)))          # A(n): drop the middle end cell
    return g, top, mid, bot, gap


def separable_violations(grid):
    """Count 2x2 parity violations of a dict {(i,j): value}."""
    keys = sorted(grid)
    violations = 0
    for (i, j) in keys:
        for (i2, j2) in keys:
            if i2 != i and j2 != j and (i, j2) in grid and (i2, j) in grid:
                if grid[(i, j)] ^ grid[(i, j2)] ^ grid[(i2, j)] ^ grid[(i2, j2)]:
                    violations += 1
    return violations, len(keys)


KNOWN_G = {1: 2, 2: 1, 3: 1, 4: 0, 5: 3, 6: 3, 7: 2, 8: 2, 9: 2, 10: 3,
           11: 3, 12: 5, 13: 2, 14: 4, 15: 1, 16: 3}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sizes", type=int, default=16)
    args = ap.parse_args()

    a, top, mid = {0: 0}, {}, {}
    for n in range(1, args.sizes + 1):
        g, tp, md, bt, gap = solve(n)
        a[n] = gap
        top[n], mid[n] = tp, md
        if n in KNOWN_G and g != KNOWN_G[n]:
            raise SystemExit(f"G(F_{n})={g} disagrees with the published {KNOWN_G[n]}")
    print("G(F_n) reproduces every published value for n <= %d"
          % min(args.sizes, max(KNOWN_G)))
    print("a(0..%d) =" % args.sizes, [a[k] for k in range(args.sizes + 1)])

    Wmid, Wtop = {}, {}
    for n in range(1, args.sizes + 1):
        for c in range(1, n + 1):
            i, j = c - 1, n - c
            if i >= 1 and j >= 1:               # interior columns only
                Wmid[(i, j)] = mid[n][c - 1]
                Wtop[(i, j)] = top[n][c - 1]

    vmid, nmid = separable_violations(Wmid)
    vtop, ntop = separable_violations(Wtop)
    print("middle row, separability 2x2 violations: %d / %d" % (vmid, nmid))
    print("edge   row, separability 2x2 violations: %d / %d" % (vtop, ntop))
    if vmid != 0 or vtop == 0:
        raise SystemExit("unexpected: middle should be separable and edge not")

    def mex(s):
        s = set(s); k = 0
        while k in s:
            k += 1
        return k

    bad = [n for n in range(1, args.sizes + 1)
           if mex([a[c - 1] ^ a[n - c] for c in range(1, n + 1)]) != KNOWN_G.get(n)]
    print("mex(middle options) != G(F_n) for n in", bad)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
