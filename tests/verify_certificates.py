"""Operational verification of pairing certificates by strategy simulation.

For a candidate certificate phi(i,j) = (tau(i), sigma(j)) this test simulates
the second player's mirror strategy against EVERY legal first move x:

  * the partner phi(x) must be playable after the opponent's move at x,
  * the set removed by the first move and the response must be phi-invariant
    (so the invariant is restored: the mirror strategy can continue).

Unlike scripts/pairing_certificates.py (which checks the algebraic conditions
on orbit representatives), this walks the strategy move by move on the whole
board - an independent implementation path for the same statement.

Run:  python3 tests/verify_certificates.py
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def closed_nbhd(i, j, m, n):
    s = {(i, j)}
    for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        a, b = i + di, j + dj
        if 0 <= a < m and 0 <= b < n:
            s.add((a, b))
    return s


def simulate(m, n, tau, sigma):
    """True iff the mirror strategy survives every legal first move."""
    phi = lambda x: (tau[x[0]], sigma[x[1]])
    # (1) no fixed cell, (2) no adjacency
    for i in range(m):
        for j in range(n):
            x = (i, j)
            y = phi(x)
            if y == x:
                return False
            if abs(y[0] - i) + abs(y[1] - j) < 2:
                return False
    # (3) strategy simulation over all legal first moves
    for i in range(m):
        for j in range(n):
            x = (i, j)
            after = closed_nbhd(i, j, m, n)  # cells removed by the first move
            y = phi(x)
            if y in after:
                return False  # partner not playable
            removed = after | closed_nbhd(y[0], y[1], m, n)
            if {phi(c) for c in removed} != removed:
                return False  # invariant not restored
    return True


def main():
    cases = [
        # (m, n, tau, sigma) as 0-based tuples; expected result
        (4, 5, (2, 3, 0, 1), (0, 1, 2, 3, 4), True),   # the 4 x n row pairing
        (4, 7, (2, 3, 0, 1), (0, 1, 2, 3, 4, 5, 6), True),
        (4, 14, (2, 3, 0, 1), tuple(range(14)), True),
        (2, 6, (1, 0), (5, 4, 3, 2, 1, 0), True),      # 2 x even rotation
        (2, 5, (1, 0), (4, 3, 2, 1, 0), False),        # 2 x odd: must fail
        (1, 4, (0,), (2, 3, 0, 1), True),              # 1 x 4 column pairing
        (3, 8, (2, 1, 0), (7, 6, 5, 4, 3, 2, 1, 0), False),  # 3 x 8: no cert (G=2)
        (6, 8, (5, 4, 3, 2, 1, 0), (7, 6, 5, 4, 3, 2, 1, 0), True),  # rotation
    ]
    fails = 0
    for m, n, tau, sigma, want in cases:
        got = simulate(m, n, tau, sigma)
        status = "ok  " if got == want else "FAIL"
        if got != want:
            fails += 1
        viol = "" if want else "no"
        print(f"  {status} {m}x{n}: strategy simulation says "
              f"{'certificate valid' if got else 'no certificate'}"
              f" (expected {'valid' if want else 'invalid/none'})")
    print()
    if fails:
        print(f"{fails} FAILURE(S)")
        return 1
    print("ALL CERTIFICATE SIMULATIONS AS EXPECTED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
