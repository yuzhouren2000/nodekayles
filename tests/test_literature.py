"""Verification suite: the *Python* reference solver against published data.

This is deliberately independent of the C solver and of the raw TSV files:
it recomputes values from scratch with `src/nodekayles/solver.py`.

Default scope is fast (~1 min).  Use --full for the big variant (slow: the
1 x 1001 and larger boards are better run with scripts/analyze.py, which
checks the same facts against the C-solver data).

Sources of the ground truth: see notes/03_literature.md.

Run:  python3 tests/test_literature.py [--full]
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))

from nodekayles import Solver  # noqa: E402

FAILS = []


def check(label, got, want):
    if got != want:
        FAILS.append(label)
        print(f"  FAIL {label}: got {got}, want {want}")
    else:
        print(f"  ok   {label}")


def load_a002187(limit):
    path = os.path.join(ROOT, "data", "A002187.txt")
    terms = {}
    if not os.path.exists(path):
        return terms
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) == 2 and int(p[0]) <= limit:
                terms[int(p[0])] = int(p[1])
    return terms


def main():
    full = "--full" in sys.argv
    s = Solver()

    npath = 120 if full else 60
    print(f"[1] 1 x n vs OEIS A002187, n = 0..{npath}")
    a = load_a002187(npath)
    bad = [n for n, want in sorted(a.items()) if s.grid_value(1, n) != want]
    check(f"1xn, {len(a)} terms match", bad, [])

    print("[2] 2 x n: G = 1,0,1,0,...  (Brown et al., Theorem 2)")
    N = 30 if full else 12
    check(f"2xn, n=1..{N}", [s.grid_value(2, n) for n in range(1, N + 1)], [1, 0] * (N // 2))

    print("[3] 3 x n vs OEIS A316632, n = 1..17")
    a316632 = [2, 1, 1, 0, 3, 3, 2, 2, 2, 3, 3, 5, 2, 4, 1, 3, 2]
    check("3xn, n=1..17", [s.grid_value(3, n) for n in range(1, 18)], a316632)

    print("[4] 4 x n = 0 (Brown et al., Theorem 6)")
    N = 14 if full else 9
    check(f"4xn, n=1..{N}", [s.grid_value(4, n) for n in range(1, N + 1)], [0] * N)

    print("[5] parity theorem on a small table (odd x odd > 0, even x even = 0)")
    lim = 8 if full else 7
    probs = []
    for m in range(1, lim + 1):
        for n in range(1, lim + 1):
            g = s.grid_value(m, n)
            if m % 2 == 0 and n % 2 == 0 and g != 0:
                probs.append((m, n, g))
            if m % 2 == 1 and n % 2 == 1 and g == 0:
                probs.append((m, n, g))
    check(f"1..{lim} x 1..{lim} parity", probs, [])

    print("[6] centre is a winning first move on odd x odd boards (3x3, 5x5)")
    probs = []
    for k in (3, 5):
        g, moves = s.winning_moves(k, k)
        if (k // 2, k // 2) not in moves:
            probs.append((k, k))
    check("centre wins on 3x3, 5x5", probs, [])

    print()
    if FAILS:
        print(f"{len(FAILS)} FAILURE(S)")
        return 1
    print("ALL CHECKS PASSED")
    print("solver stats:", s.stats())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
