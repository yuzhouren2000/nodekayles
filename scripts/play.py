"""Play Node Kayles on a small grid against an optimal opponent.

The human moves with 'r c' (0-based row and column); the opponent replies with
a move that moves the game to a P-position whenever one exists (perfect play
from then on, computed by the reference solver).

Usage:
    python3 scripts/play.py [m] [n]
eg.
    python3 scripts/play.py 5 5

Boards larger than about 6x8 are not practical for the pure-Python solver.
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))

from nodekayles import Solver  # noqa: E402


def show(m, n, occ, dead):
    print("    " + " ".join(f"{c}" for c in range(n)))
    for r in range(m):
        cells = []
        for c in range(n):
            i = r * n + c
            if occ >> i & 1:
                cells.append("X")
            elif dead >> i & 1:
                cells.append(".")
            else:
                cells.append("-")
        print(f" {r}  " + "  ".join(cells))
    print("   X = stone, . = unplayable, - = playable")


def main():
    m = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    if m * n > 48:
        print("keep m*n <= 48 for the pure-Python solver")
        return 1

    s = Solver()
    nb = s.neighbours(m, n)
    full = (1 << (m * n)) - 1
    occ = 0  # stones
    print(f"Node Kayles on {m} x {n}, you move first (normal play: last move wins).")
    print("Enter moves as: row col   (0-based).  Ctrl-D quits.")

    while True:
        # playable = full minus stones and their neighbours
        playable = full & ~occ
        for k in range(m * n):
            if occ >> k & 1:
                playable &= ~((1 << k) | nb[k])

        if playable == 0:
            print("No legal moves for you — YOU LOSE.")
            show(m, n, occ, 0)
            return 0

        show(m, n, occ, full & ~playable & ~occ)

        # human move
        try:
            line = input("your move: ").strip()
        except EOFError:
            print()
            return 0
        if not line:
            continue
        try:
            r, c = map(int, line.split())
        except ValueError:
            print("usage: row col")
            continue
        i = r * n + c
        if not (0 <= r < m and 0 <= c < n) or not (playable >> i & 1):
            print("not a legal move")
            continue
        occ |= 1 << i

        # opponent: recompute playable; move to a P-position when possible
        playable = full & ~occ
        for k in range(m * n):
            if occ >> k & 1:
                playable &= ~((1 << k) | nb[k])
        if playable == 0:
            print("The opponent has no move — YOU WIN.")
            show(m, n, occ, 0)
            return 0
        best = None
        for k in range(m * n):
            if not (playable >> k & 1):
                continue
            rest = playable & ~((1 << k) | nb[k])
            if s.value(m, n, rest) == 0:
                best = k
                break
        if best is None:  # every move loses in the long run: play anything
            best = (playable & -playable).bit_length() - 1
            print("opponent plays a losing move (you should be winning)")
        else:
            print(f"opponent plays ({best // n}, {best % n})")
        occ |= 1 << best


if __name__ == "__main__":
    raise SystemExit(main())
