"""Winning first moves for small boards.

For every board m x n (m,n <= MMAX) compute the Grundy value and the list of
first moves that move to a P-position (G = 0); for a winning board these are
exactly the first moves that win with optimal play from then on.

Writes data/winning_moves.csv and prints a compact table.

Run:  python3 scripts/winning_moves.py [MMAX]
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))

from nodekayles import Solver  # noqa: E402


def main():
    mmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    s = Solver()
    rows = []
    print("board  G  #winning  winning first moves (row,col) 0-based")
    for m in range(1, mmax + 1):
        for n in range(m, mmax + 1):
            g, moves = s.winning_moves(m, n)
            rows.append((m, n, g, len(moves), moves))
            mv = " ".join(f"({r},{c})" for r, c in moves)
            if len(mv) > 70:
                mv = mv[:70] + " ..."
            print(f"{m}x{n:<2} {g:2d}  {len(moves):3d}      {mv or '-'}")
    out = os.path.join(ROOT, "data", "winning_moves.csv")
    with open(out, "w") as f:
        f.write("rows,cols,grundy,num_winning_moves,moves\n")
        for m, n, g, k, moves in rows:
            f.write(f"{m},{n},{g},{k},\"" + ";".join(f"{r}:{c}" for r, c in moves) + "\"\n")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
