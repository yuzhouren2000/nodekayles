"""Deliberately naive independent solver: plain memoised minimax.

No connected-component decomposition, no canonicalisation, no symmetry
reduction - just the raw definition of the game on a bitmask of playable
cells.  Exponential state space; only usable for small boards.  Its purpose
is to triple-check the optimised solvers on selected boards.

Run:  python3 tests/dumb_solver.py 6 7 5 7 7 7 3 9 5 8
"""

import sys
from functools import lru_cache


def grundy(h, w, playable):
    """Grundy value of the game on 'playable' cells of an h x w board."""

    def neighbours(i):
        r, c = divmod(i, w)
        out = 0
        if r:
            out |= 1 << (i - w)
        if r + 1 < h:
            out |= 1 << (i + w)
        if c:
            out |= 1 << (i - 1)
        if c + 1 < w:
            out |= 1 << (i + 1)
        return out

    nb = [neighbours(i) for i in range(h * w)]

    @lru_cache(maxsize=None)
    def g(mask):
        seen = set()
        m = mask
        while m:
            b = m & -m
            m ^= b
            i = b.bit_length() - 1
            seen.add(g(mask & ~(b | nb[i])))
        v = 0
        while v in seen:
            v += 1
        return v

    return g(playable)


def main():
    args = [int(x) for x in sys.argv[1:]]
    assert len(args) % 2 == 0 and args, "give pairs m n"
    for i in range(0, len(args), 2):
        m, n = args[i], args[i + 1]
        full = (1 << (m * n)) - 1
        v = grundy(m, n, full)
        print(f"{m}x{n}: naive solver Grundy = {v}", flush=True)
        # a P-position has no winning move; a move is winning iff it leaves 0
        # (recompute for the board's files only when cheap)


if __name__ == "__main__":
    main()
