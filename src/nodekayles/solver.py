r"""Exact Sprague-Grundy solver for Node Kayles on rectangular grid boards.

Game (normal play), the one studied here
----------------------------------------
Players alternately put a stone on an empty cell of an m x n board such that
none of the four orthogonally adjacent cells already contains a stone.  A
player with no legal move loses.  This is exactly Node Kayles on the grid
graph P_m [] P_n (Cartesian product of two paths).

Representation
--------------
A position is the set of cells that remain *playable*, i.e. cells that are
neither occupied nor orthogonally adjacent to an occupied cell.  A move at a
playable cell v deletes v and its grid neighbours from the playable set.  The
playable set can fall apart into connected components, and by the
Sprague-Grundy theorem the value of a position is the xor of the values of its
components, so:

    G(T) = mex over v in T of  xor_{C component of T \ N[v]} G(C),   G(empty)=0

Memoisation / canonical form
----------------------------
The value of a connected component depends only on its isomorphism class as an
induced subgraph.  Every connected component met in the recursion is reduced
to a canonical form: shifted to the origin of its bounding box, transposed if
the box is taller than wide, and minimised over the symmetry group of the box
(the Klein four-group of flips/180-degree rotation; the full dihedral group of
order 8 when the box is square).  The Grundy value is then memoised under
(box height, box width, canonical mask); one memo serves every board computed
in the same process, across all shapes and translations.

Values are nimbers (non-negative ints).  value == 0 means: the player to move
loses with optimal play (a P-position).  value > 0 (an N-position) means the
player to move wins.
"""

from __future__ import annotations

import sys

sys.setrecursionlimit(1_000_000)

__all__ = ["Solver", "mex", "grid_value", "outcome"]


def mex(values) -> int:
    """Minimum excluded non-negative integer."""
    g = 0
    while g in values:
        g += 1
    return g


class Solver:
    """Exact Grundy-value solver with a memo shared across all boards."""

    def __init__(self) -> None:
        self.memo: dict[tuple[int, int, int], int] = {}
        self._neigh: dict[tuple[int, int], list[int]] = {}
        self._perms: dict[tuple[int, int], tuple[tuple[int, ...], ...]] = {}
        self._canon_cache: dict[tuple[int, int, int], tuple[int, int, int]] = {}

    # ------------------------------------------------------------------ geometry

    def neighbours(self, h: int, w: int) -> list[int]:
        """Neighbour bitmask of every cell in an h x w box (row-major, stride w)."""
        nb = self._neigh.get((h, w))
        if nb is None:
            nb = [0] * (h * w)
            for r in range(h):
                for c in range(w):
                    m = 0
                    if r:
                        m |= 1 << ((r - 1) * w + c)
                    if r + 1 < h:
                        m |= 1 << ((r + 1) * w + c)
                    if c:
                        m |= 1 << (r * w + c - 1)
                    if c + 1 < w:
                        m |= 1 << (r * w + c + 1)
                    nb[r * w + c] = m
            self._neigh[(h, w)] = nb
        return nb

    def _permutations(self, h: int, w: int) -> tuple[tuple[int, ...], ...]:
        """Graph automorphisms of the h x w box grid (with h <= w), as index maps."""
        key = (h, w)
        perms = self._perms.get(key)
        if perms is None:
            fns = [
                lambda r, c: (r, c),
                lambda r, c: (h - 1 - r, c),
                lambda r, c: (r, w - 1 - c),
                lambda r, c: (h - 1 - r, w - 1 - c),
            ]
            if h == w:
                fns += [
                    lambda r, c: (c, r),
                    lambda r, c: (c, h - 1 - r),
                    lambda r, c: (w - 1 - c, r),
                    lambda r, c: (w - 1 - c, h - 1 - r),
                ]
            perms = []
            for f in fns:
                p = [0] * (h * w)
                for r in range(h):
                    for c in range(w):
                        rr, cc = f(r, c)
                        p[r * w + c] = rr * w + cc
                perms.append(tuple(p))
            perms = tuple(perms)
            self._perms[key] = perms
        return perms

    @staticmethod
    def _apply(perm: tuple[int, ...], mask: int) -> int:
        out = 0
        while mask:
            b = mask & -mask
            out |= 1 << perm[b.bit_length() - 1]
            mask ^= b
        return out

    @staticmethod
    def _transpose(h: int, w: int, mask: int) -> int:
        out = 0
        while mask:
            b = mask & -mask
            r, c = divmod(b.bit_length() - 1, w)
            out |= 1 << (c * h + r)
            mask ^= b
        return out

    @classmethod
    def _to_bbox(cls, h: int, w: int, mask: int) -> tuple[int, int, int]:
        """Shift a mask to the origin of its bounding box."""
        rmin, rmax, cmin, cmax = h, -1, w, -1
        m = mask
        while m:
            b = m & -m
            r, c = divmod(b.bit_length() - 1, w)
            if r < rmin:
                rmin = r
            if r > rmax:
                rmax = r
            if c < cmin:
                cmin = c
            if c > cmax:
                cmax = c
            m ^= b
        bh, bw = rmax - rmin + 1, cmax - cmin + 1
        if (rmin, cmin) == (0, 0) and bh == h and bw == w:
            return bh, bw, mask
        out = 0
        m = mask
        while m:
            b = m & -m
            r, c = divmod(b.bit_length() - 1, w)
            out |= 1 << ((r - rmin) * bw + (c - cmin))
            m ^= b
        return bh, bw, out

    def canonical(self, h: int, w: int, mask: int) -> tuple[int, int, int]:
        """Canonical memo key (box height, box width, mask) of a component."""
        ck = (h, w, mask)
        res = self._canon_cache.get(ck)
        if res is None:
            bh, bw, m = self._to_bbox(h, w, mask)
            if bh > bw:
                m = self._transpose(bh, bw, m)
                bh, bw = bw, bh
            best = m
            for p in self._permutations(bh, bw):
                t = self._apply(p, m)
                if t < best:
                    best = t
            res = (bh, bw, best)
            self._canon_cache[ck] = res
        return res

    def components(self, h: int, w: int, mask: int) -> list[int]:
        """Connected components of an induced subgraph, as bitmasks."""
        nb = self.neighbours(h, w)
        comps = []
        while mask:
            comp = mask & -mask
            frontier = comp
            while frontier:
                v = frontier & -frontier
                frontier ^= v
                add = nb[v.bit_length() - 1] & mask & ~comp
                comp |= add
                frontier |= add
            comps.append(comp)
            mask ^= comp
        return comps

    # ------------------------------------------------------------------ game value

    def value(self, h: int, w: int, mask: int) -> int:
        """Grundy value of the position whose playable cells are `mask` inside an h x w box."""
        total = 0
        for comp in self.components(h, w, mask):
            total ^= self._value_of(self.canonical(h, w, comp))
        return total

    def _value_of(self, key: tuple[int, int, int]) -> int:
        """Grundy value of a connected component given by a canonical key."""
        v = self.memo.get(key)
        if v is not None:
            return v
        h, w, mask = key
        nb = self.neighbours(h, w)
        reachable = set()
        m = mask
        while m:
            b = m & -m
            m ^= b
            rest = mask & ~(b | nb[b.bit_length() - 1])
            s = 0
            for comp in self.components(h, w, rest):
                s ^= self._value_of(self.canonical(h, w, comp))
            reachable.add(s)
        v = mex(reachable)
        self.memo[key] = v
        return v

    # ------------------------------------------------------------------ convenience

    def grid_value(self, rows: int, cols: int) -> int:
        """Grundy value of the full rows x cols board."""
        if rows <= 0 or cols <= 0:
            return 0
        return self.value(rows, cols, (1 << (rows * cols)) - 1)

    def winning_moves(self, rows: int, cols: int) -> tuple[int, list[tuple[int, int]]]:
        """(value, list of first moves that move to a P-position) for the full board."""
        g = self.grid_value(rows, cols)
        nb = self.neighbours(rows, cols)
        full = (1 << (rows * cols)) - 1
        moves = []
        for i in range(rows * cols):
            rest = full & ~((1 << i) | nb[i])
            if self.value(rows, cols, rest) == 0:
                moves.append(divmod(i, cols))
        return g, moves

    def stats(self) -> dict:
        return {
            "memo_entries": len(self.memo),
            "canon_cache": len(self._canon_cache),
        }


def outcome(g: int) -> str:
    """'P' if the player to move loses (P-position), else 'N'."""
    return "P" if g == 0 else "N"


_SOLVER: Solver | None = None


def _solver() -> Solver:
    global _SOLVER
    if _SOLVER is None:
        _SOLVER = Solver()
    return _SOLVER


def grid_value(rows: int, cols: int) -> int:
    """Grundy value of the rows x cols board, using a shared memo across calls."""
    return _solver().grid_value(rows, cols)
