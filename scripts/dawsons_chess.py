"""Dawson's chess (octal game 0.137) = Node Kayles on the 1 x n board.

For a path of length L (L playable cells in a row) a move at cell i (0-based)
deletes i-1, i, i+1, leaving independent paths of lengths
    max(0, i-1)  and  max(0, L-i-2).
Hence  G(L) = mex over i of G(left) ^ G(right),  G(0) = 0.
This is exact and O(n^2); it produces the sequence of OEIS A002187.

Writes data/raw/s1_1xn.tsv ("1 n G" for n = 0..N) and checks against
data/A002187.txt when available.
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

NMAX = 1001


def dawson(nmax):
    G = [0] * (nmax + 1)
    for L in range(1, nmax + 1):
        seen = set()
        for i in range(L):
            l = max(0, i - 1)
            r = max(0, L - i - 2)
            seen.add(G[l] ^ G[r])
        g = 0
        while g in seen:
            g += 1
        G[L] = g
    return G


def main():
    G = dawson(NMAX)
    out = os.path.join(ROOT, "data", "raw", "s1_1xn.tsv")
    with open(out, "w") as f:
        for n in range(0, NMAX + 1):
            f.write(f"1\t{n}\t{G[n]}\n")
    print(f"wrote {out}  (n=0..{NMAX})")

    path = os.path.join(ROOT, "data", "A002187.txt")
    if os.path.exists(path):
        bad = 0
        checked = 0
        with open(path) as f:
            for line in f:
                parts = line.split()
                if len(parts) != 2:
                    continue
                n, v = int(parts[0]), int(parts[1])
                if n > NMAX:
                    continue
                checked += 1
                if G[n] != v:
                    bad += 1
                    print(f"  MISMATCH n={n}: ours {G[n]}, A002187 {v}")
        print(f"checked {checked} terms against OEIS A002187: {'ALL MATCH' if bad == 0 else f'{bad} mismatches'}")

    # eventual periodicity check on our own sequence
    for p in range(1, 61):
        s = None
        for start in range(0, len(G) - 2 * p):
            if all(G[i] == G[i + p] for i in range(start, len(G) - p)):
                s = start
                break
        if s is not None and s < 100:
            print(f"period {p} from index {s} (a(n+{p})=a(n) for all n>={s}, n<={NMAX})")
            break
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
