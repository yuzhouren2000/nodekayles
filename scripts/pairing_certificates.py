"""Machine-checkable pairing certificates: prove P-positions without search.

A *pairing certificate* for an m x n board is a map phi on the cells with

  (C1) phi is an involution without fixed points,
  (C2) phi(x) is never adjacent to x,
  (C3) for every cell x, the set U(x) = N[x] union N[phi(x)] is phi-invariant,

where N[x] is the closed neighbourhood.  These three conditions make the
mirror strategy "answer a move at x by a move at phi(x)" a valid second-player
strategy: after every move of the first player, the responder's partner cell is
still playable and the invariant is restored.  Hence a certificate is a
self-contained proof that the board is a P-position (Grundy value 0).

This script constructs certificates of the *product form*

    phi(i, j) = (tau(i), sigma(j)),   tau an involution of rows, sigma of columns,

using the classification proved in notes/06_structural_pairing.md: one exists
iff both dimensions are even or one dimension equals 4.  The exhaustive
search remains available with search(..., exhaustive=True) for independent
small-board checks.  Verification checks (C1)-(C3) directly.  This reproduces
Brown et al.'s Theorem 6 proof
(n x 4, tau = id, sigma = (1 3)(2 4) in 1-based columns) as a verified
certificate.

Consequences of the classification (see output):

* A certificate with tau = id proves G(m x n) = 0 for *every* m (the rows
  never enter the conditions) - it is a theorem about the whole column n.
* For two odd dimensions no product certificate can exist at all: an odd
  number of rows forces tau to have a fixed row, an odd number of columns
  forces sigma to have a fixed column, and the cell at their intersection
  would be fixed, violating (C1).
* Consequently, if a certificate with tau = id existed for n = 8, it would
  prove G(3 x 8) = 0; the value is 2, so no such certificate exists.  The
  path lemma now proves the stronger statement: tau = id permits only n = 4.

Run:  python3 scripts/pairing_certificates.py
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def involutions(n):
    """All involutions of {0..n-1} as tuples p with p(p(i)) = i."""
    result = []

    def rec(i, cur):
        if i == n:
            result.append(tuple(cur))
            return
        if cur[i] != -1:
            rec(i + 1, cur)
            return
        cur[i] = i
        rec(i + 1, cur)
        cur[i] = -1
        for k in range(i + 1, n):
            if cur[k] == -1:
                cur[i], cur[k] = k, i
                rec(i + 1, cur)
                cur[i], cur[k] = -1, -1

    rec(0, [-1] * n)
    return result


def fixed_point_free(p):
    return all(p[i] != i for i in range(len(p)))


def dmin(p):
    return min(abs(p[i] - i) for i in range(len(p)))


def check_cert(m, n, tau, sigma):
    """Direct verification of (C1)-(C3) for phi = tau (x) sigma."""
    if any(tau[tau[i]] != i for i in range(m)):
        return False
    if any(sigma[sigma[j]] != j for j in range(n)):
        return False
    fixed_cell = False
    for i in range(m):
        for j in range(n):
            if tau[i] == i and sigma[j] == j:
                fixed_cell = True
    if fixed_cell:
        return False
    # adjacency: |tau(i)-i| + |sigma(j)-j| >= 2 for all i,j (<=> dmin(tau)+dmin(sigma)>=2)
    if dmin(tau) + dmin(sigma) < 2:
        return False
    # Check every cell.  A row orbit times a column orbit can contain TWO
    # cell orbits, so taking just one Cartesian-product representative is
    # not in itself a justification for checking both.
    for i in range(m):
        A = [r for r in (i - 1, i, i + 1) if 0 <= r < m]
        i2 = tau[i]
        A2 = [r for r in (i2 - 1, i2, i2 + 1) if 0 <= r < m]
        for j in range(n):
            B = [c for c in (j - 1, j, j + 1) if 0 <= c < n]
            j2 = sigma[j]
            B2 = [c for c in (j2 - 1, j2, j2 + 1) if 0 <= c < n]
            U = set()
            for y in B:
                U.add((i, y))
            for x in A:
                U.add((x, j))
            for y in B2:
                U.add((i2, y))
            for x in A2:
                U.add((x, j2))
            V = {(tau[x], sigma[y]) for (x, y) in U}
            if V != U:
                return False
    return True


def cycles(p):
    seen = [False] * len(p)
    out = []
    for i in range(len(p)):
        if not seen[i] and p[i] != i:
            cyc = []
            j = i
            while not seen[j]:
                seen[j] = True
                cyc.append(j + 1)  # 1-based for readability
                j = p[j]
            out.append("(" + " ".join(map(str, cyc)) + ")")
    return " ".join(out) or "id"


def construct_product_certificate(m, n):
    """Return one certificate, or None if none exists (for positive sizes).

    This classifies product certificates, not winning/losing positions.
    The construction costs O(m+n); the proof is in notes/06_structural_pairing.md.
    """
    if m <= 0 or n <= 0:
        raise ValueError("board dimensions must be positive")
    if m == 4:
        return (2, 3, 0, 1), tuple(range(n))
    if n == 4:
        return tuple(range(m)), (2, 3, 0, 1)
    if m % 2 == 0 and n % 2 == 0:
        return tuple(reversed(range(m))), tuple(reversed(range(n)))
    return None


def search(m, n, limit=None, *, exhaustive=False):
    """Return up to limit product certificates (all if limit is None).

    By default the theorem answers nonexistence and limit=1 immediately.
    exhaustive=True bypasses the theorem and enumerates every candidate until
    the requested limit is reached, for independent verification.
    """
    if m <= 0 or n <= 0:
        raise ValueError("board dimensions must be positive")
    if not exhaustive:
        cert = construct_product_certificate(m, n)
        if cert is None:
            return []
        if limit == 1:
            return [cert]
    found = []
    info_m = [(p, fixed_point_free(p), dmin(p)) for p in involutions(m)]
    info_n = [(p, fixed_point_free(p), dmin(p)) for p in involutions(n)]
    # try fixed-point-free parts first - certificates tend to have them
    info_m.sort(key=lambda t: (not t[1], -t[2]))
    info_n.sort(key=lambda t: (not t[1], -t[2]))
    for tau, ft, dt in info_m:
        for sigma, fs, ds in info_n:
            if not (ft or fs):
                continue  # a fixed row and a fixed column would fix a cell
            if dt + ds < 2:
                continue  # some pair of cells at distance 1
            if check_cert(m, n, tau, sigma):
                found.append((tau, sigma))
                if limit and len(found) >= limit:
                    return found
    return found


def main():
    out_lines = []

    def log(s):
        print(s, flush=True)
        out_lines.append(s)

    log("pairing certificates, product form phi(i,j) = (tau(i), sigma(j))")
    log("classification: both dimensions even OR one dimension equals 4")
    log("proof: notes/06_structural_pairing.md; each constructed certificate is checked directly")
    log("")

    def classified(m, n):
        certs = search(m, n, limit=1)
        for tau, sigma in certs:
            if not check_cert(m, n, tau, sigma):
                raise AssertionError(f"invalid constructed certificate for {m}x{n}")
        return certs

    log("A. classification for all boards m x n, m <= 6, n <= 6:")
    for m in range(1, 7):
        for n in range(m, 7):
            certs = classified(m, n)
            if certs:
                tau, sigma = certs[0]
                log(f"  {m} x {n}: CERT  rows {cycles(tau)}  cols {cycles(sigma)}")
            else:
                log(f"  {m} x {n}: no product certificate")

    log("")
    log("A2. family 2 x n and 4 x n (m = 2 or 4 rows), n <= 14:")
    for m in (2, 4):
        for n in range(m, 15):
            certs = classified(m, n)
            if certs:
                tau, sigma = certs[0]
                log(f"  {m} x {n}: CERT  rows {cycles(tau)}  cols {cycles(sigma)}")
            else:
                log(f"  {m} x {n}: no product certificate")

    log("")
    log("B. selected boards (classified without enumerating involutions):")
    targets = [(5, 6), (5, 8), (5, 10), (5, 12), (6, 7), (7, 8), (7, 10), (8, 9),
               (9, 8), (9, 10), (6, 12), (8, 12)]
    targets += [(9, 12), (7, 12), (5, 14)]
    for m, n in targets:
        if m <= 6 and n <= 6:
            continue  # covered in A
        certs = classified(m, n)
        if certs:
            tau, sigma = certs[0]
            log(f"  {m} x {n}: CERT  rows {cycles(tau)}  cols {cycles(sigma)}")
        else:
            log(f"  {m} x {n}: no product certificate")

    log("")
    log("C. tau = id certificates: the path theorem permits only n = 4 (for every m):")
    for n in range(1, 15):
        cert = construct_product_certificate(1, n)
        if cert is not None:
            _, sigma = cert
            if not check_cert(2, n, (0, 1), sigma):
                raise AssertionError(f"invalid column certificate for n={n}")
            log(f"  n = {n}: CERT with tau = id, sigma = {cycles(sigma)}")
        else:
            log(f"  n = {n}: none")

    path = os.path.join(ROOT, "data", "pairing_certificates.txt")
    with open(path, "w") as f:
        f.write("\n".join(out_lines) + "\n")
    print(f"\nwrote {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
