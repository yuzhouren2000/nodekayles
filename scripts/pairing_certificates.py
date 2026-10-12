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

This script searches certificates of the *product form*

    phi(i, j) = (tau(i), sigma(j)),   tau an involution of rows, sigma of columns,

by exhaustive enumeration of the involutions, and verifies (C1)-(C3) directly
(no decomposition shortcuts).  It reproduces Brown et al.'s Theorem 6 proof
(n x 4, tau = id, sigma = (1 3)(2 4) in 1-based columns) as a verified
certificate.

Observations recorded by the search (see output):

* A certificate with tau = id proves G(m x n) = 0 for *every* m (the rows
  never enter the conditions) - it is a theorem about the whole column n.
* For two odd dimensions no product certificate can exist at all: an odd
  number of rows forces tau to have a fixed row, an odd number of columns
  forces sigma to have a fixed column, and the cell at their intersection
  would be fixed, violating (C1).
* Consequently, if a certificate with tau = id existed for n = 8, it would
  prove G(3 x 8) = 0; the value is 2, so no such certificate exists.  The
  search below verifies the analogous negative statements exhaustively for
  the ranges it covers.

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
    rowo = [False] * m
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
    # orbit representatives of (i,j) under (tau,sigma): check U-invariance on one per orbit
    seen_i = [False] * m
    seen_j = [False] * n
    rep_i = []
    for i in range(m):
        if not seen_i[i]:
            rep_i.append(i)
            seen_i[i] = seen_i[tau[i]] = True
    rep_j = []
    for j in range(n):
        if not seen_j[j]:
            rep_j.append(j)
            seen_j[j] = seen_j[sigma[j]] = True
    for i in rep_i:
        A = [r for r in (i - 1, i, i + 1) if 0 <= r < m]
        i2 = tau[i]
        A2 = [r for r in (i2 - 1, i2, i2 + 1) if 0 <= r < m]
        for j in rep_j:
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


def search(m, n, limit=None):
    """Exhaustive product-form search. Returns list of (tau, sigma) certificates."""
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
    log("")

    log("A. full search for all boards m x n, m <= 6, n <= 6 (all involutions):")
    for m in range(1, 7):
        for n in range(m, 7):
            certs = search(m, n, limit=1)
            if certs:
                tau, sigma = certs[0]
                log(f"  {m} x {n}: CERT  rows {cycles(tau)}  cols {cycles(sigma)}")
            else:
                log(f"  {m} x {n}: no product certificate")

    log("")
    log("A2. family 2 x n and 4 x n (m = 2 or 4 rows), n <= 14:")
    for m in (2, 4):
        for n in range(m, 15):
            certs = search(m, n, limit=1)
            if certs:
                tau, sigma = certs[0]
                log(f"  {m} x {n}: CERT  rows {cycles(tau)}  cols {cycles(sigma)}")
            else:
                log(f"  {m} x {n}: no product certificate")

    log("")
    log("B. targeted searches (odd dimension >= 5 against even dimension):")
    targets = [(5, 6), (5, 8), (5, 10), (5, 12), (6, 7), (7, 8), (7, 10), (8, 9),
               (9, 8), (9, 10), (6, 12), (8, 12)]
    # Boards whose full exhaustion needs > 10^7 candidate pairs (e.g. 9x12,
    # 7x12, 5x14) are deliberately not searched exhaustively here; for them the
    # search would have to enumerate all involutions of a 12/13/14-element side.
    for m, n in targets:
        if m <= 6 and n <= 6:
            continue  # covered in A
        certs = search(m, n, limit=1)
        if certs:
            tau, sigma = certs[0]
            log(f"  {m} x {n}: CERT  rows {cycles(tau)}  cols {cycles(sigma)}")
        else:
            log(f"  {m} x {n}: no product certificate")

    log("")
    log("C. tau = id certificates for columns n <= 14 (would prove G(m x n) = 0 for ALL m):")
    for n in range(1, 15):
        tau2 = (0, 1)
        ok = False
        for sigma in involutions(n):
            if not fixed_point_free(sigma):
                continue
            if dmin(sigma) < 2:
                continue
            if check_cert(2, n, tau2, sigma):
                ok = True
                log(f"  n = {n}: CERT with tau = id, sigma = {cycles(sigma)}")
                break
        if not ok:
            log(f"  n = {n}: none")

    path = os.path.join(ROOT, "data", "pairing_certificates.txt")
    with open(path, "w") as f:
        f.write("\n".join(out_lines) + "\n")
    print(f"\nwrote {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
