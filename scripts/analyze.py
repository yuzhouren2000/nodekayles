"""Analyse the computed Grundy-value data: strips, periodicity, P-positions,
and verification against the published values.

Reads data/raw/*.tsv (columns: m, n, G; value is symmetric, G(m x n) = G(n x m))
and produces
  data/grid_values.csv   all computed values, one line per unordered pair
  stdout                 per-family sequences, periodicity search, checks

Run:  python3 scripts/analyze.py
"""

import glob
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

# Ground truth from the literature (see notes/03_literature.md)
A002187_EXCEPTIONS = {0, 14, 16, 17, 31, 34, 51}  # OEIS: period 34 except these
A316632 = [2, 1, 1, 0, 3, 3, 2, 2, 2, 3, 3, 5, 2, 4, 1, 3, 2]  # 3 x n, n=1..17
BROWN_TABLE1 = {  # readable cells of Brown et al. 2020 Table 1
    (1, 1): 1, (1, 2): 1, (1, 3): 2, (1, 4): 0, (1, 5): 3, (1, 6): 1,
    (1, 7): 1, (1, 8): 0, (1, 9): 3, (1, 10): 3,
    (2, 1): 1, (2, 2): 0, (2, 3): 1, (2, 4): 0, (2, 5): 1, (2, 6): 0,
    (2, 7): 1, (2, 8): 0, (2, 9): 1, (2, 10): 0,
    (3, 1): 2, (3, 2): 1, (3, 3): 1, (3, 4): 0, (3, 5): 3, (3, 6): 3,
    (3, 7): 2, (3, 8): 2, (3, 9): 2, (3, 10): 3,
    (4, 1): 0, (4, 2): 0, (4, 3): 0, (4, 4): 0, (4, 5): 0, (4, 6): 0,
    (4, 7): 0, (4, 8): 0, (4, 9): 0, (4, 10): 0,
    (5, 1): 3, (5, 2): 1, (5, 3): 3, (5, 4): 0, (5, 5): 3, (5, 6): 3,
    (5, 7): 2,
    (6, 1): 1, (6, 2): 0, (6, 3): 3, (6, 4): 0, (6, 5): 3, (6, 6): 0,
    (6, 8): 0, (6, 10): 0,
    (7, 1): 1, (7, 2): 1, (7, 3): 2, (7, 4): 0, (7, 5): 2,
    (8, 1): 0, (8, 2): 0, (8, 3): 2, (8, 4): 0, (8, 6): 0, (8, 8): 0,
    (8, 10): 0,
    (9, 1): 3, (9, 2): 1, (9, 3): 2, (9, 4): 0,
    (10, 1): 3, (10, 2): 0, (10, 3): 3, (10, 4): 0, (10, 6): 0,
    (10, 8): 0, (10, 10): 0,
}


def load():
    """dict (m, n) -> G, symmetric in m <-> n."""
    g = {}
    for path in sorted(glob.glob(os.path.join(ROOT, "data", "raw", "*.tsv"))):
        with open(path) as f:
            for line in f:
                parts = line.split()
                if len(parts) != 3:
                    continue
                try:
                    m, n, v = map(int, parts)
                except ValueError:
                    continue
                if m < 1 or n < 1:
                    continue  # empty board is trivial (G = 0)
                g[(min(m, n), max(m, n))] = v
    return g


def get(g, m, n):
    return g.get((min(m, n), max(m, n)))


def find_period(seq, min_period=1, max_period=60, min_repeats=3):
    """Smallest p and s with seq[i] = seq[i+p] for all i >= s, needing
    at least min_repeats full periods in the data. Returns (s, p) or None."""
    L = len(seq)
    for p in range(min_period, max_period + 1):
        if L < p * min_repeats:
            continue
        for s in range(0, L - p * min_repeats + 1):
            if all(seq[i] == seq[i + p] for i in range(s, L - p)):
                return (s, p)
    return None


def main():
    g = load()
    ms = sorted({m for m, _ in g})
    print(f"loaded {len(g)} distinct boards, rows m in {ms[0]}..{ms[-1]}\n")

    out = os.path.join(ROOT, "data", "grid_values.csv")
    with open(out, "w") as f:
        f.write("rows,cols,grundy\n")
        for (m, n) in sorted(g):
            f.write(f"{m},{n},{g[(m, n)]}\n")
    print(f"wrote {out}\n")

    # human-readable matrix for small boards
    lim = 13
    lines = []
    width = 3
    header = "      " + "".join(f"{n:>{width}}" for n in range(1, lim + 1))
    lines.append(header)
    for m in range(1, lim + 1):
        row = f"{m:>4}: "
        for n in range(1, lim + 1):
            v = get(g, m, n)
            row += f"{v if v is not None else '?':>{width}}"
        lines.append(row)
    tbl = "\n".join(lines)
    tpath = os.path.join(ROOT, "data", "grid_table.txt")
    with open(tpath, "w") as f:
        f.write("Grundy values G(m x n), blank entries '?' not computed\n\n")
        f.write(tbl + "\n")
    print(f"wrote {tpath}\n{tbl}\n")

    print("=== strip families  G(m x n), n = 1, 2, 3, ... ===")
    for m in ms:
        ns = sorted(n for (a, b) in g for n in ((b,) if a == m else (a,) if b == m else ()))
        seq = [get(g, m, n) for n in ns]
        print(f"m={m:2d} (n={ns[0]}..{ns[-1]}): " + " ".join(str(x) for x in seq))
        consecutive = ns == list(range(ns[0], ns[-1] + 1))
        per = find_period(seq, min_repeats=4) if consecutive else None
        if per:
            s, p = per
            print(f"          observed periodic from n={ns[s]} with period {p}"
                  f" (within computed range, {len(seq)} terms)")
        elif not consecutive:
            print("          period search skipped: missing sizes in this family")
        else:
            print("          no period <= 60 observed within computed range")
    print()

    print("=== P-positions (G = 0) among computed boards ===")
    zeros = sorted((m, n) for (m, n), v in g.items() if v == 0)
    txt = ", ".join(f"{m}x{n}" for m, n in zeros)
    print(txt if txt else "none")
    print(f"({len(zeros)} of {len(g)} computed boards are P-positions)")
    print()

    print("=== verification against published values ===")
    ok = True

    # 1 x n vs OEIS A002187 (period 34 with the listed exceptions)
    bad = []
    a = None
    apath = os.path.join(ROOT, "data", "A002187.txt")
    if os.path.exists(apath):
        a = {}
        with open(apath) as f:
            for line in f:
                p = line.split()
                if len(p) == 2:
                    a[int(p[0])] = int(p[1])
    if a:
        for n, want in a.items():
            if get(g, 1, n) is not None and get(g, 1, n) != want:
                bad.append(n)
        print(f"1 x n versus OEIS A002187 ({len(a)} terms): "
              + ("ALL MATCH" if not bad else f"MISMATCHES at n={bad[:10]}"))
        ok &= not bad

    # 3 x n vs OEIS A316632
    bad = [(n, get(g, 3, n), want) for n, want in enumerate(A316632, 1)
           if get(g, 3, n) is not None and get(g, 3, n) != want]
    print(f"3 x n versus OEIS A316632 (n=1..17): "
          + ("ALL MATCH" if not bad else f"MISMATCHES {bad}"))
    ok &= not bad

    # 2 x n
    bad = [(n, get(g, 2, n)) for n in range(1, 152)
           if get(g, 2, n) is not None and get(g, 2, n) != (1 if n % 2 else 0)]
    print("2 x n versus 1,0,1,0,...: " + ("ALL MATCH" if not bad else f"MISMATCHES {bad}"))
    ok &= not bad

    # 4 x n = 0
    bad = [(n, get(g, 4, n)) for n in range(1, 200)
           if get(g, 4, n) is not None and get(g, 4, n) != 0]
    print("4 x n versus constant 0 (Theorem 6, Brown et al.): "
          + ("ALL MATCH" if not bad else f"MISMATCHES {bad}"))
    ok &= not bad

    # Brown Table 1
    bad = []
    for (m, n), want in sorted(BROWN_TABLE1.items()):
        got = get(g, m, n)
        if got is not None and got != want:
            bad.append((m, n, got, want))
    print(f"Brown et al. Table 1 readable cells ({len(BROWN_TABLE1)}): "
          + ("ALL MATCH" if not bad else f"MISMATCHES {bad}"))
    ok &= not bad

    # structural theorems
    bad = [(m, n, v) for (m, n), v in g.items()
           if m % 2 == 0 and n % 2 == 0 and v != 0]
    print("even x even always G = 0: " + ("HOLDS on data" if not bad else f"VIOLATED {bad}"))
    bad = [(m, n) for (m, n), v in g.items() if m % 2 == 1 and n % 2 == 1 and v == 0]
    print("odd x odd always G > 0:   " + ("HOLDS on data" if not bad else f"VIOLATED {bad}"))
    print()

    print("=== conjecture monitors (see RESULTS.md) ===")
    # Conjecture C: min(m,n) >= 4 and 4 | m or 4 | n  =>  G = 0
    vio, conf = [], []
    for (m, n), v in g.items():
        if min(m, n) >= 4 and (m % 4 == 0 or n % 4 == 0):
            (conf if v == 0 else vio).append((m, n, v))
    print(f"Conjecture C (min>=4 and multiple of 4 => P): {len(conf)} confirmations, "
          + ("no violations" if not vio else f"VIOLATIONS: {vio}"))
    # sporadic zeros: not explained by Conjecture C families or parity/4-side rules
    sporadic = sorted((m, n, v) for (m, n), v in g.items()
                      if v == 0 and min(m, n) >= 3
                      and m % 4 != 0 and n % 4 != 0 and not (m % 2 == 0 and n % 2 == 0))
    print(f"P-positions with both sides >= 3 and no side divisible by 4: "
          f"{sporadic if sporadic else 'none'}")
    # m = 5 family: zeros iff 4 | n
    z5 = sorted(n for n in range(1, 60) if get(g, 5, n) == 0)
    bad5 = [n for n in z5 if n % 4 != 0]
    print(f"5 x n zeros at n = {z5} (all divisible by 4: {not bad5})")
    # 6 x odd
    bad = [(n, v) for n in range(1, 60) if n % 2 == 1 and get(g, 6, n) == 0]
    print("6 x odd nonzero: " + ("HOLDS on data" if not bad else f"VIOLATED {bad}"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
