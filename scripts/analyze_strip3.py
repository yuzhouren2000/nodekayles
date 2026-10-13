"""Verify and analyse a full/end-gap/both-gap three-row first-move survey.

Run: python3 scripts/analyze_strip3.py [data/strip3_survey.jsonl]
The middle-row formula is a theorem; finite prefix matches are observations.
"""

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from analyze import A316632, find_period


def mex(values):
    g = 0
    while g in values:
        g += 1
    return g


def screen_middle(n, a, forbidden):
    """Check each reflected pair of middle moves with explicit premises.

    An unknown a(q) differs from a(p) if a known witness has the same value
    as a(p), but is separated from q by a zero-tail forbidden length gap.
    A screen never assigns an SG value to an unknown auxiliary position.
    """
    exclusions, winning, unresolved = [], [], []
    for p in range((n + 1) // 2):
        q = n - 1 - p
        if q < len(a):
            if a[p] == a[q]:
                winning.append(p + 1)
            else:
                exclusions.append({"left": p, "right": q, "reason": "known_unequal",
                                   "values": [a[p], a[q]]})
        elif q - p in forbidden:
            exclusions.append({"left": p, "right": q, "reason": "zero_tail_gap",
                               "gap": q - p, "zero_tail_length": q - p - 1})
        elif p < len(a):
            witnesses = [w for w in range(len(a)) if a[w] == a[p] and abs(q - w) in forbidden]
            if witnesses:
                w = min(witnesses)
                exclusions.append({"left": p, "right": q, "reason": "known_equal_witness",
                                   "witness": w, "equal_nimber": a[p], "gap": abs(q - w),
                                   "zero_tail_length": abs(q - w) - 1})
            else:
                unresolved.append([p, q])
        else:
            unresolved.append([p, q])
    status = "middle_win" if winning else "no_middle_win" if not unresolved else "incomplete"
    return {"n": n, "status": status, "winning_columns_in_left_half": winning,
            "exclusions": exclusions, "unresolved_split_lengths": unresolved}


def analyse(path):
    records = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    if not records or [r["n"] for r in records] != list(range(1, len(records) + 1)):
        raise ValueError("survey must be a complete consecutive prefix starting at n=1")
    full = [0] + [r["grundy"] for r in records]
    a = [0] + [r["end_gap"] for r in records]
    d = [0] + [r["both_gaps"] for r in records]
    for seq in (full, a, d):
        if any(type(g) is not int or g < 0 for g in seq):
            raise ValueError("invalid nimber")
    if full[1:len(A316632) + 1] != A316632[:len(records)]:
        raise ValueError("survey disagrees with OEIS A316632")

    known = {}
    for reference in sorted((ROOT / "data/raw").glob("*.tsv")):
        for line in reference.read_text().splitlines():
            parts = line.split()
            if len(parts) != 3:
                continue
            try:
                m, n, g = map(int, parts)
            except ValueError:
                continue
            if n == 3:
                m, n = n, m
            if m == 3:
                if n in known and known[n] != g:
                    raise ValueError(f"conflicting original data for 3x{n}")
                known[n] = g
    for n, g in known.items():
        if n < len(full) and full[n] != g:
            raise ValueError(f"survey disagrees with original C data for 3x{n}")

    rows = []
    for record in records:
        n, options = record["n"], record["options"]
        if len(options) != 3 * n or any(type(g) is not int or g < 0 for g in options):
            raise ValueError(f"invalid move list at n={n}")
        if mex(set(options)) != full[n]:
            raise ValueError(f"mex mismatch at n={n}")
        for j in range(1, n + 1):
            if options[3 * (j - 1) + 1] != (a[j - 1] ^ a[n - j]):
                raise ValueError(f"middle split mismatch at n={n}, j={j}")
            if options[3 * (j - 1)] != options[3 * (j - 1) + 2]:
                raise ValueError(f"vertical reflection mismatch at n={n}, j={j}")
        if options != [g for j in range(n - 1, -1, -1) for g in options[3*j:3*j+3]]:
            raise ValueError(f"horizontal reflection mismatch at n={n}")

        # These are legal options, so mex excludes their values. The end-gap
        # families are closed under middle moves, not under every legal move.
        for j in range(1, n):
            if a[n] == (a[j - 1] ^ d[n - j]):
                raise ValueError(f"A option-exclusion theorem violated at n={n}, j={j}")
        for j in range(2, n):
            if d[n] == (d[j - 1] ^ d[n - j]):
                raise ValueError(f"D option-exclusion theorem violated at n={n}, j={j}")
        middle = [j for j in range(1, n + 1) if options[3 * (j - 1) + 1] == 0]
        outer = [j for j in range(1, n + 1) if options[3 * (j - 1)] == 0]
        if n % 2 and (n + 1) // 2 not in middle:
            raise ValueError("odd full strip lacks its central winning move")
        if n >= 3 and n % 2 and d[n] == 0:
            raise ValueError("odd D strip lacks its central winning move")
        rows.append({"n": n, "grundy": full[n], "end_gap": a[n], "both_gaps": d[n],
                     "middle_winning_columns": middle, "outer_winning_columns": outer,
                     "max_option": max(options)})

    forbidden = [t + 1 for t in range(1, len(d)) if d[t] == 0]
    for gap in forbidden:
        if any(a[k] == a[k + gap] for k in range(len(a) - gap)):
            raise ValueError(f"zero-tail forbidden-gap theorem violated for gap={gap}")
    certificates = []
    for n in range(2, 2 * len(a), 2):
        matches = [(k, n - 1 - k) for k in range(len(a))
                   if k <= n - 1 - k < len(a) and a[k] == a[n - 1 - k]]
        if matches:
            k, ell = min(matches, key=lambda pair: max(pair))
            certificates.append({"n": n, "row": 2, "columns": [k + 1, ell + 1],
                                 "split_lengths": [k, ell], "equal_nimber": a[k],
                                 "beyond_full_survey": n > len(records)})
    screens = [screen_middle(n, a, forbidden) for n in range(2, 2 * len(a) + 1, 2)]
    impossible_a_periods = sorted({p for gap in forbidden for p in range(1, gap + 1) if gap % p == 0})
    for screen in screens:
        n = screen["n"]
        if n <= len(records):
            known_columns = [j for j in rows[n - 1]["middle_winning_columns"] if j <= n // 2]
            if screen["winning_columns_in_left_half"] != known_columns or screen["status"] == "incomplete":
                raise ValueError("middle-screen proof disagrees with a fully surveyed board")
    return records, rows, full, a, d, forbidden, certificates, screens, impossible_a_periods


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("survey", type=Path, nargs="?", default=ROOT / "data/strip3_survey.jsonl")
    parser.add_argument("--json", action="store_true", help="emit checked analysis as JSON")
    args = parser.parse_args()
    records, rows, full, a, d, forbidden, certificates, screens, impossible_a_periods = analyse(args.survey)
    if args.json:
        print(json.dumps({"n_max": len(records), "rows": rows,
                          "zero_tail_forbidden_gaps": forbidden,
                          "middle_strategy_certificates": certificates,
                          "middle_screens": screens,
                          "impossible_eventual_periods_of_A": impossible_a_periods}, ensure_ascii=False, indent=2))
        return 0
    print(f"Verified {len(records)} consecutive sizes, all first-move mex/reflections, published/original values and middle-family identities.")
    print("F(1..N):", " ".join(map(str, full[1:])))
    print("A(0..N):", " ".join(map(str, a)))
    print("D(0..N):", " ".join(map(str, d)))
    print("Full-strip P sizes:", [r["n"] for r in rows if r["grundy"] == 0])
    print("Even sizes with no middle winning move:", [r["n"] for r in rows if r["n"] % 2 == 0 and not r["middle_winning_columns"]])
    print("Among those, outer-only wins:", [r["n"] for r in rows if not r["middle_winning_columns"] and r["outer_winning_columns"]])
    print("A unequal at these length gaps for every k (using computed zero D tails):", forbidden)
    print("  gaps 2 and 3 have direct structural proofs, independent of the survey.")
    print("Eventual periods of A excluded by the zero-tail theorem (no claim about full-strip periods):", impossible_a_periods)
    print("Constructive wins beyond the full-board survey (no full-board SG claimed):")
    for certificate in certificates:
        if certificate["beyond_full_survey"]:
            print(f"  3x{certificate['n']}: middle columns {certificate['columns']}, "
                  f"A lengths {certificate['split_lengths']}, equal nimber {certificate['equal_nimber']}")
    print("Beyond the survey, every middle move excluded by zero-tail/known-value proofs:",
          [s["n"] for s in screens if s["n"] > len(records) and s["status"] == "no_middle_win"])
    print("Even-size move classification (columns start at 1; outer row may be top or bottom):")
    for row in rows:
        if row["n"] % 2 == 0:
            print(f"  n={row['n']:2} G={row['grundy']} middle={row['middle_winning_columns']} outer={row['outer_winning_columns']}")
    for label, sequence in (("F", full[1:]), ("A", a), ("D", d)):
        period = find_period(sequence, min_repeats=4)
        if period:
            start, length = period
            print(f"{label}: observed suffix repetition at index {start}, period {length}; finite data only.")
        else:
            print(f"{label}: no suffix with four complete repeats of period <=60; eventual periodicity remains open.")
    print("Maximum sampled empty-board/one-end/both-end/first-option nimbers:",
          max(full), max(a), max(d), max(r["max_option"] for r in rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
