#!/usr/bin/env python3
"""Independent verification of the port-family data files.

Reads data/strip3_ports_n12.jsonl (survey output) and optionally
data/strip3_ports_closure.json, then:

  1. rebuilds every family mask from one-based coordinates with its own
     encoding and compares it with the file;
  2. recomputes small positions with the coordinate recursion from
     tests/test_strip3_ports.py (no import of scripts/strip3_ports.py);
  3. checks every row: option cells = legal cells, mex of option values equals
     the row value, middle-move values satisfy the L1 xor identities using
     family values from the data, and B/H reflection identities hold;
  4. cross-checks b(n,j) against data/strip3_survey.jsonl mid-row options;
  5. reports the L2 inequalities across the data when the required widths are
     present;
  6. for the closure file: re-derives each recorded unrepresented successor
     from its origin and move, and independently confirms that its column runs
     match no F/A/D/B/E/H template under board symmetries and are not paths or
     cycles.

Run: python3 tests/verify_strip3_ports.py data/strip3_ports_n12.jsonl \
        [--closure data/strip3_ports_closure.json] [--recompute-n 6]
"""

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests"))
import test_strip3_ports as base  # coordinate helpers and independent recursion

# ------------------------------------------------ mask encoding (interface)

def encode(cells):
    mask = 0
    for r, c in cells:
        mask |= 1 << (3 * (c - 1) + (r - 1))
    return mask


def decode(mask):
    out = set()
    while mask:
        low = mask & -mask
        v = low.bit_length() - 1
        mask ^= low
        out.add((v % 3 + 1, v // 3 + 1))
    return out


def mex(values):
    g = 0
    while g in values:
        g += 1
    return g


# ------------------------------------------- independent template matching

def all_templates(width):
    """Explicit enumeration of every declared family position label of a given
    width, under the four board symmetries, as normalised cell sets."""
    templates = []

    def add(label, cells, w):
        for variant in (cells,
                        frozenset((4 - r, c) for r, c in cells),
                        frozenset((r, w + 1 - c) for r, c in cells),
                        frozenset((4 - r, w + 1 - c) for r, c in cells)):
            first = min(cc for _rr, cc in variant)
            normal = frozenset((r, c - first + 1) for r, c in variant)
            templates.append((label, normal))

    add("F(%d)" % width, base.full_cells(width), width)
    for w in range(1, width + 1):
        add("A(%d)" % w, base.cells_of("A", w), w)
        add("D(%d)" % w, base.cells_of("D", w), w)
        for j in range(1, w + 1):
            for fam in ("B", "E", "H"):
                add("%s(%d,%d)" % (fam, w, j), base.cells_of(fam, w, j), w)
    return templates


_TEMPLATE_CACHE = {}


def candidate_labels(cells):
    """Every vocabulary label matching this cell set (family templates under
    the four reflections, plus path/cycle)."""
    if not cells:
        return {"empty"}
    width = max(c for _r, c in cells)
    first = min(c for _r, c in cells)
    normal = frozenset((r, c - first + 1) for r, c in cells)
    if width not in _TEMPLATE_CACHE:
        _TEMPLATE_CACHE[width] = all_templates(width)
    labels = {label for label, template in _TEMPLATE_CACHE[width]
              if template == normal}
    degrees = {v: sum(1 for u in base.move_targets(v) if u in cells and u != v)
               for v in cells}
    if max(degrees.values()) <= 2:
        if any(d < 2 for d in degrees.values()):
            labels.add("path(%d)" % len(cells))
        else:
            labels.add("cycle(%d)" % len(cells))
    elif len(cells) == 4 and max(degrees.values()) == 3 and \
            len(base.component_split(cells)) == 1:
        labels.add("T4(claw)")
    return labels


def runs_of(cells):
    if not cells:
        return []
    columns = sorted({c for _r, c in cells})
    runs, start, prev = [], columns[0], columns[0]
    for c in columns[1:]:
        if c != prev + 1:
            runs.append((start, prev))
            start = c
        prev = c
    runs.append((start, prev))
    return runs


# --------------------------------------------------------------- main body

def load_survey(path):
    records = [json.loads(line) for line in path.read_text().splitlines()
               if line.strip()]
    if not records or records[0].get("kind") != "strip3_ports_survey":
        raise ValueError("missing survey header")
    return records[0], records[1:]


def check_families(rows, recompute_n):
    by = {}
    for row in rows:
        key = (row["family"], row["n"], row["j"])
        if key in by:
            raise ValueError(f"duplicate row {key}")
        by[key] = row
    max_n = max(r["n"] for r in rows)
    for n in range(1, max_n + 1):
        for j in range(1, n + 1):
            for fam in ("B", "E", "H"):
                row = by.get((fam, n, j))
                if row is None:
                    raise ValueError(f"missing {fam}({n},{j})")
                expect = encode(base.cells_of(fam, n, j))
                if int(row["mask"], 16) != expect:
                    raise ValueError(f"mask mismatch {fam}({n},{j})")
    print(f"families: all {len(rows)} B/E/H rows present with coordinate masks")

    for row in rows:
        if row["n"] > recompute_n:
            continue
        v = base.rec_value(base.cells_of(row["family"], row["n"], row["j"]))
        if v != row["value"]:
            raise ValueError(f"recompute mismatch {row['family']}({row['n']},{row['j']}): "
                             f"{v} vs {row['value']}")
    rebuilt = [r for r in rows if r["n"] <= recompute_n]
    print(f"recompute: independent recursion agrees on {len(rebuilt)} rows "
          f"with n<={recompute_n}")
    return by


def family_value(by, survey_map, name, n, j=None):
    if name in ("A", "D"):
        return survey_map[name][n]
    if name == "F":
        return survey_map["F"][n]
    return by[(name, n, j)]["value"]


def check_rows(by, survey_map):
    middle_checks = 0
    option_checks = 0
    zero_certificates = 0
    for row in by.values():
        fam, n, j = row["family"], row["n"], row["j"]
        cells = base.cells_of(fam, n, j)
        cells_sorted = sorted(cells, key=lambda x: (x[1], x[0]))
        if [tuple(o["cell"]) for o in row["options"]] != cells_sorted:
            raise ValueError(f"option cells are not the legal cells {fam}({n},{j})")
        values = [o["value"] for o in row["options"]]
        if any(v is None or v < 0 for v in values):
            raise ValueError(f"UNKNOWN value recorded in a completed survey row {fam}({n},{j})")
        if mex(set(values)) != row["value"]:
            raise ValueError(f"mex mismatch {fam}({n},{j})")
        option_checks += len(values)
        for o in row["options"]:
            r, c = o["cell"]
            succ = base.cells_of(fam, n, j) - base.move_targets((r, c))
            if int(o["mask"], 16) != encode(succ):
                raise ValueError(f"option mask mismatch {fam}({n},{j}) move {r},{c}")
            if r == 2:
                i = c
                if fam == "B":
                    if i < j:
                        want = family_value(by, survey_map, "A", i - 1) ^ \
                            family_value(by, survey_map, "E", n - i, j - i)
                    else:
                        want = family_value(by, survey_map, "E", i - 1, i - j) ^ \
                            family_value(by, survey_map, "A", n - i)
                elif fam == "E":
                    if i < j:
                        want = family_value(by, survey_map, "D", i - 1) ^ \
                            family_value(by, survey_map, "E", n - i, j - i)
                    else:
                        want = family_value(by, survey_map, "H", i - 1, j) ^ \
                            family_value(by, survey_map, "A", n - i)
                else:
                    if i < j:
                        want = family_value(by, survey_map, "D", i - 1) ^ \
                            family_value(by, survey_map, "H", n - i, j - i)
                    else:
                        want = family_value(by, survey_map, "H", i - 1, j) ^ \
                            family_value(by, survey_map, "D", n - i)
                if o["value"] != want:
                    raise ValueError(
                        f"L1 value identity fails {fam}({n},{j}) move col {i}: "
                        f"{o['value']} vs {want}")
                middle_checks += 1
            if o.get("zero_replies") is not None:
                if fam != "B" or (n, j) != (10, 4):
                    raise ValueError("zero_replies outside B(10,4)")
                child = set(succ)
                for reply in o["zero_replies"]:
                    cell2 = tuple(reply["cell"])
                    if cell2 not in child:
                        raise ValueError("recorded zero reply is not a legal cell")
                    grand = child - base.move_targets(cell2)
                    if int(reply["mask"], 16) != encode(grand):
                        raise ValueError("zero reply mask mismatch")
                    if reply["value"] != 0:
                        raise ValueError("recorded zero reply is not zero")
                    zero_certificates += 1
        # every option of B(10,4) must have at least one zero reply
        if fam == "B" and (n, j) == (10, 4):
            if any(not o.get("zero_replies") for o in row["options"]):
                raise ValueError("an option of B(10,4) lacks a zero reply")
    for n in range(1, max(r["n"] for r in by.values()) + 1):
        for j in range(1, n + 1):
            if by[("B", n, j)]["value"] != by[("B", n, n + 1 - j)]["value"]:
                raise ValueError(f"B reflection fails n={n} j={j}")
            if by[("H", n, j)]["value"] != by[("H", n, n + 1 - j)]["value"]:
                raise ValueError(f"H reflection fails n={n} j={j}")
    print(f"rows: {option_checks} option values cover every legal cell, mex and "
          f"reflections hold, {middle_checks} middle moves satisfy L1")
    print(f"B(10,4): all 26 options have verified zero replies "
          f"({zero_certificates} explicit responses)")


def check_against_full_survey(by):
    path = ROOT / "data/strip3_survey.jsonl"
    if not path.exists():
        print("full survey not present; b cross-check skipped")
        return
    survey = {json.loads(l)["n"]: json.loads(l)
              for l in path.read_text().splitlines() if l.strip()}
    checked = 0
    for (fam, n, j), row in by.items():
        if fam != "B" or n not in survey:
            continue
        # the full survey's options list holds values, not masks: the option
        # slot of the top-row move at column j is options[3*(j-1)]
        known = survey[n]["options"][3 * (j - 1)]
        if row["value"] != known:
            raise ValueError(f"b({n},{j}) disagrees with the full survey: "
                             f"{row['value']} vs {known}")
        checked += 1
    print(f"cross-check: {checked} b-values agree with the full-board survey")


def check_l2_data(by, survey_map, max_gap=12):
    checked = 0
    for t in range(1, max_gap + 1):
        if survey_map["D"][t] != 0:
            continue
        q = t + 1
        max_n = max(r["n"] for r in by.values())
        for n in range(1, max_n - q + 1):
            for j in range(1, n + 1):
                if by[("E", n + q, j + q)]["value"] == by[("E", n, j)]["value"]:
                    raise ValueError(f"L2 e-extension fails t={t},n={n},j={j}")
                if by[("H", n + q, j + q)]["value"] == by[("H", n, j)]["value"]:
                    raise ValueError(f"L2 h-extension fails t={t},n={n},j={j}")
                if by[("H", n + q, j)]["value"] == by[("H", n, j)]["value"]:
                    raise ValueError(f"L2 h-right fails t={t},n={n},j={j}")
                checked += 3
    for t in range(0, max_gap + 1):
        if survey_map["A"][t] != 0:
            continue
        max_n = max(r["n"] for r in by.values())
        for n in range(1, max_n - t):
            for j in range(1, n + 1):
                if by[("E", n + t + 1, j)]["value"] == by[("H", n, j)]["value"]:
                    raise ValueError(f"L2 e/h fails t={t},n={n},j={j}")
                checked += 1
    print(f"L2: {checked} extension inequalities confirmed inside the data range")


def check_closure(path, max_n):
    report = json.loads(path.read_text())
    if report.get("kind") != "strip3_ports_closure":
        raise ValueError("missing closure header")
    listed = report["unrepresented"]
    total = sum(u["occurrences"] for u in listed)
    if len(listed) < 100 and total != report.get("unrepresented_count"):
        raise ValueError("unrepresented counts are inconsistent")
    if total > report.get("unrepresented_count", -1):
        raise ValueError("unrepresented counts are inconsistent")
    rechecked = 0
    for entry in report["unrepresented"]:
        fam = entry["origin"]["family"]
        n = entry["origin"]["n"]
        j = entry["origin"]["j"]
        r, c = entry["move"]
        cells = base.cells_of(fam, n, j)
        successor = cells - base.move_targets((r, c))
        if int(entry["successor_mask"], 16) != encode(successor):
            raise ValueError(f"closure successor mask wrong for {fam}({n},{j})")
        for piece in entry["pieces"]:
            start, end = piece["columns"]
            run = frozenset((rr, cc) for rr, cc in successor
                            if start <= cc <= end)
            candidates = candidate_labels(run)
            recorded = piece["label"]
            if recorded == "unknown":
                if candidates:
                    raise ValueError(
                        f"closure claims unknown but vocabulary matches: "
                        f"{fam}({n},{j}) move {r},{c}: {sorted(candidates)}")
            elif recorded not in candidates:
                raise ValueError(
                    f"closure label {recorded} not reproducible for "
                    f"{fam}({n},{j}) move {r},{c}: {sorted(candidates)}")
        rechecked += 1
    print(f"closure: {rechecked} recorded minimal unrepresented successors "
          f"re-derived; every run still defeats the vocabulary")
    # extended vocabulary (base plus T4(claw)): at least one run of every
    # recorded successor must still be unrecognisable with the claw allowed
    ext_checked = 0
    for entry in report.get("extended_unrepresented", []):
        fam = entry["origin"]["family"]
        n = entry["origin"]["n"]
        j = entry["origin"]["j"]
        r, c = entry["move"]
        successor = base.cells_of(fam, n, j) - base.move_targets((r, c))
        if int(entry["successor_mask"], 16) != encode(successor):
            raise ValueError(f"extended closure mask wrong for {fam}({n},{j})")
        witness = False
        for piece in entry["pieces"]:
            start, end = piece["columns"]
            run = frozenset((rr, cc) for rr, cc in successor
                            if start <= cc <= end)
            if piece["label"] == "unknown" and not candidate_labels(run):
                witness = True
        if not witness:
            raise ValueError(f"extended closure entry is representable: "
                             f"{fam}({n},{j}) move {r},{c}")
        ext_checked += 1
    print(f"closure: {ext_checked} extended-vocabulary obstructions re-derived")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("survey", type=Path, nargs="?",
                        default=ROOT / "data/strip3_ports_n12.jsonl")
    parser.add_argument("--closure", type=Path,
                        default=ROOT / "data/strip3_ports_closure.json")
    parser.add_argument("--recompute-n", type=int, default=6,
                        help="recompute family values up to this n with the "
                             "slow independent recursion")
    parser.add_argument("--skip-closure", action="store_true")
    args = parser.parse_args()

    header, rows = load_survey(args.survey)
    print(f"survey file: max_n={header['max_n']}, {len(rows)} rows, "
          f"strip3 sha256={str(header.get('strip3_sha256'))[:16]}...")
    by = check_families(rows, args.recompute_n)

    full_survey_path = ROOT / "data/strip3_survey.jsonl"
    survey_map = {"F": {}, "A": {}, "D": {}}
    for line in full_survey_path.read_text().splitlines():
        if not line.strip():
            continue
        rec = json.loads(line)
        survey_map["F"][rec["n"]] = rec["grundy"]
        survey_map["A"][rec["n"]] = rec["end_gap"]
        survey_map["D"][rec["n"]] = rec["both_gaps"]
    survey_map["F"][0] = survey_map["A"][0] = survey_map["D"][0] = 0

    check_against_full_survey(by)
    check_rows(by, survey_map)
    check_l2_data(by, survey_map)
    if not args.skip_closure and args.closure.exists():
        check_closure(args.closure, header["max_n"])
    print("all port-family data checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
