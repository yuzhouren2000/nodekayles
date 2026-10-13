#!/usr/bin/env python3
"""Build P-position certificates from short two-move separators.

No full-board or B(n,j) value is requested.  The saved certificate consists of
legal responses, disjoint geometric pieces, and the SG premises for the pieces.
See notes/12_two_move_separators.md for the parameter-independent identities.
"""

import argparse
import json
from pathlib import Path

from strip3_ports import (cells_F, closed_cells, family_cells, run_masks,
                          sha256_file, to_mask)

ROOT = Path(__file__).resolve().parents[1]
FORMAT = "nodekayles-strip3-two-move-v1"


def shape(family, n, j=None):
    if family == "point":
        return {(2, 1)}
    if family == "T":
        return cells_F(n) - {(1, n)}
    if family == "C":
        return cells_F(n) - {(1, n), (3, n)}
    return family_cells(family, n, j)


def load_values(extension):
    values = {("A", 0, None): 0, ("D", 0, None): 0,
              ("T", 0, None): 0, ("C", 0, None): 0,
              ("point", 1, None): 1}
    for line in (ROOT / "data/strip3_survey.jsonl").read_text().splitlines():
        row = json.loads(line)
        n = row["n"]
        values["A", n, None] = row["end_gap"]
        values["D", n, None] = row["both_gaps"]
        values["B", n, 1] = row["options"][0]
    for line in (ROOT / "data/strip3_ports_n12.jsonl").read_text().splitlines():
        row = json.loads(line)
        if row.get("family") in ("E", "H"):
            values[row["family"], row["n"], row["j"]] = row["value"]
    for path in (extension, ROOT / "data/strip3_two_move_auxiliary.json"):
        if not path.exists():
            continue
        for row in json.loads(path.read_text())["rows"]:
            if type(row["value"]) is int and row["value"] >= 0:
                values[row["family"], row["n"], row.get("j")] = row["value"]
    return values


def descriptor(family, n, j=None):
    return family, n, j


def middle_pair(n, i, k):
    if abs(i - k) < 2:
        return None
    lo, hi = sorted((i, k))
    return [(descriptor("A", lo - 1), 1, lo - 1),
            (descriptor("D", hi - lo - 1), lo + 1, hi - 1),
            (descriptor("A", n - hi), hi + 1, n)]


def mixed_pair(n, j, i):
    if i == j:
        return None
    if i < j:
        return [(descriptor("A", i - 1), 1, i - 1),
                (descriptor("E", n - i, j - i), i + 1, n)]
    return [(descriptor("E", i - 1, i - j), 1, i - 1),
            (descriptor("A", n - i), i + 1, n)]


def opposite_pair(n, j, i):
    """All five opposite-edge separators.  No formula for a distant reply."""
    delta = abs(j - i)
    if delta > 2:
        return None
    lo, hi = sorted((i, j))
    if delta == 0:
        return [(descriptor("C", j - 1), 1, j - 1),
                (descriptor("C", n - j), j + 1, n)]
    if delta == 1:
        return [(descriptor("T", lo - 1), 1, lo - 1),
                (descriptor("T", n - hi), hi + 1, n)]
    return [(descriptor("B", lo, 1), 1, lo),
            (descriptor("point", 1), lo + 1, lo + 1),
            (descriptor("B", n - hi + 1, 1), hi, n)]


def reflected_match(cells, model):
    """Find the explicit isomorphism from the model into an actual piece."""
    if not cells and not model:
        return {"column_offset": 0, "horizontal": False, "vertical": False}
    if not cells or not model:
        return None
    offset = min(c for _, c in cells) - 1
    width = max(c for _, c in cells) - offset
    for horizontal in (False, True):
        for vertical in (False, True):
            mapped = {(4 - r if vertical else r,
                       offset + (width + 1 - c if horizontal else c))
                      for r, c in model}
            if mapped == cells:
                return {"column_offset": offset, "horizontal": horizontal,
                        "vertical": vertical}
    return None


def candidate(n, first, reply, parts, rule, values):
    if parts is None:
        return None
    total = 0
    for key, _, _ in parts:
        if key not in values:
            return None
        total ^= values[key]
    if total:
        return None
    full = cells_F(n)
    child = full - closed_cells(first)
    if reply not in child:
        return None
    tail = child - closed_cells(reply)
    used, pieces = set(), []
    for key, left, right in parts:
        family, length, j = key
        cells = {v for v in tail if left <= v[1] <= right}
        mapping = reflected_match(cells, shape(family, length, j))
        if mapping is None:
            raise AssertionError((n, first, reply, key, cells))
        used |= cells
        if not cells:
            continue
        pieces.append({"family": family, "n": length, "j": j,
                       "value": values[key], "mask": f"{to_mask(cells):x}",
                       "isomorphism": mapping})
    if used != tail:
        raise AssertionError("separator missed vertices")
    if any(any(v in closed_cells(u) for v in b) for k, a in enumerate(
            [{v for v in tail if l <= v[1] <= r} for _, l, r in parts])
           for b in [{v for v in tail if l <= v[1] <= r}
                     for _, l, r in parts[k + 1:]] for u in a):
        raise AssertionError("separator has a crossing edge")
    return {"first": list(first), "reply": list(reply), "rule": rule,
            "tail_mask": f"{to_mask(tail):x}", "pieces": pieces}


def safe_rotation(n, first):
    reply = (4 - first[0], n + 1 - first[1])
    child = cells_F(n) - closed_cells(first)
    if reply not in child:
        return None
    tail = child - closed_cells(reply)
    rotate = lambda v: (4 - v[0], n + 1 - v[1])
    if any(rotate(v) not in tail or rotate(v) in closed_cells(v) for v in tail):
        return None
    return {"first": list(first), "reply": list(reply), "rule": "safe_rotation",
            "tail_mask": f"{to_mask(tail):x}", "pieces": []}


def cost(case):
    widths = [p["n"] for p in case["pieces"]]
    return (max(widths, default=0), sum(widths),
            sum(p["family"] == "E" for p in case["pieces"]),
            case["reply"])


def build(n, values):
    if n < 2 or n % 2:
        raise ValueError("this certificate builder expects a positive even width")
    cases, unresolved = [], []
    for row in (1, 2):
        for j in range(1, n // 2 + 1):
            first = (row, j)
            options = []
            rot = safe_rotation(n, first)
            if rot:
                options.append(rot)
            for i in range(1, n + 1):
                if row == 1:
                    specs = [((2, i), mixed_pair(n, j, i), "outer_middle"),
                             ((3, i), opposite_pair(n, j, i), "opposite_outer")]
                else:
                    specs = [((2, i), middle_pair(n, j, i), "middle_middle"),
                             ((1, i), mixed_pair(n, i, j), "middle_outer")]
                for reply, parts, rule in specs:
                    result = candidate(n, first, reply, parts, rule, values)
                    if result:
                        options.append(result)
            if options:
                cases.append(min(options, key=cost))
            else:
                unresolved.append(list(first))
    premises = {}
    for case in cases:
        for piece in case["pieces"]:
            key = piece["family"], piece["n"], piece["j"]
            premises[key] = {"family": key[0], "n": key[1], "j": key[2],
                             "value": piece["value"],
                             "mask": f"{to_mask(shape(*key)):x}"}
    return {"format": FORMAT, "n": n,
            "status": "P_CERTIFICATE" if not unresolved else "PARTIAL",
            "cases": cases, "unresolved": unresolved,
            "premises": sorted(premises.values(), key=lambda p: (
                p["n"], p["family"], p["j"] or 0)),
            "max_premise_width": max((k[1] for k in premises), default=0)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--n", type=int)
    mode.add_argument("--recompute-values", type=Path,
                      help="recompute every listed auxiliary descriptor with strip3 masks")
    parser.add_argument("--extension", type=Path,
                        default=ROOT / "data/strip3_ports_extension.json")
    parser.add_argument("--ht-bits", type=int, default=26)
    parser.add_argument("--timeout", type=float, default=300)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.recompute_values:
        rows = [{k: p.get(k) for k in ("family", "n", "j")}
                for p in json.loads(args.recompute_values.read_text())["rows"]]
        for row in rows:
            row["mask"] = f"{to_mask(shape(row['family'], row['n'], row['j'])):x}"
        result, info = run_masks([(p["n"], int(p["mask"], 16)) for p in rows],
                                 ht_bits=args.ht_bits, timeout=args.timeout)
        for row, value in zip(rows, result):
            row["value"] = value
        report = {"format": "nodekayles-strip3-auxiliary-values-v1",
                  "strip3_sha256": sha256_file(ROOT / "src/cfast/strip3.c"),
                  "input_sha256": sha256_file(args.recompute_values),
                  "info": info, "rows": rows}
        args.output.write_text(json.dumps(report, indent=2) + "\n")
        print(f"auxiliary recomputation: {info['emitted']}/{len(rows)} exact values; "
              f"{info['status']}; {info['elapsed_s']}s")
        return
    report = build(args.n, load_values(args.extension))
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(f"3x{args.n}: {report['status']}, {len(report['cases'])}/{args.n} "
          f"first-move orbits have zero responses; "
          f"{len(report['premises'])} premises, width <= {report['max_premise_width']}")
    if report["unresolved"]:
        print("unresolved:", report["unresolved"])


if __name__ == "__main__":
    main()
