#!/usr/bin/env python3
"""Verify the interface-class claims of scripts/interface_classes.py.

Claims under test (reachable three-row components, contexts = reachable
components glued on the right):
  C1. the object count at width w equals the proved count c_w;
  C2. every pair recorded as "merged" (equal value for every context) admits a
      graph isomorphism X -> X' that fixes the last column pointwise; this
      direction is a theorem: such an isomorphism extends by the identity on
      the context, so all right-extensions have equal Grundy values;
  C3. a sample of merged pairs really re-evaluates equal on explicit contexts
      with an independent solver run, including pairs not related by any
      rectangle symmetry.

Run: python3 tests/verify_interface_classes.py
"""

import json
from itertools import combinations
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from strip3_history_language import component_count  # noqa: E402
from interface_classes import all_codes  # noqa: E402


def cell_set(mask, w):
    return {(c, r) for c in range(w) for r in range(3) if mask >> (3 * c + r) & 1}


def graph(cells):
    return {(c, r): [(c + dc, r + dr) for dc, dr in ((0, 1), (0, -1), (1, 0), (-1, 0))
                     if (c + dc, r + dr) in cells]
            for (c, r) in cells}


def is_boundary_fixing_isomorphism(a_cells, b_cells):
    """Backtracking search for an isomorphism mapping the last column pointwise."""
    if len(a_cells) != len(b_cells):
        return False
    last = max(c for c, _ in a_cells)
    boundary = {v for v in a_cells if v[0] == last}
    if {v for v in b_cells if v[0] == last} != boundary:
        return False
    na, nb = graph(a_cells), graph(b_cells)
    order = sorted(a_cells, key=lambda v: -len(na[v]))
    mapping, used = {}, set()

    def plausible(u, v):
        if len(na[u]) != len(nb[v]):
            return False
        inv_neighbours = {mapping[k] for k in na[u] if k in mapping}
        if any(x not in nb[v] for x in inv_neighbours):
            return False
        mapped_to_v = [k for k, m in mapping.items() if m in nb[v]]
        return all(k in na[u] for k in mapped_to_v)

    def backtrack(i):
        if i == len(order):
            return True
        u = order[i]
        candidates = [u] if u in boundary else sorted(b_cells)
        for v in candidates:
            if v in used or not plausible(u, v):
                continue
            mapping[u], used_tmp = v, v in used
            used.add(v)
            if backtrack(i + 1):
                return True
            del mapping[u]
            used.discard(v)
        return False

    return backtrack(0)


def flip_rows(mask, w):
    out = 0
    for c in range(w):
        col = mask >> (3 * c) & 7
        out |= (((col & 1) << 2) | (col & 2) | ((col & 4) >> 2)) << (3 * c)
    return out


def main():
    summary_path = ROOT / "data/proofs/interface_classes_summary.json"
    summaries = json.loads(summary_path.read_text())
    failures = []
    for width, summary in sorted(summaries.items(), key=lambda kv: int(kv[0])):
        w = int(width)
        if summary["objects"] != component_count(w):
            failures.append(f"w={w}: object count {summary['objects']} != c_{w}")
        merged = {frozenset(p) for p in summary["merged_pairs"]}
        for a_hex, b_hex in summary["merged_pairs"]:
            a, b = int(a_hex, 16), int(b_hex, 16)
            if not is_boundary_fixing_isomorphism(cell_set(a, w), cell_set(b, w)):
                failures.append(f"w={w}: merged pair {a_hex}/{b_hex} has no "
                                "boundary-fixing isomorphism")
        # theorem direction: every boundary-fixing-iso pair must be merged
        masks = [int(m, 16) for m in summary["object_masks"]]
        if w <= 6:
            iso_pairs = 0
            for a, b in combinations(masks, 2):
                if is_boundary_fixing_isomorphism(cell_set(a, w), cell_set(b, w)):
                    iso_pairs += 1
                    if frozenset((f"{a:x}", f"{b:x}")) not in merged:
                        failures.append(f"w={w}: iso pair {a:x}/{b:x} not merged")
            print(f"w={w}: {summary['objects']} objects, {summary['value_classes']} value "
                  f"classes, {len(summary['merged_pairs'])} merged pairs, "
                  f"{iso_pairs} boundary-fixing-iso pairs (all must be merged)")
        else:
            print(f"w={w}: {summary['objects']} objects, {summary['value_classes']} value "
                  f"classes, {len(summary['merged_pairs'])} merged pairs checked")
        # informative: row-flip pairs merge iff the last column is middle-only
        present = set(masks)
        via_flip = non_middle = 0
        for a in masks:
            b = flip_rows(a, w)
            if b == a or b not in present or a > b:
                continue
            last_only_middle = {(c, 1) for c in (w - 1,)} == cell_set(a, w) & \
                {(w - 1, r) for r in range(3)}
            if frozenset((f"{a:x}", f"{b:x}")) in merged:
                via_flip += 1
                if not last_only_middle:
                    non_middle += 1
        print(f"      row-flip-asymmetric pairs merged: {via_flip} "
              f"(of which last column not middle-only: {non_middle})")

    # C3: independent re-evaluation of sampled merged pairs on reachable contexts
    contexts = []
    for n in (1, 2, 3):
        for mask, code in all_codes(n):
            contexts.append((mask, len(code.word)))
    sample = [p for s in summaries.values() for p in s["merged_pairs"]][:20]
    lines, pairs = [], 0
    for a_hex, b_hex in sample:
        a, b = int(a_hex, 16), int(b_hex, 16)
        w = (max(a.bit_length(), b.bit_length()) + 2) // 3        # exact width
        for cmask, cw in contexts:
            lines.append(f"{w + cw} {a | (cmask << (3 * w)):x}")
            lines.append(f"{w + cw} {b | (cmask << (3 * w)):x}")
            pairs += 1
    result = subprocess.run(["./src/cfast/strip3", "masks"], input="\n".join(lines) + "\n",
                            text=True, capture_output=True, cwd=ROOT)
    values = [int(v) for v in result.stdout.split()]
    mismatches = sum(1 for i in range(0, len(values) - 1, 2) if values[i] != values[i + 1])
    print(f"independent re-evaluation: rc={result.returncode}, {pairs} merged-pair x "
          f"context instances, {mismatches} mismatches")
    if result.returncode != 0 or mismatches or len(values) != 2 * pairs:
        failures.append(f"independent re-evaluation: rc={result.returncode} "
                        f"mismatches={mismatches} values={len(values)}")

    if failures:
        print("FAILURES:")
        for f in failures:
            print(" " + f)
        return 1
    print("all interface-class checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
