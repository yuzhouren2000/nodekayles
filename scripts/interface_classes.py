#!/usr/bin/env python3
"""Interface classes of the reachable three-row language, measured on masks.

Question: is the Grundy value enough to predict how a connected component
behaves when more columns are attached?  notes/13 answered "no" with one
reachable pair; this tool measures the whole picture at small widths.

For every pair of reachable components X, X' of the same width w with equal
value g(X) = g(X'), the script searches for a right-context Y (again a
reachable component, width <= K) with g(X | Y<<3w) != g(X' | Y<<3w).  Two
components that survive every context of width <= K are "K-equivalent".  The
number of distinct signature rows (g(X | Y<<3w))_Y over contexts of width <= K
is the number of "K-classes"; classes can only refine as K grows, so if the
count keeps growing with K there is no bounded interface on this family.

Values are computed by the ORIGINAL hash-table solver, which accepts arbitrary
induced subgraphs (masks mode) and is the longest cross-validated oracle in
the repository; this keeps the measurement independent of the language
solver's indexing.

Run: python3 scripts/interface_classes.py --objects 5 --contexts 8 --json out.json
"""

import argparse
from collections import defaultdict
from itertools import combinations
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from strip3_history_language import (Code, words, ENDS,  # noqa: E402
                                     endpoint_choices)


def all_codes(n):
    """Every distinct reachable component mask of exactly n columns.

    The family is exactly {(H, L, R)} with L in D_{d_L(H)} and R in
    D_{d_R(H)} (notes/17): combinations outside endpoint_choices describe
    connected masks that are not reachable, so filter them out.
    """
    if n < 3:
        from strip3_history_language import small_catalogue
        return [(mask, code) for mask, code in small_catalogue(n)]
    found = {}
    for word in words(n):
        for left in endpoint_choices(word):
            for right in endpoint_choices(word, True):
                code = Code(word, left, right)
                mask = code.mask()
                if not mask & 7 or mask.bit_length() <= 3 * (n - 1):
                    continue
                found[mask] = code
    return [(mask, found[mask]) for mask in sorted(found)]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--objects", type=int, required=True)
    parser.add_argument("--contexts", type=int, required=True)
    parser.add_argument("--solver", default="src/cfast/strip3")
    parser.add_argument("--json", type=Path, default=None,
                        help="full behaviour matrix (large)")
    parser.add_argument("--summary", type=Path, default=None,
                        help="compact summary: counts, histogram, merged pairs")
    args = parser.parse_args()

    w0 = args.objects
    objects = all_codes(w0)
    contexts = []
    for n in range(1, args.contexts + 1):
        contexts.extend(all_codes(n))
    widths = [len(code.word) for _, code in contexts]
    print(f"objects: {len(objects)} of width {w0}; "
          f"contexts: {len(contexts)} of widths 1..{args.contexts}")

    def values_for(masks, widths_):
        stdin = "".join(f"{n} {m:x}\n" for n, m in zip(widths_, masks))
        result = subprocess.run([args.solver, "masks"], input=stdin, text=True,
                                capture_output=True, cwd=ROOT)
        if result.returncode != 0:
            raise SystemExit(f"solver rc={result.returncode}: {result.stderr[-2000:]}")
        out = [int(line) for line in result.stdout.split()]
        if len(out) != len(masks):
            raise SystemExit("solver output count mismatch")
        return out

    started = time.time()
    object_values = values_for([m for m, _ in objects], [w0] * len(objects))
    print(f"object values: {sorted(set(object_values))}  ({time.time()-started:.1f}s)")

    # full behaviour matrix: row per object over every context
    jobs = [(xi, yi, Xmask | (Ymask << (3 * w0)), w0 + widths[yi])
            for xi, (Xmask, _) in enumerate(objects)
            for yi, (Ymask, _) in enumerate(contexts)]
    matrix_values = values_for([m for _, _, m, _ in jobs], [w for _, _, _, w in jobs])
    print(f"behaviour matrix: {len(jobs)} entries ({time.time()-started:.1f}s)")

    row = {}
    for (xi, yi, _, _), v in zip(jobs, matrix_values):
        row.setdefault(xi, {})[yi] = v

    # signature classes for each K
    class_counts = {}
    for K in range(1, args.contexts + 1):
        subset = [yi for yi in range(len(contexts)) if widths[yi] <= K]
        table = defaultdict(list)
        for xi in range(len(objects)):
            key = tuple(row[xi][yi] for yi in subset)
            table[key].append(xi)
        class_counts[K] = len(table)
    print("K : signature classes (objects of width %d, contexts <= K)" % w0)
    for K in sorted(class_counts):
        print(f"  {K} : {class_counts[K]}")

    by_value = defaultdict(list)
    for xi, v in enumerate(object_values):
        by_value[v].append(xi)
    print(f"value classes: {len(by_value)}")

    # pair stability: how many value-equal pairs survive all contexts <= K
    stable = {}
    for K in range(1, args.contexts + 1):
        subset = [yi for yi in range(len(contexts)) if widths[yi] <= K]
        count = 0
        examples = []
        for v, group in by_value.items():
            for xi, xj in combinations(group, 2):
                if all(row[xi][yi] == row[xj][yi] for yi in subset):
                    count += 1
                    if len(examples) < 3:
                        examples.append((xi, xj))
        stable[K] = {"count": count, "examples": examples}
    print("K : value-equal pairs stable under all contexts <= K")
    for K in sorted(stable):
        print(f"  {K} : {stable[K]['count']}")

    if args.json:
        out = {
            "objects_width": w0,
            "contexts_max_width": args.contexts,
            "objects": [{"mask": f"{m:x}", "word": list(c.word),
                         "left": c.left, "right": c.right,
                         "value": object_values[xi]}
                        for xi, (m, c) in enumerate(objects)],
            "contexts": [{"mask": f"{m:x}", "width": widths[yi]}
                         for yi, (m, _) in enumerate(contexts)],
            "matrix": [[row[xi][yi] for yi in range(len(contexts))]
                       for xi in range(len(objects))],
            "class_counts": {str(k): v for k, v in class_counts.items()},
            "stable_pairs": {str(k): v["count"] for k, v in stable.items()},
            "value_classes": len(by_value),
        }
        args.json.write_text(json.dumps(out) + "\n")
        print(f"wrote {args.json} ({args.json.stat().st_size/1e6:.1f} MB)")

    if args.summary:
        # final classes over all contexts + merged pair witnesses
        final = defaultdict(list)
        for xi in range(len(objects)):
            final[tuple(row[xi][yi] for yi in range(len(contexts)))].append(xi)
        histogram = defaultdict(int)
        merged_pairs = []
        for v in final.values():
            histogram[len(v)] += 1
            for pair in combinations(v, 2):
                merged_pairs.append([objects[pair[0]][0], objects[pair[1]][0]])
        summary = {
            "kind": "interface_classes_summary",
            "objects_width": w0,
            "contexts_max_width": args.contexts,
            "objects": len(objects),
            "object_masks": [f"{m:x}" for m, _ in objects],
            "object_values": object_values,
            "contexts": len(contexts),
            "value_classes": len(by_value),
            "class_counts": {str(k): v for k, v in class_counts.items()},
            "stable_value_pairs": {str(k): v["count"] for k, v in stable.items()},
            "class_size_histogram": {str(k): v for k, v in sorted(histogram.items())},
            "merged_pairs": [[f"{a:x}", f"{b:x}"] for a, b in merged_pairs],
        }
        args.summary.write_text(json.dumps(summary, indent=1) + "\n")
        print(f"wrote {args.summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
