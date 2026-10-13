#!/usr/bin/env python3
"""Search exact interface equalities in the actual reachable three-row family.

Only complete endpoint columns are exposed. A matching recursive interface
is a sufficient arbitrary-context certificate; failure to match is not a
proof of semantic inequivalence. Label-preserving graph isomorphisms are
filtered before reporting candidates. No third-party packages are required.
"""

import argparse
from collections import Counter
from itertools import product
import json
from pathlib import Path
import subprocess
import time

from kayles_boundary_algebra import BoundaryAlgebra
from kayles_boundary_transfer import column_graph
from strip3_history_language import (Code, endpoint_choices, small_catalogue,
                                    words)


def reachable_columns(width):
    if width < 3:
        masks = (mask for mask, _ in small_catalogue(width))
    else:
        masks = (Code(word, left, right).mask() for word in words(width)
                 for left, right in product(endpoint_choices(word),
                                            endpoint_choices(word, True)))
    for mask in masks:
        yield tuple(mask >> (3 * c) & 7 for c in range(width))


def boundary_graph(columns, boundary):
    graph, cells = column_graph(columns)
    last = len(columns) - 1
    labels = tuple((r + 1 if c == last else 0) if boundary == "right"
                   else (r + 1 if c == 0 else r + 4 if c == last else 0)
                   for r, c in cells)
    return graph, labels


def labelled_isomorphic(graph_a, labels_a, graph_b, labels_b):
    """Joint color refinement followed by complete adjacency backtracking."""
    n = len(graph_a)
    if n != len(graph_b) or Counter(labels_a) != Counter(labels_b):
        return False
    graphs = (graph_a, graph_b)
    keys = [(labels_a[v], graph_a[v].bit_count()) for v in range(n)]
    keys += [(labels_b[v], graph_b[v].bit_count()) for v in range(n)]
    palette = {key: i for i, key in enumerate(sorted(set(keys)))}
    colors = [[palette[key] for key in keys[:n]],
              [palette[key] for key in keys[n:]]]
    while True:
        keys = [(colors[s][v], tuple(sorted(colors[s][u] for u in range(n)
                                           if graphs[s][v] >> u & 1)))
                for s in range(2) for v in range(n)]
        palette = {key: i for i, key in enumerate(sorted(set(keys)))}
        refined = [[palette[key] for key in keys[:n]],
                   [palette[key] for key in keys[n:]]]
        if Counter(refined[0]) != Counter(refined[1]):
            return False
        stable = len(set(refined[0] + refined[1])) == len(set(colors[0] + colors[1]))
        colors = refined
        if stable:
            break
    mapping, used = {}, set()

    def candidates(v):
        return [u for u in range(n) if u not in used and colors[0][v] == colors[1][u]
                and all(bool(graph_a[v] >> x & 1) == bool(graph_b[u] >> y & 1)
                        for x, y in mapping.items())]

    def visit():
        if len(mapping) == n:
            return True
        v, choices = min(((v, candidates(v)) for v in range(n) if v not in mapping),
                         key=lambda item: len(item[1]))
        for u in choices:
            mapping[v] = u
            used.add(u)
            if visit():
                return True
            used.remove(u)
            del mapping[v]
        return False

    return visit()


def search(max_width, boundary, proof_dir=None, progress=True):
    algebra = BoundaryAlgebra()
    representatives, matches, rows = {}, [], []
    examined = isomorphic = 0
    started = time.monotonic()
    for width in range(2, max_width + 1):
        count = 0
        for columns in reachable_columns(width):
            if columns[-1] != 7 or boundary == "both" and columns[0] != 7:
                continue
            graph, labels = boundary_graph(columns, boundary)
            state = algebra.state(graph, labels)
            count += 1
            group = representatives.setdefault(state, [])
            if any(labelled_isomorphic(graph, labels, *other) for _, other in group):
                isomorphic += 1
                continue
            if group:
                previous, _ = group[0]
                match = {"columns_a": previous, "columns_b": columns,
                         "vertices_a": len(group[0][1][0]), "vertices_b": len(graph),
                         "state": state}
                if proof_dir is not None:
                    proof_dir.mkdir(parents=True, exist_ok=True)
                    certificate = BoundaryAlgebra().proof((group[0][1], (graph, labels)))
                    name = f"reachable_{boundary}_{len(matches):03d}.json"
                    (proof_dir / name).write_text(json.dumps(certificate) + "\n")
                    match["certificate"] = name
                matches.append(match)
            group.append((columns, (graph, labels)))
        examined += count
        row = {"width": width, "graphs": count, "total_graphs": examined,
               "interface_nodes": len(algebra.nodes), "isomorphic_duplicates": isomorphic,
               "nonisomorphic_matches": len(matches),
               "elapsed_seconds": round(time.monotonic() - started, 3)}
        rows.append(row)
        if progress:
            print(json.dumps(row), flush=True)
    return {"format": "reachable-full-boundary-search-v1", "boundary": boundary,
            "max_width": max_width, "rows": rows, "matches": matches,
            "scope": "exact recursive equality only; no converse equivalence claim"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-width", type=int, default=8)
    parser.add_argument("--boundary", choices=("right", "both"), default="right")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--proof-dir", type=Path)
    args = parser.parse_args()
    if args.max_width < 2:
        parser.error("--max-width must be at least 2")
    report = search(args.max_width, args.boundary, args.proof_dir)
    try:
        report["git_revision"] = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True, stderr=subprocess.DEVNULL).strip()
    except (FileNotFoundError, subprocess.CalledProcessError):
        report["git_revision"] = None
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
