#!/usr/bin/env python3
"""Cancel paired Node Kayles exteriors using deletion-neutral core ports.

A finite certificate covers every core state reachable by actual moves or
direct port deletion. Deleting any surviving port must preserve its SG.
The exterior needs a valid fixed pairing, not a graph automorphism.
Only the core is evaluated; no exterior SG computation is used.
"""

import argparse
from functools import cache
import json
from pathlib import Path

from strip3_defect_algebra import validate_core

FORMAT = "node-kayles-neutral-core-v1"
PAIRED_FORMAT = "node-kayles-paired-neutral-core-v1"


def _graph(adjacency):
    graph = tuple(adjacency)
    if graph:
        validate_core(graph)
    return graph


def neutral_core_certificate(adjacency, ports):
    """Return a complete finite certificate, or None if neutrality fails.

    A port is an ordinary playable vertex. Its direct deletion is a query
    used to model attacks from the exterior, never an additional game move.
    Certificate construction can be exponential in the core size.
    """
    graph, ports = _graph(adjacency), tuple(ports)
    n = len(graph)
    if (any(type(v) is not int or not 0 <= v < n for v in ports)
            or len(set(ports)) != len(ports)):
        raise ValueError("ports must be distinct graph vertices")
    port_mask = sum(1 << v for v in ports)
    closed = tuple(neighbours | 1 << v for v, neighbours in enumerate(graph))

    @cache
    def value(mask):
        outcomes = {value(mask & ~closed[v]) for v in range(n) if mask >> v & 1}
        result = 0
        while result in outcomes:
            result += 1
        return result

    full = (1 << n) - 1
    pending, seen = [full], set()
    while pending:
        mask = pending.pop()
        if mask in seen:
            continue
        seen.add(mask)
        for v in range(n):
            if not mask >> v & 1:
                continue
            pending.append(mask & ~closed[v])
            if port_mask >> v & 1:
                child = mask & ~(1 << v)
                if value(child) != value(mask):
                    return None
                pending.append(child)
    return {"format": FORMAT, "adjacency": graph, "ports": tuple(sorted(ports)),
            "grundy": value(full),
            "states": tuple((mask, value(mask)) for mask in sorted(seen))}


def paired_core_certificate(adjacency, involution):
    """Recognise a neutral fixed core with a compatible paired exterior.

    The permutation need not preserve individual edges. Each pair must be
    nonadjacent, its closed-neighbourhood union must be invariant, and each
    fixed core vertex must see either both or neither member of every pair.
    """
    graph, sigma = _graph(adjacency), tuple(involution)
    n = len(graph)
    if (len(sigma) != n or any(type(v) is not int for v in sigma)
            or set(sigma) != set(range(n))
            or any(sigma[sigma[v]] != v for v in range(n))):
        raise ValueError("an involutive permutation is required")

    def invariant(mask):
        return all(bool(mask >> v & 1) == bool(mask >> sigma[v] & 1)
                   for v in range(n))

    fixed = tuple(v for v in range(n) if sigma[v] == v)
    fixed_mask = sum(1 << v for v in fixed)
    for v, neighbours in enumerate(graph):
        if sigma[v] == v:
            if not invariant(neighbours):
                return None
        elif neighbours >> sigma[v] & 1:
            return None
        else:
            union = neighbours | graph[sigma[v]] | 1 << v | 1 << sigma[v]
            if not invariant(union):
                return None
    core = tuple(sum(1 << j for j, u in enumerate(fixed) if graph[v] >> u & 1)
                 for v in fixed)
    ports = tuple(i for i, v in enumerate(fixed) if graph[v] & ~fixed_mask)
    certificate = neutral_core_certificate(core, ports)
    if certificate is None:
        return None
    return {"format": PAIRED_FORMAT, "adjacency": graph, "involution": sigma,
            "grundy": certificate["grundy"], "fixed_vertices": fixed,
            "core_ports": tuple(fixed[i] for i in ports), "core": certificate}


def leaf_support_reduction(adjacency):
    """Find a pendant odd star whose sole attachment has another leaf.

    Return (remaining graph, xor offset 1, deleted vertices), or None. This
    structural reduction evaluates no SG, including that of the remainder.
    """
    graph = _graph(adjacency)
    for root, neighbours in enumerate(graph):
        leaves = tuple(v for v in range(len(graph)) if graph[v] == 1 << root)
        if not leaves or len(leaves) % 2 == 0:
            continue
        attachments = neighbours & ~sum(1 << v for v in leaves)
        if attachments.bit_count() != 1:
            continue
        support = attachments.bit_length() - 1
        if not any(v != root and graph[v] == 1 << support for v in range(len(graph))):
            continue
        deleted = (root,) + leaves
        keep = tuple(v for v in range(len(graph)) if v not in deleted)
        reduced = tuple(sum(1 << j for j, u in enumerate(keep) if graph[v] >> u & 1)
                        for v in keep)
        return reduced, 1, deleted
    return None


def reflection_involution(columns, axis="horizontal"):
    """Return graph and a grid reflection, or None for a noninvariant mask."""
    from kayles_boundary_transfer import column_graph
    if axis not in ("horizontal", "vertical"):
        raise ValueError("unknown reflection axis")
    columns = tuple(columns)
    graph, cells = column_graph(columns)
    indices = {cell: i for i, cell in enumerate(cells)}
    mates = tuple((r, len(columns) - c - 1) if axis == "horizontal" else (2 - r, c)
                  for r, c in cells)
    if any(cell not in indices for cell in mates):
        return None
    return graph, tuple(indices[cell] for cell in mates)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--core-word", help="three-row core, column digits 0..7")
    mode.add_argument("--paired-word", help="check a reflected three-row position")
    parser.add_argument("--ports", type=int, nargs="*", default=(),
                        help="zero-based core vertex indices, ordered by column then row")
    parser.add_argument("--axis", choices=("horizontal", "vertical"), default="horizontal")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    word = args.core_word if args.core_word is not None else args.paired_word
    if any(ch not in "01234567" for ch in word):
        parser.error("column words must use digits 0..7")
    columns = tuple(map(int, word))
    if args.core_word is not None:
        from kayles_boundary_transfer import column_graph
        graph, _ = column_graph(columns)
        certificate = neutral_core_certificate(graph, args.ports)
    else:
        if args.ports:
            parser.error("--ports applies only to --core-word")
        reflected = reflection_involution(columns, args.axis)
        certificate = paired_core_certificate(*reflected) if reflected is not None else None
    if args.output:
        args.output.write_text(json.dumps(certificate, indent=2) + "\n")
    core = certificate.get("core", certificate) if certificate is not None else None
    print(json.dumps({"certified": certificate is not None,
                      "grundy": core["grundy"] if core else None,
                      "core_states": len(core["states"]) if core else None}))


if __name__ == "__main__":
    main()
