#!/usr/bin/env python3
"""Check a neutral-core proof without importing its builder or SG solver."""

import json
from pathlib import Path
import sys


def verify(proof):
    if proof.get("format") == "node-kayles-paired-neutral-core-v1":
        return verify_paired(proof)
    if proof.get("format") != "node-kayles-neutral-core-v1":
        raise ValueError("unknown certificate format")
    graph, ports = proof["adjacency"], proof["ports"]
    n = len(graph)
    full = (1 << n) - 1
    if (len(set(ports)) != len(ports)
            or any(type(v) is not int or not 0 <= v < n for v in ports)):
        raise ValueError("invalid ports")
    for v, neighbours in enumerate(graph):
        if (type(neighbours) is not int or neighbours < 0 or neighbours & ~full
                or neighbours >> v & 1):
            raise ValueError("invalid simple graph")
        if any(bool(neighbours >> u & 1) != bool(graph[u] >> v & 1) for u in range(n)):
            raise ValueError("asymmetric graph")
    states = dict(proof["states"])
    if type(proof["grundy"]) is not int or proof["grundy"] < 0:
        raise ValueError("invalid root value")
    if len(states) != len(proof["states"]):
        raise ValueError("duplicate states")
    if any(type(m) is not int or not 0 <= m <= full or type(g) is not int or g < 0
           for m, g in states.items()):
        raise ValueError("invalid state/value")
    pending, seen = [full], set()
    while pending:
        mask = pending.pop()
        if mask in seen:
            continue
        seen.add(mask)
        if mask not in states:
            raise ValueError("missing reachable state")
        options = set()
        for v in range(n):
            if mask >> v & 1:
                child = mask & ~(graph[v] | 1 << v)
                if child not in states:
                    raise ValueError("missing actual option")
                options.add(states[child])
                pending.append(child)
        actual = 0
        while actual in options:
            actual += 1
        if states[mask] != actual:
            raise ValueError("incorrect mex")
        for port in ports:
            if mask >> port & 1:
                child = mask & ~(1 << port)
                if child not in states or states[child] != states[mask]:
                    raise ValueError("port deletion is not neutral")
                pending.append(child)
    if seen != set(states) or proof["grundy"] != states[full]:
        raise ValueError("extra states or incorrect root")
    return proof["grundy"]


def verify_paired(proof):
    """Check the exterior geometry and extract its core independently."""
    graph, sigma = proof["adjacency"], proof["involution"]
    n = len(graph)
    full = (1 << n) - 1
    if (len(sigma) != n or any(type(v) is not int for v in sigma)
            or set(sigma) != set(range(n))
            or any(sigma[sigma[v]] != v for v in range(n))):
        raise ValueError("invalid pairing")
    if any(type(row) is not int or not 0 <= row <= full for row in graph):
        raise ValueError("invalid graph")
    for v, row in enumerate(graph):
        if row >> v & 1 or any(bool(row >> u & 1) != bool(graph[u] >> v & 1)
                               for u in range(n)):
            raise ValueError("invalid simple graph")
    fixed = [v for v in range(n) if sigma[v] == v]
    if list(proof["fixed_vertices"]) != fixed:
        raise ValueError("incorrect fixed core")
    for v in range(n):
        mate = sigma[v]
        if mate == v:
            mask = graph[v]
        else:
            if graph[v] >> mate & 1:
                raise ValueError("adjacent partners")
            mask = graph[v] | graph[mate] | 1 << v | 1 << mate
        if any(bool(mask >> u & 1) != bool(mask >> sigma[u] & 1)
               for u in range(n)):
            raise ValueError("incompatible pairing")
    outside = full ^ sum(1 << v for v in fixed)
    ports = [i for i, v in enumerate(fixed) if graph[v] & outside]
    expected_graph = [sum(1 << j for j, u in enumerate(fixed) if graph[v] >> u & 1)
                      for v in fixed]
    core = proof["core"]
    if (core.get("format") != "node-kayles-neutral-core-v1"
            or list(core["adjacency"]) != expected_graph
            or list(core["ports"]) != ports
            or list(proof["core_ports"]) != [fixed[i] for i in ports]):
        raise ValueError("incorrect extracted core")
    value = verify(core)
    if type(proof["grundy"]) is not int or proof["grundy"] != value:
        raise ValueError("incorrect paired value")
    return value


if __name__ == "__main__":
    print(verify(json.loads(Path(sys.argv[1]).read_text())))
