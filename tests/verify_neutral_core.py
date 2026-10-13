#!/usr/bin/env python3
"""Check a neutral-core proof without importing its builder or SG solver."""

import json
from pathlib import Path
import sys


def verify(proof):
    if proof.get("format") == "node-kayles-neutral-obstruction-v1":
        return verify_obstruction(proof)
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


def verify_obstruction(proof):
    """Replay an actual history and check a star counterexample, without SG search."""
    graph, witness = proof["adjacency"], proof["witness"]
    n = len(graph)
    full = (1 << n) - 1
    if any(type(row) is not int or not 0 <= row <= full for row in graph):
        raise ValueError("invalid graph")
    for v, row in enumerate(graph):
        if row >> v & 1 or any(bool(row >> u & 1) != bool(graph[u] >> v & 1)
                               for u in range(n)):
            raise ValueError("invalid simple graph")
    port = witness["port"]
    if type(port) is not int or not 0 <= port < n:
        raise ValueError("invalid port")
    mask = full
    for v in witness["legal_history"]:
        if type(v) is not int or not 0 <= v < n or not mask >> v & 1:
            raise ValueError("illegal history")
        mask &= ~(graph[v] | 1 << v)
        if not mask >> port & 1:
            raise ValueError("history deletes the port")
    if type(witness["terminal_mask"]) is not int or mask != witness["terminal_mask"]:
        raise ValueError("incorrect terminal state")
    leaves = [v for v in range(n) if graph[v] == 1 << port]
    if list(witness["original_leaves"]) != leaves or len(leaves) % 2:
        raise ValueError("expected an even original pendant set")
    leaf_mask = sum(1 << v for v in leaves)
    if graph[port] & mask != leaf_mask or leaf_mask & ~mask:
        raise ValueError("terminal port component is not the claimed star")
    others = mask & ~(leaf_mask | 1 << port)
    if any(graph[v] & mask for v in range(n) if others >> v & 1):
        raise ValueError("remaining vertices are not isolates")
    offset = others.bit_count() % 2
    before = (2 if leaves else 1) ^ offset
    after = offset
    if (type(witness["terminal_grundy"]) is not int
            or type(witness["deleted_port_grundy"]) is not int
            or witness["terminal_grundy"] != before
            or witness["deleted_port_grundy"] != after):
        raise ValueError("incorrect star values")
    return {"terminal_grundy": before, "deleted_port_grundy": after}


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
    print(json.dumps(verify(json.loads(Path(sys.argv[1]).read_text()))))
