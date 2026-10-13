#!/usr/bin/env python3
"""Exhaust the size-independent four-value-block lemma in notes/15.

This verifier imports no production solver. Its transition graph allows
arbitrary moves outside the tested block, so it covers every core whose
direct options miss that block, rather than a finite sample of cores.
"""

import argparse
from collections import deque
import hashlib
import json
from pathlib import Path


def c(phase, n):
    return ((0, 1, 2, 1), (0, 1, 2, 3))[phase][n % 4]


def offsets(phase):
    table = []
    for r in range(4):
        row = []
        for s in range(4):
            ds = {c(phase ^ ((s + 3) % 2), r - s - 3)}
            if phase ^ (s % 2):
                ds.add(c(0, r - s - 1))
            row.append(ds)
        table.append(row)
    return table


def block_lemma(phase):
    # Four 4-bit history masks, q(n-2), q(n-1), n mod 4.
    # Symbol 4 means outside the block; xor 0..3 cannot enter the block.
    table = offsets(phase)
    initial = (0, 0, 0, 0, 4, 4, 3)
    reached = {initial}
    queue = deque([initial])
    edges = 0
    while queue:
        state = queue.popleft()
        history = state[:4]
        a, b, r = state[4:]
        forbidden = set()
        for s, mask in enumerate(history):
            for v in range(4):
                if mask >> v & 1:
                    forbidden.update(v ^ d for d in table[r][s])
        delta = phase ^ (r % 2)
        if a < 4:
            forbidden.add(a ^ delta)
        if not delta and b < 4:
            forbidden.add(b)
        if len(forbidden) == 4:
            raise AssertionError(("the entire block is forbidden", phase, state))
        first_free = min(set(range(4)) - forbidden)
        # The core supplies NO option in the block. Global mex is either
        # outside (below it), or its least tail-unforbidden member.
        for value in (4, first_free):
            updated = list(history)
            if a < 4:
                updated[(r - 2) % 4] |= 1 << a
            successor = (*updated, b, value, (r + 1) % 4)
            edges += 1
            if successor not in reached:
                reached.add(successor)
                queue.append(successor)
    serialized = json.dumps(sorted(reached), separators=(",", ":")).encode()
    return {"phase": phase, "states": len(reached), "transitions": edges,
            "reachable_states_sha256": hashlib.sha256(serialized).hexdigest(),
            "full_forbidden_blocks": 0}


def check_tail_identity():
    """Compare coefficient sets for EVERY prefix, not random SG values.

    Equality of these maps proves full option-set equality for arbitrary
    sequences. It includes both endpoint corrections and the two initial
    recent-prefix cases.
    """
    checked = 0
    for phase in (0, 1):
        for n in range(3, 131):
            literal = [set() for _ in range(n)]
            for i in range(n):
                tooth = phase ^ (i % 2)
                if tooth:
                    literal[i].add(c(phase ^ ((i + 1) % 2), n - i - 1))
                if i:
                    isolated = (1 + int(i + 1 < n)) * (1 - tooth) % 2
                    literal[i - 1].add(c(phase ^ ((i + 2) % 2), max(n - i - 2, 0))
                                       ^ isolated)
            compressed = [set() for _ in range(n)]
            for j in range(n - 2):
                compressed[j].add(c(phase ^ ((j + 3) % 2), n - j - 3))
                if phase ^ (j % 2):
                    compressed[j].add(c(0, n - j - 1))
            delta = phase ^ (n % 2)
            compressed[n - 2].add(delta)
            if not delta:
                compressed[n - 1].add(0)
            assert literal == compressed, (phase, n)
            checked += 1
        table = offsets(phase)
        for r in range(4):
            delta = phase ^ (r % 2)
            if delta:
                assert 1 in table[r][(r - 2) % 4]
            else:
                assert 0 in table[r][(r - 1) % 4]
    return checked


def check_mex_maps():
    # Only the relative order/equality of four integers matters. Testing
    # all ordered pairs from 0..3 covers every weak-order type of them.
    checked = 0
    for a in range(4):
        for b in range(a + 1, 4):
            for u in range(4):
                for v in range(u + 1, 4):
                    def f(x):
                        return b if x == a else a

                    def h(x):
                        return v if x == u else u

                    for x in (a, b):
                        assert f(h(f(h(x)))) == f(h(x)), (a, b, u, v, x)
                    checked += 1
    return checked


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = {"format": "nodekayles-comb-stabilization-lemma-v1",
              "tail_coefficient_identities": check_tail_identity(),
              "mex_map_order_types_checked": check_mex_maps(),
              "block_lemma": [block_lemma(p) for p in (0, 1)],
              "scope": "Universal finite block lemma; infinite consequences follow "
                       "by the proof in notes/15_comb_stabilization.md."}
    raw = json.dumps(report, ensure_ascii=False, indent=2) + "\n"
    if args.output:
        args.output.write_text(raw)
    print(raw, end="")


if __name__ == "__main__":
    main()
