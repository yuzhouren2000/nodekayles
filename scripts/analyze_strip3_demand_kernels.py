#!/usr/bin/env python3
"""Count saved proof regions without rerunning a game solver.

Semantic correctness is checked separately by the standalone C checker. This
analysis checks recorded hashes, counts satisfying assignments, and compares
width-six counts with independent vertex-set enumeration, not BDD traversal.
"""
import functools
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import resource
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def inspect(data):
    lines = data.decode("ascii").splitlines()
    history = lines[0] == "node-kayles-history-kernel-v1"
    assert history or lines[0] == "node-kayles-symbolic-kernel-v1"
    n, target, size = map(int, lines[1].split())
    cells = 3*n
    graph = [(cells, 0, 0)]+[tuple(map(int, line.split())) for line in lines[4:]]
    domain = int(lines[2]); kernels = list(map(int, lines[3].split()))
    assert len(graph) == size and len(kernels) == target+1

    @functools.lru_cache(None)
    def suffix(ref):
        variable, low, high = graph[ref >> 1]
        return count(low, variable+1)+count(high, variable+1)

    def count(ref, start=0):
        if ref < 2:
            return ref << (cells-start)
        if ref & 1:
            return (1 << (cells-start))-count(ref ^ 1, start)
        return suffix(ref) << (graph[ref >> 1][0]-start)

    def member(ref, assignment):
        while ref >= 2:
            variable, low, high = graph[ref >> 1]
            ref = (high if assignment >> variable & 1 else low) ^ (ref & 1)
        return ref

    domain_size = count(domain)
    kernel_sizes = [count(ref) for ref in kernels]
    if history:
        # Derive valid single-column choices from vertical edges.
        choices = [s for s in range(8) if not s & (s << 1)]
        ways = {0: 1}
        for _ in range(n):
            ways = {s: sum(total for left, total in ways.items() if not left & s) for s in choices}
        assert sum(ways.values()) == domain_size
    enumerated = None
    if n <= 6:
        # Independent chosen-vertex enumeration. No column language or BDD
        # counting is used here; evaluate the saved roots on actual histories.
        closed = []
        for v in range(cells):
            closed.append(sum(1 << u for u in range(cells)
                              if abs(u//3-v//3)+abs(u%3-v%3) <= 1))
        members = set()

        def visit(candidates, chosen, available):
            members.add(chosen if history else available)
            while candidates:
                bit = candidates & -candidates
                candidates ^= bit
                v = bit.bit_length()-1
                visit(candidates & ~closed[v], chosen | bit, available & ~closed[v])

        visit((1 << cells)-1, 0, (1 << cells)-1)
        assert len(members) == domain_size
        assert all(member(domain, mask) for mask in members)
        assert [sum(member(ref, mask) for mask in members) for ref in kernels] == kernel_sizes
        enumerated = len(members)
    return {"n": n, "SG": target, "domain": "independent_histories" if history else "reachable_masks",
            "domain_size": domain_size, "kernel_sizes": kernel_sizes, "covered_size": sum(kernel_sizes),
            "covered_fraction": round(sum(kernel_sizes)/domain_size, 8),
            "small_domain_independent_enumeration": enumerated}


def main():
    os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
    resource.setrlimit(resource.RLIMIT_AS, (512*1024**2,)*2)
    resource.setrlimit(resource.RLIMIT_CPU, (2,)*2)
    report = {"format": "strip3-demand-kernel-coverage-v1", "script_sha256": sha(Path(__file__).read_bytes()),
              "method": "exact BDD model counts, with independent width-six vertex enumeration",
              "semantics": "verified separately in data/proofs/strip3_demand_kernel_audit.json", "samples": []}
    for name in ("formula", "availability", "history"):
        screen_path = ROOT/f"data/strip3_{name}_demand_screen.json"
        raw = screen_path.read_bytes(); screen = json.loads(raw)
        for path, digest in screen["source_sha256"].items():
            frozen = subprocess.run(["git", "show", f"{screen['revision']}:{path}"], cwd=ROOT,
                                    capture_output=True, check=True).stdout
            assert sha(frozen) == digest, path
        for sample in screen["samples"]:
            if "certificate" not in sample:
                continue
            record = sample["certificate"]
            packed = (ROOT/record["path"]).read_bytes(); data = gzip.decompress(packed)
            assert sha(packed) == record["file_sha256"] and sha(data) == record["sha256"]
            result = inspect(data)
            # Last line gives completed counters, earlier lines are checkpoints.
            counts = dict((key, int(value)) for key, value in re.findall(r"(\w+)=(\d+)", sample["solve"]["stderr"]))
            result.update(certificate=record["path"], certificate_sha256=record["sha256"],
                          screen=str(screen_path.relative_to(ROOT)), screen_sha256=sha(raw),
                          solve_counts=counts)
            report["samples"].append(result)
    (ROOT/"data/strip3_demand_kernel_coverage.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
