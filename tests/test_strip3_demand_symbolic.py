"""Independent subset-mex and certificate checks for demand proof learning.

Histories are enumerated as vertex sets and converted to remaining graphs;
no BDD code or producer transition helper is imported into the reference.
Every child process has the local 2 s / 512 MiB limit, on one inherited CPU.
"""
import json
import os
from pathlib import Path
import resource
import subprocess
import tempfile
import time
from functools import lru_cache

from test_strip3 import neighbourhoods, subset_table

ROOT = Path(__file__).resolve().parents[1]


def limits():
    resource.setrlimit(resource.RLIMIT_AS, (512*1024**2,)*2)
    resource.setrlimit(resource.RLIMIT_CPU, (2,)*2)


def history_table(n, values):
    closed = neighbourhoods(n)
    result = []
    for s in range(1 << (3*n)):
        if any(s >> v & 1 and s & (closed[v] ^ (1 << v)) for v in range(3*n)):
            result.append(-1)
            continue
        remaining = (1 << (3*n))-1
        for v in range(3*n):
            if s >> v & 1:
                remaining &= ~closed[v]
        result.append(values[remaining])
    return result


def main():
    os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
    resource.setrlimit(resource.RLIMIT_AS, (512*1024**2,)*2)
    started = time.monotonic()
    base_env = {k: v for k, v in os.environ.items() if not k.startswith("STRIP3_")}

    def run(command, text=None, **settings):
        remaining = 20-(time.monotonic()-started)
        assert remaining > 0, "20 second batch budget exhausted"
        return subprocess.run(list(map(str, command)), input=text, capture_output=True, text=True,
                              timeout=min(2, remaining), preexec_fn=limits, env={**base_env, **settings})

    comparisons = certificates = 0
    with tempfile.TemporaryDirectory(prefix="strip3-demand-test-") as folder:
        folder = Path(folder)
        paths = {}
        for name, source in (("availability", "src/cfast/strip3_demand_symbolic.c"),
                             ("history", "src/cfast/strip3_history_symbolic.c"),
                             ("checker", "tests/verify_strip3_symbolic_kernel.c")):
            paths[name] = folder/name
            build = run(["gcc", "-O3", "-Wall", "-Wextra", "-Werror", "-o", paths[name], ROOT/source])
            assert build.returncode == 0, build.stderr
        for n in range(1, 6):
            values = subset_table(n)
            for name in ("availability", "history"):
                expected = values if name == "availability" else history_table(n, values)
                for mode in (("0", "1", "2") if n <= 4 else ("2",)):
                    result = run([paths[name], "table", n], STRIP3_DEMAND_GENERALISE=mode)
                    assert result.returncode == 0, result.stderr
                    assert list(map(int, result.stdout.split())) == expected, (name, n, mode)
                    comparisons += sum(g >= 0 for g in expected)
        for name in ("availability", "history"):
            for mode in ("0", "1", "2"):
                proof = folder/"kernel.txt"
                result = run([paths[name], "board", 6, "--certificate", proof], STRIP3_DEMAND_GENERALISE=mode)
                assert result.returncode == 0 and result.stdout.split() == ["3", "6", "3"], result.stderr
                verified = run([paths["checker"]], proof.read_text())
                assert verified.returncode == 0, verified.stderr
                metrics = json.loads(verified.stdout)
                assert metrics["SG"] == 3 and metrics["n"] == 6
                if mode != "0":
                    assert metrics["complete_lower_kernels"] < 3
                certificates += 1
            # Force collections while concrete query frames remain active.
            for gc, bits in (("0", "16"), ("1", "14")):
                result = run([paths[name], "board", 6], STRIP3_SYMBOLIC_BITS=bits, STRIP3_SYMBOLIC_GC=gc)
                assert result.returncode == 0 and result.stdout.split() == ["3", "6", "3"], result.stderr
                if gc == "1":
                    assert "gc=0 " not in result.stderr, result.stderr
            result = run([paths[name], "board", 6], STRIP3_SYMBOLIC_BITS="8")
            assert result.returncode == 3 and not result.stdout and "incomplete" in result.stderr
            for option in (dict(STRIP3_DEMAND_GENERALISE="3"), dict(STRIP3_DEMAND_TRACE="2"),
                           dict(STRIP3_SYMBOLIC_BITS="7"), dict(STRIP3_SYMBOLIC_GC="2")):
                result = run([paths[name], "board", 4], **option)
                assert result.returncode == 2 and not result.stdout
            for args in (("board", 0), ("board", 101), ("table", 7), ("position", 2, "40"),
                         ("position", 2, "-1"), ("position", 2, "garbage")):
                result = run([paths[name], *args])
                assert result.returncode == 2 and not result.stdout
            result = run([paths[name], "board", 4, "--certificate", folder/"absent"/"file"])
            assert result.returncode == 3 and not result.stdout and "write failed" in result.stderr
        # Check all five words, including moves across 64-bit boundaries.
        n = 100
        closed = neighbourhoods(n)

        @lru_cache(None)
        def sparse_sg(mask):
            options = {sparse_sg(mask & ~closed[v]) for v in range(3*n) if mask >> v & 1}
            g = 0
            while g in options:
                g += 1
            return g

        available = sum(1 << v for v in (0, 63, 64, 127, 128, 191, 192, 255, 256, 299))
        result = run([paths["availability"], "position", n, f"{available:x}"])
        assert result.returncode == 0 and int(result.stdout.split()[-1]) == sparse_sg(available), result.stderr
        history = sum(1 << v for v in range(3*n) if (v//3+v%3) % 2 == 0)
        for removed in ((), (0, 64, 128, 192, 256, 298)):
            s = history
            for v in removed:
                s &= ~(1 << v)
            mask = (1 << (3*n))-1
            for v in range(3*n):
                if s >> v & 1:
                    mask &= ~closed[v]
            assert mask.bit_count() <= 10
            result = run([paths["history"], "position", n, f"{s:x}"])
            assert result.returncode == 0 and int(result.stdout.split()[-1]) == sparse_sg(mask), result.stderr
        result = run([paths["history"], "position", 2, "3"])
        assert result.returncode == 2 and not result.stdout
        # History certificates have a different domain, transition and root.
        original = proof.read_text().splitlines()
        mutations = []
        changed = original[:]; changed[0] = "node-kayles-symbolic-kernel-v1"; mutations.append(changed)
        changed = original[:]; changed[2] = "1"; mutations.append(changed)
        changed = original[:]; roots = changed[3].split(); roots[0] = "0"
        changed[3] = " ".join(roots); mutations.append(changed)
        changed = original[:]; roots = changed[3].split(); roots[-1] = changed[2]
        changed[3] = " ".join(roots); mutations.append(changed)
        changed = original[:]; roots = changed[3].split(); roots[-1] = "0"
        changed[3] = " ".join(roots); mutations.append(changed)
        # Empty history singleton falsely claimed SG=0 must fail actual moves.
        changed = original[:]; width, _, count = map(int, changed[1].split()); child = 1
        for v in reversed(range(3*width)):
            # The producer normalises complement edges to an even low edge.
            low, high = child, 0
            flip = low & 1
            changed.append(f"{v} {low ^ flip} {high ^ flip}")
            child = 2*count ^ flip; count += 1
        changed[1] = f"{width} 0 {count}"; changed[3] = str(child)
        bad = run([paths["checker"]], "\n".join(changed)+"\n")
        assert bad.returncode != 0 and "outside its winning witnesses" in bad.stderr, bad
        for changed in mutations:
            bad = run([paths["checker"]], "\n".join(changed)+"\n")
            assert bad.returncode != 0 and not bad.stdout
        for text in ("\n".join(original[:-1])+"\n", "\n".join(original)+"\ntrailing\n"):
            assert run([paths["checker"]], text).returncode != 0
        assert run([paths["checker"]], proof.read_text(), STRIP3_KERNEL_BITS="8").returncode != 0
    print({"independent_SG_comparisons": comparisons, "standalone_certificates": certificates,
           "partial_lower_kernels": "verified", "semantic_tampering": "rejected",
           "word_boundaries_and_sparse_100_column_geometry": "passed",
           "full_3x100": "not solved", "wall_s": round(time.monotonic()-started, 6)})


if __name__ == "__main__":
    main()
