"""Check exact history-image language, symbolic SG and BDD compaction.

The reference enumerates independent vertex sets and uses direct subset mex;
it shares neither the column transducer nor any BDD operation with the solver.
"""
import os
from pathlib import Path
import resource
import subprocess
import tempfile

from test_strip3 import neighbourhoods, subset_table

ROOT = Path(__file__).resolve().parents[1]


def limits():
    resource.setrlimit(resource.RLIMIT_AS, (512*1024**2,)*2)
    resource.setrlimit(resource.RLIMIT_CPU, (2,)*2)


def run(binary, mode, n, **settings):
    env = {k: v for k, v in os.environ.items() if not k.startswith("STRIP3_")}
    env.update(settings)
    return subprocess.run([str(binary), mode, str(n)], capture_output=True, text=True,
                          timeout=2, env=env, preexec_fn=limits)


def independent_history_images(n):
    closed = neighbourhoods(n)
    full = (1 << (3*n))-1
    images = set()

    def visit(candidates, available):
        images.add(available)
        while candidates:
            bit = candidates & -candidates
            candidates ^= bit
            v = bit.bit_length()-1
            visit(candidates & ~closed[v], available & ~closed[v])

    visit(full, full)
    return images


def main():
    os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
    resource.setrlimit(resource.RLIMIT_AS, (512*1024**2,)*2)
    expected_counts = (4, 14, 54, 172, 600, 2037)
    checked_values, checked_moves = 0, 0
    with tempfile.TemporaryDirectory(prefix="strip3-reachable-symbolic-test-") as folder:
        binary, baseline = [Path(folder)/name for name in ("reachable", "induced")]
        for path, source in ((binary, "strip3_reachable_symbolic.c"), (baseline, "strip3_symbolic.c")):
            subprocess.run(["gcc", "-O3", "-Wall", "-Wextra", "-Werror", "-o", str(path),
                            str(ROOT/"src/cfast"/source)], check=True, timeout=2, preexec_fn=limits)
        for n in range(1, 7):
            images = independent_history_images(n)
            assert len(images) == expected_counts[n-1]
            result = run(binary, "domain", n)
            assert result.returncode == 0, result.stderr
            assert list(map(int, result.stdout.split())) == [int(m in images) for m in range(1 << (3*n))]
            # Every real successor remains in the independent-history image.
            closed = neighbourhoods(n)
            for mask in images:
                for v in range(3*n):
                    if mask >> v & 1:
                        assert mask & ~closed[v] in images
                        checked_moves += 1
            expected = subset_table(n)
            for gc, bits in (("0", "20"), ("1", "20"), ("1", "18")):
                result = run(binary, "table", n, STRIP3_SYMBOLIC_GC=gc, STRIP3_SYMBOLIC_BITS=bits)
                assert result.returncode == 0, result.stderr
                assert list(map(int, result.stdout.split())) == [g if m in images else -1
                                                                 for m, g in enumerate(expected)]
                checked_values += len(images)
                if n == 6 and gc == "1" and bits == "18":
                    assert "gc=0 " not in result.stderr, result.stderr
            if n == 5:
                result = run(baseline, "table", n)
                assert result.returncode == 0, result.stderr
                assert list(map(int, result.stdout.split())) == expected
        for n in (0, 101, "6bad"):
            assert run(binary, "board", n).returncode == 2
        for mode in ("domain", "table"):
            assert run(binary, mode, 7).returncode == 2
        for setting in (dict(STRIP3_SYMBOLIC_BITS="7"), dict(STRIP3_SYMBOLIC_GC="2")):
            result = run(binary, "board", 4, **setting)
            assert result.returncode == 2 and not result.stdout
        for gc in ("0", "1"):
            result = run(binary, "board", 12, STRIP3_SYMBOLIC_BITS="8", STRIP3_SYMBOLIC_GC=gc)
            assert result.returncode == 3 and not result.stdout and "incomplete" in result.stderr
        assert run(baseline, "domain", 2).returncode == 2
    print({"domain_counts": expected_counts, "independent_SG_comparisons": checked_values,
           "real_move_closure_checks": checked_moves, "baseline_induced_SG_comparisons": 32768,
           "forced_compaction_and_budget_failures": "passed"})


if __name__ == "__main__":
    main()
