"""Independent subset-mex checks for absent-option upper-bound search."""
import os
from pathlib import Path
import random
import resource
import subprocess
import tempfile

from test_strip3 import subset_table, sparse_value

ROOT = Path(__file__).resolve().parents[1]


def limits():
    resource.setrlimit(resource.RLIMIT_AS, (512*1024**2,)*2)
    resource.setrlimit(resource.RLIMIT_CPU, (3,)*2)


def run(binary, mode, text=None, *args, **extra):
    env = {k: v for k, v in os.environ.items() if not k.startswith("STRIP3_")}
    env["STRIP3_TARGET_BITS"] = "18"
    env.update(extra)
    return subprocess.run([str(binary), mode, *map(str, args)], input=text, text=True,
                          capture_output=True, timeout=2, env=env, preexec_fn=limits)


def main():
    os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
    expected = subset_table(5)
    rng = random.Random(3100)
    masks = list(range(len(expected))); rng.shuffle(masks)
    queries = []
    for mask in masks:
        targets = list(range(10)) + [127]; rng.shuffle(targets)
        queries.extend((mask, target) for target in targets)
    text = "".join(f"5 {mask:x} {target}\n" for mask, target in queries)
    query_expected = [int(expected[m] == k) for m, k in queries]
    with tempfile.TemporaryDirectory(prefix="strip3-absence-test-") as folder:
        binary = Path(folder) / "absence"
        subprocess.run(["gcc", "-O3", "-Wall", "-Wextra", "-pthread", "-o", str(binary),
                        str(ROOT / "src/cfast/strip3_absence.c")], check=True, timeout=2)
        configs = [dict(STRIP3_ABSENCE_PROBE=str(probe)) for probe in range(3)]
        configs += [dict(STRIP3_ABSENCE_PROBE="1", STRIP3_TARGET_CACHE_FIRST="1"),
                    dict(STRIP3_ABSENCE_PROBE="2", STRIP3_TARGET_ORDER="3", STRIP3_TARGET_FRONTIER="3"),
                    dict(STRIP3_ABSENCE_PROBE="0", STRIP3_TARGET_CACHE_FIRST="1"),
                    dict(STRIP3_ABSENCE_PROBE="2", STRIP3_TARGET_CACHE_FIRST="1")]
        for config in configs:
            result = run(binary, "queries", text, **config)
            assert result.returncode == 0, result.stderr
            assert list(map(int, result.stdout.split())) == query_expected, config
            result = run(binary, "masks", "".join(f"5 {m:x}\n" for m in masks), **config)
            assert result.returncode == 0, result.stderr
            assert list(map(int, result.stdout.split())) == [expected[m] for m in masks], config
        print(f"{len(configs)} configurations: all 32768 induced 3x5 SG values and {len(queries)*len(configs)} shuffled equality queries agree")
        sparse = [0, 1 << 125, 7 << 120, 15 << 61]
        for start in (18, 20, 21, 37, 38):
            sparse.extend(rng.randrange(1 << 12) << (3*start) for _ in range(6))
        sparse.extend(sum(1 << v for v in rng.sample(range(126), 12)) for _ in range(32))
        result = run(binary, "masks", "".join(f"42 {m:x}\n" for m in sparse))
        assert result.returncode == 0, result.stderr
        assert list(map(int, result.stdout.split())) == [sparse_value(m) for m in sparse]
        result = run(binary, "board", None, 16, STRIP3_TARGET_BITS="8")
        assert result.returncode == 3 and not result.stdout and "incomplete" in result.stderr
        result = run(binary, "board", None, 16, STRIP3_TARGET_MOVE_LIMIT="10")
        assert result.returncode == 3 and not result.stdout and "moves=10" in result.stderr
        for probe in ("-1", "3", "", "bad"):
            assert run(binary, "board", None, 4, STRIP3_ABSENCE_PROBE=probe).returncode == 2
        print(f"{len(sparse)} cross-bit sparse positions agree; invalid modes and exhausted budgets produce no SG")


if __name__ == "__main__":
    main()
