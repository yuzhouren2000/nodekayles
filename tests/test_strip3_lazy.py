"""Independent SG and shuffled-query checks for suspended three-valued search."""
import os
from pathlib import Path
import random
import resource
import subprocess
import tempfile
import json

from test_strip3 import subset_table, sparse_value
from verify_strip3_equality_proof import verify

ROOT = Path(__file__).resolve().parents[1]


def limits():
    resource.setrlimit(resource.RLIMIT_AS, (512*1024**2,)*2)
    resource.setrlimit(resource.RLIMIT_CPU, (3,)*2)


def run(binary, mode, text=None, *args, **extra):
    env = {k: v for k, v in os.environ.items() if not k.startswith("STRIP3_")}
    env.update(STRIP3_TARGET_BITS="18", STRIP3_LAZY_BITS="16")
    env.update(extra)
    return subprocess.run([str(binary), mode, *map(str, args)], input=text, text=True, capture_output=True,
                          env=env, timeout=2, preexec_fn=limits)


def main():
    os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
    expected = subset_table(5); rng = random.Random(3400)
    masks = list(range(len(expected))); rng.shuffle(masks)
    queries = []
    for m in masks:
        targets = list(range(10))+[127]; rng.shuffle(targets)
        queries.extend((m, k) for k in targets)
    with tempfile.TemporaryDirectory(prefix="strip3-lazy-test-") as folder:
        binary = Path(folder) / "lazy"
        subprocess.run(["gcc", "-O3", "-Wall", "-Wextra", "-pthread", "-o", str(binary),
                        str(ROOT / "src/cfast/strip3_lazy.c")], check=True, timeout=2)
        configs = [(1, 1, 2, 4, 0), (16, 2, 2, 4, 0), (256, 2, 2, 4, 0),
                   (16, 1, 3, 4, 0), (16, 2, 0, 128, 0), (4096, 2, 2, 4, 0), (256, 2, 2, 4, 1)]
        for quantum, growth, order, frontier, cache in configs:
            env = dict(STRIP3_LAZY_QUANTUM=str(quantum), STRIP3_LAZY_GROWTH_BITS=str(growth),
                       STRIP3_TARGET_ORDER=str(order), STRIP3_TARGET_FRONTIER=str(frontier),
                       STRIP3_TARGET_CACHE_FIRST=str(cache))
            result = run(binary, "queries", "".join(f"5 {m:x} {k}\n" for m, k in queries), **env)
            assert result.returncode == 0, result.stderr
            assert list(map(int, result.stdout.split())) == [int(expected[m] == k) for m, k in queries], env
            result = run(binary, "masks", "".join(f"5 {m:x}\n" for m in masks), **env)
            assert result.returncode == 0, result.stderr
            assert list(map(int, result.stdout.split())) == [expected[m] for m in masks], env
        sparse = [0, 1 << 125, 7 << 120, 15 << 61]
        for start in (18, 20, 21, 37, 38):
            sparse.extend(rng.randrange(1 << 12) << (3*start) for _ in range(6))
        sparse.extend(sum(1 << v for v in rng.sample(range(126), 12)) for _ in range(32))
        result = run(binary, "masks", "".join(f"42 {m:x}\n" for m in sparse))
        assert result.returncode == 0, result.stderr
        assert list(map(int, result.stdout.split())) == [sparse_value(m) for m in sparse]
        for setting in (dict(STRIP3_TARGET_BITS="8"), dict(STRIP3_LAZY_BITS="8"),
                        dict(STRIP3_LAZY_WORK_LIMIT="10"), dict(STRIP3_TARGET_MOVE_LIMIT="10")):
            result = run(binary, "board", None, 16, **setting)
            assert result.returncode == 3 and not result.stdout and "incomplete" in result.stderr, result
        for setting in (dict(STRIP3_LAZY_QUANTUM="0"), dict(STRIP3_LAZY_GROWTH_BITS="0"), dict(STRIP3_LAZY_WORK_LIMIT="-1")):
            assert run(binary, "board", None, 4, **setting).returncode == 2
        # A separate certificate checker must reconstruct every logical step
        # from completed facts; suspended obligations cannot be proof leaves.
        proof_binary = Path(folder) / "lazy_proof"
        subprocess.run(["gcc", "-O3", "-Wall", "-Wextra", "-pthread", "-o", str(proof_binary),
                        str(ROOT / "src/cfast/strip3_lazy_proof.c")], check=True, timeout=2)
        cert = Path(folder) / "proof.json"
        result = run(proof_binary, "board", None, 12, "--certificate", str(cert),
                     STRIP3_TARGET_CACHE_FIRST="1", STRIP3_LAZY_QUANTUM="256", STRIP3_PROOF_BITS="18")
        assert result.returncode == 0, result.stderr
        metrics = json.loads(result.stdout.splitlines()[-1])
        checked = verify(json.loads(cert.read_text()))
        assert checked["SG"] == 5 and metrics["oracle_calls"] == 0
        for key in ("n", "SG", "tree_cost", "sum_rules", "edges"): assert checked[key] == metrics[key]
    print(f"{len(configs)} configurations agree on all 32768 induced 3x5 SG values and {len(configs)*len(queries)} shuffled equality queries")
    print(f"{len(sparse)} cross-bit sparse graphs agree; suspended and exhausted work never become false SG answers")
    print("a complete n=12 proof from resumed search verifies without an exact-SG oracle")


if __name__ == "__main__":
    main()
