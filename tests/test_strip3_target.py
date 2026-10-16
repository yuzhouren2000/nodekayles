"""Independent subset-mex and out-of-order equality checks for the prototype."""
import os
from pathlib import Path
import random
import subprocess
import tempfile

from test_strip3 import subset_table, sparse_value

ROOT = Path(__file__).resolve().parents[1]


def run(binary, mode, text=None, order=1, *args, **env):
    return subprocess.run([str(binary), mode, *map(str, args)], input=text,
                          text=True, capture_output=True, timeout=10,
                          env={**os.environ, "STRIP3_TARGET_BITS": "18",
                               "STRIP3_TARGET_ORDER": str(order), **env})


def main():
    expected = subset_table(5)  # no components, symmetry, or target recursion
    rng = random.Random(100)
    with tempfile.TemporaryDirectory(prefix="strip3-target-test-") as folder:
        binary = Path(folder) / "target"
        subprocess.run(["gcc", "-O3", "-Wall", "-Wextra", "-pthread", "-o", str(binary),
                        str(ROOT / "src/cfast/strip3_target.c")], check=True, timeout=10)
        masks = list(range(len(expected)))
        rng.shuffle(masks)
        for order, low_first in ((0, 0), (1, 0), (2, 0), (2, 1)):
            # Query targets before asking for the exact SG of each position.
            queries = []
            for m in masks:
                targets = list(range(9)) + [127]
                rng.shuffle(targets)
                queries.extend((m, k) for k in targets)
            env = {"STRIP3_TARGET_LOW_FIRST": str(low_first)}
            result = run(binary, "queries", "".join(f"5 {m:x} {k}\n" for m, k in queries), order, **env)
            assert result.returncode == 0, result.stderr
            assert list(map(int, result.stdout.split())) == [int(expected[m] == k) for m, k in queries]
            result = run(binary, "masks", "".join(f"5 {m:x}\n" for m in masks), order, **env)
            assert result.returncode == 0, result.stderr
            assert list(map(int, result.stdout.split())) == [expected[m] for m in masks]
        print("all 32768 induced 3x5 positions: exact SG and 1310720 out-of-order equality queries agree")

        sparse = [0, (1 << 125), (7 << 120), (15 << 61)]
        for start in (18, 20, 21, 37, 38):
            sparse += [rng.randrange(1 << 12) << (3*start) for _ in range(6)]
        sparse += [sum(1 << v for v in rng.sample(range(126), 12)) for _ in range(32)]
        result = run(binary, "masks", "".join(f"42 {m:x}\n" for m in sparse))
        assert result.returncode == 0, result.stderr
        assert list(map(int, result.stdout.split())) == [sparse_value(m) for m in sparse]
        print(f"{len(sparse)} sparse positions across bit 64/126 agree with independent graph mex")

        for mode, text, args in (("board", None, (43,)), ("equal", None, (3, 128)),
                                 ("queries", "3 3ff 0\n", ()), ("masks", "3 xyz\n", ()),
                                 ("queries", "3 7 0 extra\n", ())):
            assert run(binary, mode, text, 1, *args).returncode == 2
        result = run(binary, "board", None, 1, 16, STRIP3_TARGET_BITS="8")
        assert result.returncode == 3 and not result.stdout and "incomplete" in result.stderr
        print("invalid input rejected; exhausted memo emits no SG")


if __name__ == "__main__":
    main()
