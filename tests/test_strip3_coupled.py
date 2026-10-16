"""Check unresolved sums directly against independent component subset mex."""
from pathlib import Path
import os
import random
import re
import subprocess
import tempfile

from test_strip3 import subset_table

ROOT = Path(__file__).resolve().parents[1]


def run(binary, mode, text, threshold, nim_first):
    result = subprocess.run([str(binary), mode], input=text, text=True,
                            capture_output=True, timeout=10,
                            env={**os.environ, "STRIP3_TARGET_BITS": "18",
                                 "STRIP3_TARGET_ORDER": "2", "STRIP3_TARGET_LOW_FIRST": "0",
                                 "STRIP3_COUPLED_THRESHOLD": str(threshold),
                                 "STRIP3_COUPLED_NIM_FIRST": str(nim_first)})
    assert result.returncode == 0, result.stderr
    return list(map(int, result.stdout.split())), result.stderr


def main():
    small = subset_table(4)
    full = subset_table(5)
    rng = random.Random(10002)
    positions = []
    for count in (2, 3):
        for _ in range(96):
            masks = [rng.randrange(1 << 12) for _ in range(count)]
            m = sum(mask << (18*i) for i, mask in enumerate(masks))
            g = 0
            for mask in masks:
                g ^= small[mask]
            positions.append((m, g))
    with tempfile.TemporaryDirectory(prefix="strip3-coupled-test-") as folder:
        binary = Path(folder) / "coupled"
        subprocess.run(["gcc", "-O3", "-Wall", "-Wextra", "-pthread", "-o", str(binary),
                        str(ROOT / "src/cfast/strip3_coupled.c")], check=True, timeout=10)
        for threshold, nim_first in ((0, 0), (0, 1), (4, 0), (12, 0)):
            # Fresh process: component SGs have not already populated its memo.
            queries = [(m, k, g) for m, g in positions for k in range(16)]
            rng.shuffle(queries)
            got, stats = run(binary, "queries", "".join(f"18 {m:x} {k}\n" for m, k, _ in queries),
                             threshold, nim_first)
            assert got == [int(g == k) for _, k, g in queries]
            if not threshold:
                assert int(re.search(r"sum_memo=(\d+)", stats)[1]) > 0, stats
            got, _ = run(binary, "masks", "".join(f"5 {m:x}\n" for m in range(len(full))), threshold, nim_first)
            assert got == full
        print("4 sum strategies: all 32768 induced 3x5 SGs and 12288 two/three-sum target queries agree")

        # A tuple contains repeated components: cancellation must cross gaps.
        mask = 0xfff | (0xfff << 18) | (0xfff << 36)
        got, _ = run(binary, "masks", f"18 {mask:x}\n", 0, 0)
        assert got == [small[0xfff]]
        print("equal components cancel across independent spatial embeddings")


if __name__ == "__main__":
    main()
