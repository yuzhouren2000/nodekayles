"""Independent SG, arbitrary-target, sum and budget checks of search prototypes.

Finite correctness checks do not establish scalable complexity. No experiment
in this file computes a board larger than 3x18 (made of small separated pieces).
"""
import os
from pathlib import Path
import random
import subprocess
import tempfile

from test_strip3 import subset_table, sparse_value

ROOT = Path(__file__).resolve().parents[1]


def main():
    five, four = subset_table(5), subset_table(4)
    rng = random.Random(10030)
    with tempfile.TemporaryDirectory(prefix="strip3-search-experiments-") as folder:
        symbolic, proof = [Path(folder) / name for name in ("symbolic", "dfpn")]
        for binary, source in ((symbolic, "strip3_symbolic.c"), (proof, "strip3_dfpn.c")):
            subprocess.run(["gcc", "-O3", "-Wall", "-Wextra", "-pthread", "-o", str(binary),
                            str(ROOT / "src/cfast" / source)], check=True, timeout=10)
        result = subprocess.run([str(symbolic), "table", "5"], capture_output=True, text=True, timeout=3)
        assert result.returncode == 0, result.stderr
        assert list(map(int, result.stdout.split())) == five
        print("symbolic fixed point: all 32768 induced 3x5 SG values agree with subset mex")

        positions = list(range(len(five))); rng.shuffle(positions)
        sum_positions = []
        for count in (2, 3):
            for _ in range(64):
                masks = [rng.randrange(1 << 12) for _ in range(count)]
                g, mask = 0, 0
                for i, m in enumerate(masks):
                    mask |= m << (18*i); g ^= four[m]
                sum_positions.append((18, mask, g))
        sparse = [(42, 1 << 125), (42, (15 << 61))]
        sparse += [(42, sum(1 << v for v in rng.sample(range(126), 12))) for _ in range(32)]
        sum_positions += [(n, m, sparse_value(m)) for n, m in sparse]
        for threshold in (0, 10):
            env = {**os.environ, "STRIP3_TARGET_BITS": "18", "STRIP3_DFPN_THRESHOLD": str(threshold)}
            result = subprocess.run([str(proof), "masks"], input="".join(f"5 {m:x}\n" for m in positions),
                                    capture_output=True, text=True, timeout=3, env=env)
            assert result.returncode == 0, result.stderr
            assert list(map(int, result.stdout.split())) == [five[m] for m in positions]
            queries = [(5, m, k, five[m]) for m in rng.sample(positions, 2048) for k in range(12)]
            queries += [(n, m, k, g) for n, m, g in sum_positions for k in range(16)]
            rng.shuffle(queries)
            result = subprocess.run([str(proof), "queries"], input="".join(f"{n} {m:x} {k}\n" for n,m,k,_ in queries),
                                    capture_output=True, text=True, timeout=3, env=env)
            assert result.returncode == 0, result.stderr
            assert list(map(int, result.stdout.split())) == [int(g == k) for _,_,k,g in queries]
        print("DFPN: all 3x5 SGs in shuffled order; independent arbitrary-target, sum and cross-word queries agree")
        for binary, setting, n in ((symbolic, "STRIP3_SYMBOLIC_BITS", 8), (proof, "STRIP3_DFPN_BITS", 12)):
            env = {**os.environ, setting: "8"}
            result = subprocess.run([str(binary), "board", str(n)], capture_output=True, text=True, timeout=3, env=env)
            assert result.returncode == 3 and not result.stdout and "incomplete" in result.stderr, result
        print("symbolic/DFPN budgets terminate without any false completed SG output")


if __name__ == "__main__":
    main()
