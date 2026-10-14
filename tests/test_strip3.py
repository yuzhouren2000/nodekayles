"""Independent checks of the three-row solver and its structural shortcuts.

Exhaust every induced subgraph of 3x5 using a direct subset recurrence (no
components, geometry normalisation, or theorem shortcuts). Exercise all 32
reduction settings, sparse positions across bit 64 and bit 125, input/budget
errors, then compare a first-move survey to the original multiword C solver.

Run: python3 tests/test_strip3.py [--survey data/strip3_survey.jsonl]
"""

import argparse
from functools import lru_cache
import json
import os
from pathlib import Path
import random
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def neighbourhoods(n):
    result = []
    for v in range(3 * n):
        c, r = divmod(v, 3)
        result.append(sum(1 << (3 * cc + rr) for rr, cc in
                          ((r, c), (r - 1, c), (r + 1, c), (r, c - 1), (r, c + 1))
                          if 0 <= rr < 3 and 0 <= cc < n))
    return result


def mex(values):
    g = 0
    while g in values:
        g += 1
    return g


def subset_table(n):
    nb = neighbourhoods(n)
    table = [0] * (1 << (3 * n))
    for mask in range(1, len(table)):
        options = set()
        bits = mask
        while bits:
            b = bits & -bits
            bits ^= b
            options.add(table[mask & ~nb[b.bit_length() - 1]])
        table[mask] = mex(options)
    return table


def sparse_value(mask):
    vertices = [v for v in range(126) if mask >> v & 1]
    closed = []
    for v in vertices:
        c, r = divmod(v, 3)
        closed.append(sum(1 << i for i, u in enumerate(vertices)
                          if abs(r - u % 3) + abs(c - u // 3) <= 1))

    @lru_cache(None)
    def value(subset):
        return mex({value(subset & ~closed[i]) for i in range(len(vertices)) if subset >> i & 1})

    return value((1 << len(vertices)) - 1)


def run(binary, args, input_text=None, flags=7, **env):
    return subprocess.run([str(binary), *args], input=input_text, text=True,
                          capture_output=True, timeout=90,
                          env={**os.environ, "STRIP3_REDUCTIONS": str(flags), **env})


def check_masks(binary, positions, expected, flags):
    result = run(binary, ["masks"], "".join(f"{n} {m:x}\n" for n, m in positions), flags)
    if result.returncode:
        raise AssertionError(result.stderr)
    got = list(map(int, result.stdout.splitlines()))
    if len(got) != len(expected):
        raise AssertionError("incomplete mask output")
    for (n, mask), a, b in zip(positions, got, expected):
        if a != b:
            raise AssertionError(f"flags={flags}, n={n}, mask={mask:x}: strip3={a}, direct={b}")


def check_original(temp, records, limit):
    # Separate adjacency layout, component code, canonicalisation and hash table.
    source = '''#define main kayles_cli_main
#include "%s"
#undef main
int main(void) {
    int n; unsigned long long hi, lo; ht_init(18);
    while (scanf("%%d %%llx %%llx", &n, &hi, &lo) == 3) {
        uint64_t position[WORDS] = {0};
        for (int r = 0; r < 3; r++) for (int c = 0; c < n; c++) {
            int v = 3*c+r, target = r*n+c;
            uint64_t word = v < 64 ? lo : hi;
            if ((word >> (v %% 64)) & 1) position[target / 64] |= 1ULL << (target %% 64);
        }
        printf("%%d\\n", xor_components(3, n, position));
    }
    return 0;
}
''' % (ROOT / "src/cfast/kayles.c")
    (temp / "original.c").write_text(source)
    subprocess.run(["gcc", "-O3", "-march=native", "-o", str(temp / "original"), str(temp / "original.c")], check=True)
    positions = {}
    for record in records:
        n = record["n"]
        if n > limit:
            continue
        full = (1 << (3 * n)) - 1
        gap = full & ~(1 << (3 * (n - 1) + 1))
        positions[n, full] = record["grundy"]
        positions[n, gap] = record["end_gap"]
        positions[n, gap & ~2] = record["both_gaps"]
        for v, closed in enumerate(neighbourhoods(n)):
            positions[n, full & ~closed] = record["options"][v]
    input_text = "".join(f"{n} {m >> 64:x} {m & ((1 << 64) - 1):x}\n" for n, m in positions)
    result = subprocess.run([str(temp / "original")], input=input_text, text=True,
                            capture_output=True, timeout=90)
    if result.returncode:
        raise AssertionError(result.stderr)
    got = list(map(int, result.stdout.splitlines()))
    if got != list(positions.values()):
        raise AssertionError("original C solver disagrees with survey")
    print(f"original C: {len(positions)} full/end-gap/both-gap/first-move values agree through n={limit}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--survey", type=Path)
    parser.add_argument("--verify-limit", type=int, default=20)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="strip3-check-") as folder:
        temp = Path(folder)
        binary = temp / "strip3"
        subprocess.run(["gcc", "-O3", "-march=native", "-Wall", "-Wextra", "-o", str(binary),
                        str(ROOT / "src/cfast/strip3.c")], check=True)
        expected = subset_table(5)
        positions = [(5, m) for m in range(len(expected))]
        for flags in range(32):
            check_masks(binary, positions, expected, flags)
        print(f"direct subset recurrence: all {len(expected)} induced 3x5 positions agree under all 32 reduction settings")

        rng = random.Random(1729)
        sparse = [(42, 0)] + [(42, 1 << v) for v in range(126)]
        for start in (18, 19, 20, 21, 22, 37):
            for _ in range(12):
                sparse.append((42, rng.randrange(1 << 12) << (3 * start)))
        for _ in range(128):
            sparse.append((42, sum(1 << v for v in rng.sample(range(126), rng.randrange(1, 13)))))
        expected_sparse = [sparse_value(mask) for _, mask in sparse]
        for flags in (0, 7, 31):
            check_masks(binary, sparse, expected_sparse, flags)
        print(f"direct sparse-graph recurrence: {len(sparse)} positions across the 64/126-bit boundaries agree")

        for command, text in ((["board", "43"], None), (["board", ""], None),
                              (["masks"], "42 40000000000000000000000000000000\n"),
                              (["masks"], "2 100\n"), (["masks"], "3 xyz\n"),
                              (["masks"], "bad\n"), (["masks"], "3 1 extra\n")):
            assert run(binary, command, text).returncode == 2, (command, text)
        for variable, value in (("STRIP3_HT_BITS", "7"), ("STRIP3_HT_BITS", ""),
                                ("STRIP3_REDUCTIONS", "32")):
            assert run(binary, ["board", "3"], **{variable: value}).returncode == 2
        budget = run(binary, ["strip", "16"], STRIP3_HT_BITS="8")
        assert budget.returncode == 3 and "memo budget reached" in budget.stderr
        completed = [list(map(int, line.split())) for line in budget.stdout.splitlines()]
        assert completed and completed[-1][1] < 16
        print("invalid inputs rejected; bounded memo exits with code 3 and preserves only completed rows")

        if args.survey:
            records = [json.loads(line) for line in args.survey.read_text().splitlines()]
        else:
            survey = run(binary, ["survey", str(args.verify_limit)])
            assert survey.returncode == 0, survey.stderr
            records = [json.loads(line) for line in survey.stdout.splitlines()]
        check_original(temp, records, args.verify_limit)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
