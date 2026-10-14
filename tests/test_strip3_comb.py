#!/usr/bin/env python3
"""Check complete comb recognition, reflected embeddings and SG preservation."""

import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from export_strip3_comb import MAXCOLS


def arm(n):
    return tuple(4 if i % 2 == 0 else 6 for i in range(n))


def canonical(columns):
    cells = {(r, c) for c, bits in enumerate(columns) for r in range(3) if bits >> r & 1}
    if not cells:
        return 0
    r0, c0 = min(r for r, c in cells), min(c for r, c in cells)
    cells = {(r - r0, c - c0) for r, c in cells}
    h, n = max(r for r, c in cells) + 1, max(c for r, c in cells) + 1
    shapes = [cells]
    if h == n == 3:
        shapes.append({(c, r) for r, c in cells})
    return min(sum(1 << (3 * (n - 1 - c if horizontal else c)
                              + (h - 1 - r if vertical else r)) for r, c in shape)
               for shape in shapes for horizontal in (False, True) for vertical in (False, True))


def position(columns):
    return len(columns), sum(bits << (3 * c) for c, bits in enumerate(columns))


def families():
    for width in range(2, 6):
        for left in range(MAXCOLS - width - 1):
            for right in range(MAXCOLS - width - 1 - left):
                yield arm(left)[::-1] + (6,) + (7,) * width + (6,) + arm(right), 8, 0
    for bridge in range(1, MAXCOLS - 8, 2):
        for left in range(MAXCOLS - 8 - bridge):
            for right in range(MAXCOLS - 8 - bridge - left):
                columns = (arm(left)[::-1] + (6, 7, 7, 6) + arm(bridge)
                           + (6, 7, 7, 7, 6) + arm(right))
                yield columns, 16, int(left >= 5 and right >= 5 and bridge >= 11)


def run(binary, positions, flags, timeout=60):
    return subprocess.run([str(binary), "masks"], text=True, capture_output=True, timeout=timeout,
                          input="".join(f"{n} {mask:x}\n" for n, mask in positions),
                          env={**os.environ, "STRIP3_REDUCTIONS": str(flags),
                               "STRIP3_THREADS": "1", "STRIP3_HT_BITS": "24"})


def main():
    examples = list(families())
    catalog = {canonical(cols): (kind, formula) for cols, kind, formula in examples}
    positions, expected = [], []
    for cols, kind, formula in examples:
        reflected = tuple(((bits & 1) << 2) | (bits & 2) | ((bits & 4) >> 2) for bits in cols)
        for word in (cols, cols[::-1], reflected, reflected[::-1]):
            positions.append(position(word))
            expected.append((kind, formula))
    # Translate selected families to the last column, including bit 125.
    for cols, kind, formula in examples[::41]:
        positions.append(position((0,) * (MAXCOLS - len(cols)) + cols))
        expected.append((kind, formula))
    rng = random.Random(20261014)
    modified = []
    for cols, _, _ in rng.sample(examples, 1000):
        n, mask = position(cols)
        v = rng.randrange(3 * n)
        changed = mask ^ (1 << v)
        modified.append((n, changed))
        positions.append((n, changed))
        word = tuple((changed >> (3 * c)) & 7 for c in range(n))
        expected.append(catalog.get(canonical(word), (0, 0)))
    with tempfile.TemporaryDirectory(prefix="strip3-comb-check-") as folder:
        temp = Path(folder)
        binary = temp / "strip3"
        subprocess.run(["gcc", "-O3", "-march=native", "-Wall", "-Wextra", "-pthread",
                        "-o", str(binary), str(ROOT / "src/cfast/strip3.c")], check=True)
        # Probe the recognizer directly, including cases handled by old
        # symmetry reductions before the recognizer in the normal recursion.
        wrapper = temp / "probe.c"
        wrapper.write_text(f'''#define main strip3_cli_main
#include "{ROOT / 'src/cfast/strip3.c'}"
#undef main
int main(void) {{
    initialise(); reductions = 31;
    int n; unsigned long long hi, lo;
    while (scanf("%d %llx %llx", &n, &hi, &lo) == 3) {{
        Mask m = ((Mask)hi << 64) | lo;
        int g = -1, formula = 0, kind = 0;
        if (m) {{ m = canonical(m); kind = certified_comb(m, width(m), &g, &formula); }}
        printf("%d %d %d\\n", kind, g, formula);
    }}
    return 0;
}}
''')
        probe = temp / "probe"
        subprocess.run(["gcc", "-O3", "-march=native", "-Wall", "-Wextra", "-pthread",
                        "-o", str(probe), str(wrapper)], check=True)
        result = subprocess.run([str(probe)], text=True, capture_output=True, timeout=60,
                                input="".join(f"{n} {m >> 64:x} {m & ((1 << 64) - 1):x}\n" for n, m in positions),
                                env={**os.environ, "STRIP3_THREADS": "1", "STRIP3_HT_BITS": "24"})
        assert result.returncode == 0, result.stderr
        recognized = [tuple(map(int, line.split())) for line in result.stdout.splitlines()]
        assert len(recognized) == len(expected)
        for pos, want, got in zip(positions, expected, recognized):
            assert (got[0], got[2]) == want, (pos, want, got)
        # All bounded parameter pairs/triples, plus one-cell changes, are
        # independently computed by the old C recursion with these rules OFF.
        originals = [position(cols) for cols, _, _ in examples]
        off = run(binary, originals + modified, 7)
        assert off.returncode == 0, off.stderr
        sg = list(map(int, off.stdout.splitlines()))
        assert len(sg) == len(originals) + len(modified)
        for i, value in enumerate(sg[:len(originals)]):
            assert all(recognized[4 * i + j][1] == value for j in range(4)), (i, value)
        for flags in (15, 23, 31):
            on = run(binary, originals + modified, flags)
            assert on.returncode == 0 and list(map(int, on.stdout.splitlines())) == sg, on.stderr
        # The sharded memo also uses immutable comb tables.
        threaded = run(binary, originals[:1000] + modified[:100], 31)
        assert threaded.returncode == 0, threaded.stderr
        parallel = subprocess.run([str(binary), "masks"], text=True, capture_output=True, timeout=60,
                                  input="".join(f"{n} {m:x}\n" for n, m in originals[:1000] + modified[:100]),
                                  env={**os.environ, "STRIP3_THREADS": "2", "STRIP3_HT_PREALLOC": "1",
                                       "STRIP3_HT_BITS": "24", "STRIP3_REDUCTIONS": "31"})
        assert parallel.returncode == 0 and parallel.stdout == threaded.stdout, parallel.stderr
    print(f"{len(examples)} full-core parameter cases; {len(positions)} reflected/translated/modified recognition checks; "
          "SG agrees with recursive C under 7/15/23/31; 2-thread output agrees")


if __name__ == "__main__":
    main()
