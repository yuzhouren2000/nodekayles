"""Independently check finite connection ranks, witnesses, and reachable joins.

Column elimination below does not call the experiment's rank routine. Direct
subset mex checks a complete small connection matrix; a separate, row-major
multiword solver replays selected entries from the larger saved matrices.
Neither finite ranks nor entry sampling prove an asymptotic lower bound.
"""
from functools import cache
import gzip
from itertools import product
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from strip3_connection_rank import catalogue, gf2_witness, measure, ranks
from test_strip3 import sparse_value, subset_table


def column_rank(rows, ncols):
    """Build column vectors, then reduce at LOW row pivots."""
    columns = [0] * ncols
    for i, row in enumerate(rows):
        while row:
            low = row & -row
            columns[low.bit_length()-1] |= 1 << i
            row ^= low
    pivots = [0] * len(rows)
    rank = 0
    for column in columns:
        while column:
            low = column & -column
            i = low.bit_length()-1
            if not pivots[i]:
                pivots[i] = column
                rank += 1
                break
            column ^= pivots[i]
    return rank


def bit_rows(matrix, name):
    ncols = len(matrix[0])
    nbits = max(1, max(map(max, matrix)).bit_length())
    if name == "SG_shared_features":
        return [sum(((g >> b) & 1) << (b*ncols+j)
                    for b in range(nbits) for j, g in enumerate(row)) for row in matrix], ncols*nbits
    if name == "P":
        return [sum((g == 0) << j for j, g in enumerate(row)) for row in matrix], ncols
    bit = int(name.removeprefix("SG_bit_"))
    return [sum(((g >> bit) & 1) << j for j, g in enumerate(row)) for row in matrix], ncols


def verify_matrix(sample):
    matrix = [bytes.fromhex(row) for row in sample["matrix_hex"]]
    assert len(matrix) == sample["objects"]
    assert all(len(row) == sample["contexts"] for row in matrix)
    assert len(set(matrix)) == sample["distinct_behaviours"]
    assert max(map(max, matrix)) == sample["maximum_SG"]
    for name, witness in sample["ranks"].items():
        rows, ncols = bit_rows(matrix, name)
        r, selected, columns = witness["rank"], witness["rows"], witness["columns"]
        assert len(set(selected)) == len(set(columns)) == r
        assert all(0 <= i < len(rows) for i in selected)
        assert all(0 <= j < ncols for j in columns)
        assert column_rank(rows, ncols) == r, (sample["width"], name)
        minor = [sum(((rows[i] >> j) & 1) << k for k, j in enumerate(columns)) for i in selected]
        assert column_rank(minor, r) == r, (sample["width"], name, "singular witness")
    return matrix


def remaining(history, width):
    """Coordinate deletion, without the generator's bit neighbourhoods."""
    assert all(abs(c-d)+abs(r-s) > 1 for c, r in history for d, s in history if (c, r) != (d, s))
    points = {(c, r) for c in range(width) for r in range(3)}
    for c, r in history:
        points.difference_update(((c, r), (c-1, r), (c+1, r), (c, r-1), (c, r+1)))
    # All guarded joins must be one component, spanning both blocks.
    seen = {next(iter(points))}; todo = list(seen)
    while todo:
        c, r = todo.pop()
        new = {(c-1, r), (c+1, r), (c, r-1), (c, r+1)} & (points-seen)
        seen.update(new); todo.extend(new)
    assert seen == points
    return sum(1 << (3*c+r) for c, r in points)


@cache
def guarded_histories(width):
    found = {}
    # Exhaust candidates independently, rather than importing words/regular.
    for choices in product((-1, 0, 2), repeat=width-4):
        history = tuple((i+2, row) for i, row in enumerate(choices) if row >= 0)
        if any(d-c < (2 if r == s else 3) for (c, r), (d, s) in zip(history, history[1:])):
            continue
        mask = remaining(history, width)
        assert mask not in found
        found[mask] = history
    return found


def verify_guarded(sample):
    width = sample["width"]
    objects = guarded_histories(width)
    assert sorted(objects) == [int(s, 16) for s in sample["objects_hex"]]
    for left, history in objects.items():
        for right_hex, size in sample["contexts_hex"]:
            right = int(right_hex, 16)
            joined = history + tuple((c+width, r) for c, r in guarded_histories(size)[right])
            assert remaining(joined, width+size) == left | (right << (3*width))


def original_samples(directory, positions):
    # Separate geometry, component recursion, canonicalisation and memo layout.
    source = r'''#define main kayles_cli_main
#include "%s"
#undef main
int main(void) {
    int n; unsigned long long hi, lo; ht_init(18);
    while (scanf("%%d %%llx %%llx", &n, &hi, &lo) == 3) {
        uint64_t position[WORDS] = {0};
        for (int r = 0; r < 3; r++) for (int c = 0; c < n; c++) {
            int v = 3*c+r, target = r*n+c;
            uint64_t word = v < 64 ? lo : hi;
            if ((word >> (v %% 64)) & 1) position[target/64] |= 1ULL << (target%%64);
        }
        printf("%%d\n", xor_components(3, n, position));
    }
    return 0;
}
''' % (ROOT / "src/cfast/kayles.c")
    driver = directory / "original.c"; driver.write_text(source)
    binary = directory / "original"
    subprocess.run(["gcc", "-O3", "-Wall", "-Wextra", "-o", str(binary), str(driver)], check=True, timeout=2)
    stdin = "".join(f"{n} {mask >> 64:x} {mask & ((1 << 64)-1):x}\n" for n, mask in positions)
    # Same hard limits as the experiment. A timeout fails verification.
    from strip3_connection_rank import resource_limits
    result = subprocess.run([str(binary)], input=stdin, text=True, capture_output=True, timeout=2,
                            env={**os.environ, "KAYLES_HT_MAXBITS": "20"},
                            preexec_fn=lambda: resource_limits(2))
    assert result.returncode == 0, result.stderr
    assert list(map(int, result.stdout.split())) == list(positions.values())


def main():
    os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
    rng = random.Random(100)
    for height, width in ((0, 0), (7, 4), (4, 9), (20, 24)):
        for _ in range(12):
            rows = [rng.getrandbits(width) for _ in range(height)]
            witness = gf2_witness(rows)
            assert witness["rank"] == column_rank(rows, width)
    reports = []
    for path in sys.argv[1:]:
        data = gzip.decompress(Path(path).read_bytes()) if path.endswith(".gz") else Path(path).read_bytes()
        reports.append(json.loads(data))
    positions = {}
    small_checks = 0
    matrices = 0
    with tempfile.TemporaryDirectory(prefix="strip3-rank-verify-") as folder:
        directory = Path(folder)
        # A full matrix of all reachable width 1..3 objects and width <=2
        # contexts: all joins fit in 3x5, covered by the direct subset table.
        binary = directory / "strip3"
        subprocess.run(["gcc", "-O3", "-pthread", "-o", str(binary), str(ROOT / "src/cfast/strip3.c")],
                       check=True, timeout=2)
        expected = subset_table(5)
        contexts = [(m, n) for n in (1, 2) for m in catalogue(n)]
        for width in (1, 2, 3):
            objects = catalogue(width)
            info, matrix = measure(binary, objects, contexts, width, 2, directory)
            assert info["status"] == "complete", info
            for i, left in enumerate(objects):
                for j, (right, size) in enumerate(contexts):
                    mask = left | (right << (3*width))
                    assert matrix[i][j] == expected[mask]
                    small_checks += 1
            sample = dict(objects=len(objects), contexts=len(contexts), width=width,
                          ranks=ranks(matrix), matrix_hex=[row.hex() for row in matrix],
                          distinct_behaviours=len(set(map(bytes, matrix))), maximum_SG=max(map(max, matrix)))
            verify_matrix(sample)
        # Failure/timeout must never return a partial matrix as a rank result.
        info, matrix = measure(Path("/bin/false"), [7], [(7, 1)], 1, 2, directory)
        assert info["status"] == "failed" and matrix is None
        spin = directory / "spin"
        spin.write_text("#!/usr/bin/env python3\nwhile True: pass\n"); spin.chmod(0o700)
        info, matrix = measure(spin, [7], [(7, 1)], 1, .02, directory)
        assert info["status"] == "timeout" and matrix is None
        for report in reports:
            for sample in report["samples"]:
                if sample["status"] != "complete":
                    assert not any(k in sample for k in ("ranks", "matrix_hex", "maximum_SG"))
                    continue
                matrix = verify_matrix(sample); matrices += 1
                if report["configuration"]["guarded"]:
                    verify_guarded(sample)
                # Select sparse representatives from EACH width, including the
                # widest saved matrices, without replaying an expensive full scan.
                candidates = []
                for i, left_hex in enumerate(sample["objects_hex"]):
                    for j, (right_hex, size) in enumerate(sample["contexts_hex"]):
                        left, right = int(left_hex, 16), int(right_hex, 16)
                        mask = left | (right << (3*sample["width"]))
                        candidates.append((mask.bit_count(), sample["width"]+size, mask, matrix[i][j]))
                # Keep the contexts of MAXIMUM width, so the check includes
                # long joins, not only the smallest right contexts.
                last = max(n for _, n, _, _ in candidates)
                for _, n, mask, value in sorted(c for c in candidates if c[1] == last)[:6]:
                    positions[n, mask] = value
        if positions:
            original_samples(directory, positions)
            # A second independent coordinate/subset recurrence for sparse cases.
            for (n, mask), value in positions.items():
                if mask.bit_count() <= 17:
                    assert sparse_value(mask) == value
    print(f"{small_checks} small entries agree with direct subset mex; {matrices} saved matrices and every rank minor checked")
    print(f"{len(positions)} saved entries agree with the independent row-major solver; guarded catalogues and every join are reachable")
    print("failed and timed-out jobs return no matrix; all conclusions remain finite")


if __name__ == "__main__":
    main()
