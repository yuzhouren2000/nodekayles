"""Independent subset-mex, raw recognition and integrated equality checks."""
import os
from pathlib import Path
import random
import re
import resource
import subprocess
import tempfile
import time

from test_strip3 import subset_table, sparse_value
from verify_strip3_caterpillar import adjacency, caterpillar_spine, interval_value, local_check, audit_corpus

ROOT = Path(__file__).resolve().parents[1]


def main():
    os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
    started = time.monotonic()

    def run(command, text=None, **extra):
        assert time.monotonic()-started < 20, "test batch budget"
        def limits():
            resource.setrlimit(resource.RLIMIT_AS, (512*1024**2,)*2)
            resource.setrlimit(resource.RLIMIT_CPU, (3,)*2)
        env = {k: v for k, v in os.environ.items() if not k.startswith("STRIP3_")}
        env.update(STRIP3_TARGET_BITS="18", STRIP3_TARGET_FRONTIER="4",
                   STRIP3_TARGET_SHARE_OPTIONS="1", STRIP3_TARGET_CACHE_FIRST="1")
        return subprocess.run(list(map(str, command)), input=text, text=True, capture_output=True,
                              timeout=2, preexec_fn=limits, env={**env, **extra})

    expected = subset_table(5)
    rng = random.Random(3800)
    masks = list(range(len(expected))); rng.shuffle(masks)
    direct = []
    for mask in masks:
        spine = caterpillar_spine(adjacency(mask))
        direct.append(-1 if spine is None else interval_value(spine))
        if spine is not None:
            assert direct[-1] == expected[mask]
    with tempfile.TemporaryDirectory(prefix="strip3-caterpillar-") as folder:
        folder = Path(folder)
        binary, raw_binary = folder / "candidate", folder / "recognition"
        result = run(["gcc", "-O3", "-Wall", "-Wextra", "-pthread", "-o", binary,
                      ROOT / "src/cfast/strip3_caterpillar.c"])
        assert result.returncode == 0, result.stderr
        raw_source = folder / "raw.c"
        raw_source.write_text('''#define STRIP3_CATERPILLAR_EXTERNAL_MAIN
#include "%s"
int main(void) {
    for (int c=0; c<MAXCOLS; c++) for (int r=0; r<3; r++) rowmask[r] |= (Mask)1 << (3*c+r);
    unsigned long long hi, lo;
    while (scanf("%%llx %%llx", &hi, &lo) == 2) printf("%%d\\n", caterpillar_value(((Mask)hi << 64)|lo));
    return 0;
}
''' % (ROOT / "src/cfast/strip3_caterpillar.c"))
        result = run(["gcc", "-O3", "-Wall", "-Wextra", "-pthread", "-o", raw_binary, raw_source])
        assert result.returncode == 0, result.stderr
        text = "".join(f"{mask >> 64:x} {mask & ((1 << 64)-1):x}\n" for mask in masks)
        result = run([raw_binary], text)
        assert result.returncode == 0 and list(map(int, result.stdout.split())) == direct, result.stderr
        text = "".join(f"5 {mask:x}\n" for mask in masks)
        result = run([binary, "masks"], text)
        assert result.returncode == 0 and list(map(int, result.stdout.split())) == [expected[m] for m in masks], result.stderr
        queries = [(mask, k) for mask in masks for k in (8, expected[mask], 0, 127)]
        rng.shuffle(queries)
        result = run([binary, "queries"], "".join(f"5 {m:x} {k}\n" for m, k in queries))
        assert result.returncode == 0 and list(map(int, result.stdout.split())) == [int(expected[m] == k) for m, k in queries], result.stderr

        sparse = [rng.randrange(1 << 15) << offset for offset in (54, 60, 63, 102, 111) for _ in range(5)]
        # A C8 and an unreachable subdivided claw must be rejected by the raw
        # tree evaluator. A disconnected cycle+isolated vertex has E=V-1.
        c8 = ((1 << 9)-1) & ~(1 << 4)
        claw = sum(1 << (3*c) for c in range(5)) | (6 << 6)
        negatives = [c8, claw, c8 | (1 << 15), 0]
        result = run([raw_binary], "".join(f"0 {m:x}\n" for m in negatives))
        assert result.returncode == 0 and list(map(int, result.stdout.split())) == [-1]*len(negatives)
        result = run([binary, "masks"], "".join(f"42 {m:x}\n" for m in sparse))
        assert result.returncode == 0 and list(map(int, result.stdout.split())) == [sparse_value(m) for m in sparse]
        # A long tree exercises the maximum spine/storage, independently of
        # full-board search. Both endpoints have one attached middle-row leaf.
        long_tree = sum(1 << (3*c) for c in range(42)) | 2 | (2 << 123)
        long_spine = caterpillar_spine(adjacency(long_tree))
        result = run([raw_binary], f"{long_tree >> 64:x} {long_tree & ((1 << 64)-1):x}\n")
        assert result.returncode == 0 and int(result.stdout) == interval_value(long_spine)
        # Turning the experiment off preserves the established engine exactly.
        result = run([binary, "board", 24], STRIP3_CATERPILLAR="0")
        assert result.returncode == 0 and result.stdout.split() == ["3", "24", "1"]
        stats = dict(re.findall(r"(\w+)=(\d+)", result.stderr))
        assert int(stats["memo"]) == 5544 and int(stats["moves"]) == 50220
        assert run([binary, "board", 3], STRIP3_CATERPILLAR="bad").returncode == 2
        result = run([binary, "board", 16], STRIP3_TARGET_MOVE_LIMIT="10")
        assert result.returncode == 3 and not result.stdout and "incomplete" in result.stderr
    local = local_check()
    corpus = audit_corpus(ROOT / "data/proofs/strip3_equality/strip3_equality_target_tree_min_26.json.gz")
    print(f"32768 arbitrary induced masks: raw recognition, interval values and integrated SG agree; {len(queries)} shuffled equality queries agree")
    print(f"{len(sparse)} cross-bit masks, long spine, C8, subdivided claw, disconnected false positive, disabled counts and budget errors checked")
    print(f"local proof coverage: {local}; saved-corpus audit: {corpus}")
    print(f"single-core test batch wall_s={time.monotonic()-started:.3f}")


if __name__ == "__main__":
    main()
