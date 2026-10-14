#!/usr/bin/env python3
"""Independent geometric checks and actual C-cache membership for notes/27.

The coordinate history generator checks only independence, then tests full
spanning connectivity. It does not use the regular-language gap conditions.
C probes verify that the constructed states were ALREADY memoized when full
board SG search finished, with the optional comb rules both on and off.
"""

import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from protected_strip_search import family, protected_word, seed_count, seed_words


def closed(v):
    r, c = v
    return {(r, c), (r - 1, c), (r + 1, c), (r, c - 1), (r, c + 1)}


def independent_seeds(length):
    def visit(prefix):
        if len(prefix) == length:
            yield prefix
            return
        for bit in (0, 1, 4):
            if not bit & prefix[-1]:
                yield from visit(prefix + (bit,))
    yield from visit((1,))


def after(word):
    position = {(r, c) for c in range(len(word)) for r in range(3)}
    for c, bit in enumerate(word):
        if bit:
            v = (bit.bit_length() - 1, c)
            assert v in position
            position -= closed(v)
    return position


def spans_connected(position, n):
    if {c for r, c in position} != set(range(n)):
        return False
    seen, pending = {min(position)}, [min(position)]
    while pending:
        new = (closed(pending.pop()) - seen) & position
        seen |= new
        pending.extend(new)
    return seen == position


def geometry(position, n, cutoff):
    assert spans_connected(position, n)
    assert {r for r, c in position if c == 0} == {2}
    assert all((r, c) in position for c in range(n - cutoff - 1, n) for r in range(3))
    assert len((closed((1, n - 2)) - {(1, n - 2)}) & position) == 4
    transforms = [{(2-r if v else r, n-1-c if h else c) for r, c in position}
                  for v, h in ((False, False), (False, True), (True, False), (True, True))]
    assert all(shape != position for shape in transforms[1:])
    return min(sum(1 << (3*c+r) for r, c in shape) for shape in transforms)


def check_geometry():
    accepted = prefixes = rejected = 0
    for length in range(1, 11):
        expected_seeds = set(seed_words(length))
        actual_seeds = set()
        for seed in independent_seeds(length):
            word = seed + (0,) * 7
            if spans_connected(after(word), len(word)):
                actual_seeds.add(seed)
            else:
                rejected += 1
        assert actual_seeds == expected_seeds
        assert len(actual_seeds) == seed_count(length)
        for cutoff in range(2, 9):
            keys = set()
            n = length + cutoff + 2
            constructor = dict(family(n, cutoff))
            for seed in actual_seeds:
                word = protected_word(seed, cutoff)
                position = after(word)
                key = geometry(position, n, cutoff)
                assert key not in keys
                keys.add(key)
                bits = sum(1 << (3*c+r) for r,c in position)
                assert bits == constructor[word]
                accepted += 1
                # Every ascending prefix of this legal history keeps the
                # shield and the asymmetry, not just the final residual.
                partial = [0] * n
                for c, bit in enumerate(word):
                    if bit:
                        partial[c] = bit
                        geometry(after(partial), n, cutoff)
                        prefixes += 1
    assert [seed_count(n) for n in range(1, 7)] == [1, 1, 2, 4, 7, 13]
    assert all(seed_count(n) == sum(seed_count(n-i) for i in (1,2,3)) for n in range(4, 101))
    block_checks = 0
    for width in range(3, 10):
        block = {(r,c) for c in range(width) for r in range(3)}
        for v in block:
            assert spans_connected(block - {v}, width)
            block_checks += 1
    return {"coordinate_final_states": accepted, "coordinate_history_prefixes": prefixes,
            "independent_but_disconnected_seeds_rejected": rejected,
            "biconnected_shield_vertex_deletion_checks": block_checks,
            "cutoffs_checked": list(range(2,9)), "seed_lengths_checked": list(range(1,11))}


def check_cache(temp):
    source = temp / "cache_probe.c"
    source.write_text(f'''#define main strip3_cli_main
#include "{ROOT / 'src/cfast/strip3.c'}"
#undef main
int main(int argc, char **argv) {{
    if (argc != 2) return 2;
    int n = atoi(argv[1]); initialise();
    int g = value(full_board(n));
    size_t memo = entries;
    unsigned long long moves = tst.moves, calls = tst.calls;
    unsigned long long hi, lo; int checked = 0;
    while (scanf("%llx %llx", &hi, &lo) == 2) {{
        value(((Mask)hi << 64) | lo); checked++;
    }}
    printf("%d %zu %d %zu %llu %llu\\n", g, memo, checked,
           (size_t)entries-memo, (unsigned long long)tst.moves-moves,
           (unsigned long long)tst.calls-calls);
    return 0;
}}
''')
    binary = temp / "cache_probe"
    subprocess.run(["gcc", "-O3", "-march=native", "-Wall", "-Wextra", "-pthread",
                    "-o", str(binary), str(source)], check=True)
    rows = []
    for n in range(8, 21, 2):
        targets = list(family(n))
        text = "".join(f"{mask >> 64:x} {mask & ((1 << 64)-1):x}\n" for _, mask in targets)
        expected = seed_count(n-7)
        assert len(targets) == expected
        for flags in (0, 7, 15, 23, 31):
            result = subprocess.run([str(binary), str(n)], input=text, text=True, capture_output=True, timeout=60,
                                    env={**os.environ, "STRIP3_REDUCTIONS": str(flags), "STRIP3_THREADS": "1",
                                         "STRIP3_HT_BITS": "24", "STRIP3_HT_PREALLOC": "1"})
            assert result.returncode == 0, result.stderr
            g, memo, checked, added, moves, calls = map(int, result.stdout.split())
            assert checked == expected and added == moves == 0 and calls == checked, (n, flags, result.stdout)
            rows.append({"n": n, "reductions": flags, "SG": g, "board_memo": memo,
                         "protected_states_already_cached": checked, "additional_memo": added,
                         "additional_moves": moves})
        assert len({row['SG'] for row in rows if row['n']==n}) == 1
    return rows


def main():
    audit = {"format": "node-kayles-bounded-thickness-search-lower-bound-audit-v1",
             **check_geometry()}
    with tempfile.TemporaryDirectory(prefix="protected-strip-check-") as folder:
        audit["C_cache_probes"] = check_cache(Path(folder))
    audit["C_source_sha256"] = hashlib.sha256((ROOT / "src/cfast/strip3.c").read_bytes()).hexdigest()
    audit["total_C_cache_membership_checks"] = sum(row["protected_states_already_cached"] for row in audit["C_cache_probes"])
    output = ROOT / "data/proofs/protected_strip_search_audit.json"
    output.write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps({key:value for key,value in audit.items() if key!='C_cache_probes'}))


if __name__ == "__main__":
    main()
