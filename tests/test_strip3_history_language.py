#!/usr/bin/env python3
"""Independent graph/history checks of the complete connected-state language.

Histories, closed-neighbourhood deletion and connectivity below use coordinate
sets. They do not import the constructor's bit operations or gap recognizer.
"""

from functools import cache
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

from strip3_history_language import (LanguageIndex, LanguageSolver, component_count,
                                    decode, history_count, words)


def board(n):
    return {(r, c) for r in range(3) for c in range(n)}


def closed(vertex):
    r, c = vertex
    return {(r, c), (r - 1, c), (r + 1, c), (r, c - 1), (r, c + 1)}


def pieces(position):
    position = set(position)
    while position:
        seed = min(position)
        reached, stack = {seed}, [seed]
        position.remove(seed)
        while stack:
            vertex = stack.pop()
            new = (closed(vertex) - {vertex}) & position
            reached |= new
            position -= new
            stack.extend(new)
        yield reached


def mask(position, shift=True):
    offset = min(c for _, c in position) if shift and position else 0
    return sum(1 << (3 * (c - offset) + r) for r, c in position)


def vertices(bits):
    return {(i % 3, i // 3) for i in range(bits.bit_length()) if bits >> i & 1}


def independent_histories(n):
    # The possible independent subsets of one vertical P3, with horizontal
    # adjacency checked between columns. No connected-language gap rules.
    column_sets = ((), (0,), (1,), (2,), (0, 2))
    def extend(prefix, previous):
        if len(prefix) == n:
            yield tuple((r, c) for c, rows in enumerate(prefix) for r in rows)
            return
        for rows in column_sets:
            if not set(rows) & set(previous):
                yield from extend(prefix + (rows,), rows)
    yield from extend((), ())


def after(n, history):
    selected = set()
    removed = set()
    for vertex in history:
        assert vertex in board(n)
        assert not closed(vertex) & selected
        selected.add(vertex)
        removed |= closed(vertex)
    return board(n) - removed


def exhaustive_histories():
    seen = {w: set() for w in range(1, 10)}
    histories = checked = spanning = 0
    for n in range(1, 10):
        full_span = set()
        for history in independent_histories(n):
            histories += 1
            position = after(n, history)
            parts = tuple(pieces(position))
            if len(parts) == 1 and {c for _, c in position} == set(range(n)):
                full_span.add(mask(position))
            for part in parts:
                bits = mask(part)
                w = (bits.bit_length() + 2) // 3
                code = decode(bits)
                assert code.mask() == bits
                seen[w].add(bits)
                checked += 1
        assert len(full_span) == history_count(n), (n, len(full_span))
        assert full_span == {mask(after(n, tuple((b.bit_length() - 1, i)
                                               for i, b in enumerate(word) if b)))
                             for word in words(n)}
        spanning += len(full_span)
    # Width+4 suffices for all constructor witnesses, hence the exhaustive
    # union of actual histories through 9 is the entire family through 5.
    for w in range(1, 6):
        index = LanguageIndex(w)
        assert seen[w] == {index.unrank(i).mask() for i in range(index.size)}
    print(f"all legal histories through 3x9: {histories} histories, "
          f"{checked} actual components, {spanning} distinct spanning residuals")
    return seen, {"histories": histories, "components": checked, "spanning": spanning}


def bijection_and_witnesses():
    checked = moves = 0
    for n in range(1, 13):
        index = LanguageIndex(n)
        masks = set()
        for rank in range(index.size):
            code = index.unrank(rank)
            bits = code.mask()
            assert index.rank(bits) == rank
            masks.add(bits)
            assert len(tuple(pieces(vertices(bits)))) == 1
            width, history, target = code.witness()
            assert width <= n + 4
            assert target in {mask(part, shift=False) for part in pieces(after(width, history))}
            checked += 1
            if n <= 9:
                position = vertices(bits)
                for vertex in position:
                    for part in pieces(position - closed(vertex)):
                        child = mask(part)
                        assert decode(child).mask() == child
                    moves += 1
        assert len(masks) == component_count(n)
        if n >= 3:
            # A separate weighted sum over the unindexed outer language.
            def weight(word):
                return next((1 << i for i, bit in enumerate(word[:3]) if bit), 5)
            assert index.size == sum(weight(w) * weight(w[::-1]) for w in words(n))
    print(f"perfect index and explicit legal witnesses: {checked} masks; "
          f"all {moves} legal moves through width 9 remain in the family")
    return {"indexed_masks_and_witnesses": checked, "closed_moves": moves}


def rejection_checks(actual):
    checked = rejected = 0
    for bits in range(1, 1 << 15):
        if not bits & 7:
            continue
        if len(tuple(pieces(vertices(bits)))) != 1:
            continue
        w = (bits.bit_length() + 2) // 3
        try:
            code = decode(bits)
        except ValueError:
            assert bits not in actual[w]
            rejected += 1
        else:
            assert bits in actual[w] and code.mask() == bits
        checked += 1
    for w in (1, 2, 3, 12):
        index = LanguageIndex(w)
        for bad in (-1, index.size, True):
            try:
                index.unrank(bad)
            except ValueError:
                pass
            else:
                raise AssertionError("invalid index accepted")
    for bad in (0, -1, True, 8, 0o333):
        try:
            decode(bad)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid or unreachable mask accepted")
    print(f"classification of all {checked} connected column-trimmed 3x5 subsets: "
          f"{rejected} unreachable masks rejected")
    return {"connected_subsets_classified": checked, "unreachable_rejected": rejected}


def random_large():
    rng = random.Random(839286755)
    checked = 0
    for n in (10, 18, 33, 40, 80):
        index = LanguageIndex(n)
        for _ in range(120):
            rank = rng.randrange(index.size)
            code = index.unrank(rank)
            assert index.rank(code.mask()) == rank
            width, history, target = code.witness()
            assert target in {mask(p, shift=False) for p in pieces(after(width, history))}
            available = board(n)
            history = []
            # Also cover arbitrary mixed-row histories, including middle cuts.
            for _ in range(rng.randrange(3 * n + 1)):
                if not available:
                    break
                vertex = rng.choice(sorted(available))
                history.append(vertex)
                available -= closed(vertex)
            for part in pieces(available):
                bits = mask(part)
                assert decode(bits).mask() == bits
                checked += 1
    print(f"600 large index/witness checks and {checked} random actual components through width 80")
    return {"large_index_checks": 600, "random_large_components": checked}


def independent_sg():
    @cache
    def raw(position):
        if not position:
            return 0
        options = {raw(frozenset(set(position) - closed(vertex))) for vertex in position}
        g = 0
        while g in options:
            g += 1
        return g

    dense, sparse = LanguageSolver(9, dense=True), LanguageSolver(9)
    checked = option_checks = 0
    for n in range(1, 10):
        index = LanguageIndex(n)
        for rank in range(index.size):
            bits = index.unrank(rank).mask()
            position = frozenset(vertices(bits))
            expected = raw(position)
            assert dense.value(bits) == sparse.value(bits) == expected
            for vertex in position:
                remaining = set(position) - closed(vertex)
                assert dense.value(mask(remaining, shift=False)) == raw(frozenset(remaining))
                option_checks += 1
            checked += 1
    values = [2, 1, 1, 0, 3, 3, 2, 2, 2, 3, 3, 5]
    survey = LanguageSolver(12, dense=True)
    assert [survey.value((1 << (3 * n)) - 1) for n in range(1, 13)] == values
    print(f"independent coordinate-only subset SG: {checked} component values and "
          f"{option_checks} actual options; empty-board prefix through 12")
    return {"independent_sg_components": checked, "independent_sg_options": option_checks,
            "empty_board_verified_through": 12}


def c_implementation():
    rng = random.Random(1839286)
    positions, expected, ranks = [], [], []
    solver = LanguageSolver(10, dense=True)
    for n in range(1, 11):
        index = LanguageIndex(n)
        for rank in range(index.size):
            bits = index.unrank(rank).mask()
            positions.append((n, bits))
            expected.append(solver.value(bits))
            ranks.append(rank)
    large = []
    for n in (18, 29, 32, 33, 40, 42):
        index = LanguageIndex(n)
        for _ in range(100):
            rank = rng.randrange(index.size)
            large.append((n, index.unrank(rank).mask(), rank))
    with tempfile.TemporaryDirectory(prefix="strip3-language-check-") as folder:
        binary = Path(folder) / "language"
        subprocess.run(["gcc", "-O3", "-march=native", "-Wall", "-Wextra", "-pthread",
                        "-o", str(binary), str(ROOT / "src/cfast/strip3_language.c")], check=True)

        def run(command, input_text=None, **env):
            return subprocess.run([str(binary), *command], input=input_text,
                                  text=True, capture_output=True, timeout=90,
                                  env={**os.environ, **env})

        inputs = "".join(f"{n} {bits:x}\n" for n, bits in positions)
        for flags in range(8):
            result = run(["masks"], inputs, STRIP3_LANGUAGE_REDUCTIONS=str(flags))
            assert result.returncode == 0, result.stderr
            assert list(map(int, result.stdout.split())) == expected
        inputs += "".join(f"{n} {bits:x}\n" for n, bits, _ in large)
        result = run(["ranks"], inputs)
        assert result.returncode == 0, result.stderr
        assert list(map(int, result.stdout.split())) == ranks + [r for _, _, r in large]
        result = run(["strip", "22"])
        assert result.returncode == 0, result.stderr
        assert [list(map(int, line.split()))[1] for line in result.stdout.splitlines()] == [
            2, 1, 1, 0, 3, 3, 2, 2, 2, 3, 3, 5, 2, 4, 1, 3, 2, 2, 1, 1, 5, 0]
        negatives = [(["board", "1"], None, {"STRIP3_LANGUAGE_BYTES": "1"}, 3),
                     (["board", "33"], None, {"STRIP3_LANGUAGE_BYTES": "1024"}, 3),
                     (["board", "1"], None, {"STRIP3_LANGUAGE_BYTES": "-1"}, 2),
                     (["board", "1"], None, {"STRIP3_LANGUAGE_BYTES": "18446744073709551616"}, 2),
                     (["board", "1"], None, {"STRIP3_LANGUAGE_REDUCTIONS": "8"}, 2),
                     (["board", "43"], None, {}, 2),
                     (["masks"], "3 db\n", {}, 2),
                     (["ranks"], "2 2d\n", {}, 2),
                     (["masks"], "42 80000000000000000000000000000000\n", {}, 2)]
        for command, inputs, env, status in negatives:
            result = run(command, inputs, **env)
            assert result.returncode == status and not result.stdout, (command, result)
        sparse = run(["masks"], "42 20000000000000000000000000000000\n42 0\n")
        assert sparse.returncode == 0 and sparse.stdout.split() == ["1", "0"]
    print(f"C byte table: {len(positions)} SG values under all 8 reduction settings; "
          f"{len(ranks) + len(large)} perfect ranks; empty boards through 22; 9 rejection checks")
    return {"c_sg_positions_per_reduction_setting": len(positions), "c_reduction_settings": 8,
            "c_rank_checks": len(ranks) + len(large), "c_empty_board_through": 22,
            "c_negative_checks": len(negatives)}


def main():
    actual, report = exhaustive_histories()
    report.update(bijection_and_witnesses())
    report.update(rejection_checks(actual))
    report.update(random_large())
    report.update(independent_sg())
    report.update(c_implementation())
    report["counts_through_40"] = [component_count(n) for n in range(41)]
    report["histories_through_40"] = [history_count(n) for n in range(41)]
    report["status"] = "passed"
    output = ROOT / "data/proofs/strip3_history_language_audit.json"
    output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    print("all connected-language checks passed")


if __name__ == "__main__":
    main()
