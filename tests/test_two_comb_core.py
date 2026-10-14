#!/usr/bin/env python3
"""Independent original-game and finite-spectrum tests for two-arm cores."""

from copy import deepcopy
from functools import cache
import gzip
from itertools import combinations
import json
from pathlib import Path
import random
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from two_comb_core import CertifiedTwoComb, TwoCombCore, tail_graph
from kayles_boundary_transfer import column_graph
from strip3_history_language import decode
from verify_two_comb_core import Checker


def original(core, roots, phases, lengths):
    """Named-vertex graph and raw legal moves, with no comb formula."""
    edges = {("core", v): set() for v in range(len(core))}
    def edge(v, u):
        edges.setdefault(v, set()).add(u)
        edges.setdefault(u, set()).add(v)
    for v, ns in enumerate(core):
        for u in range(v):
            if ns >> u & 1:
                edge(("core", v), ("core", u))
    for side, (root, phase, length) in enumerate(zip(roots, phases, lengths)):
        if length:
            edge(("core", root), ("spine", side, 0))
        for i in range(length):
            if i:
                edge(("spine", side, i - 1), ("spine", side, i))
            if phase ^ (i % 2):
                edge(("spine", side, i), ("leaf", side, i))
    @cache
    def sg(position):
        if not position:
            return 0
        reached = {next(iter(position))}
        pending = list(reached)
        while pending:
            v = pending.pop()
            new = (edges[v] & position) - reached
            reached |= new
            pending.extend(new)
        if len(reached) != len(position):
            return sg(frozenset(reached)) ^ sg(position - reached)
        values = {sg(position - edges[v] - {v}) for v in position}
        result = 0
        while result in values:
            result += 1
        return result
    initial = frozenset(edges)
    return sg(initial), {v: sg(initial - edges[v] - {v}) for v in initial}


def load(name):
    return json.loads(gzip.decompress((ROOT / f"data/proofs/two_comb_{name}.json.gz").read_bytes()))


def original_checks():
    count = 0
    for k in (2, 3):
        possible = list(combinations(range(k), 2))
        for code in range(1 << len(possible)):
            graph = [0] * k
            for i, (u, v) in enumerate(possible):
                if code >> i & 1:
                    graph[u] |= 1 << v
                    graph[v] |= 1 << u
            roots = (0, k - 1)
            for phases in ((0, 0), (0, 1), (1, 0), (1, 1)):
                model = TwoCombCore(graph, roots, phases).build(6)
                for n in range(7):
                    for m in range(7):
                        expected, options = original(graph, roots, phases, (n, m))
                        assert model.value(model.full, (n, m)) == expected
                        actual = model.literal_options(model.full, (n, m))
                        assert {value for _, value in actual} == set(options.values())
                        for move, value in actual:
                            assert options[move] == value
                        count += 1
    rng = random.Random(231026)
    for _ in range(40):
        k = rng.randrange(2, 7)
        graph = [0] * k
        for u, v in combinations(range(k), 2):
            if rng.randrange(3) == 0:
                graph[u] |= 1 << v
                graph[v] |= 1 << u
        roots = tuple(rng.sample(range(k), 2))
        phases = tuple(rng.randrange(2) for _ in range(2))
        model = TwoCombCore(graph, roots, phases).build(8)
        for _ in range(10):
            lengths = tuple(rng.randrange(9) for _ in range(2))
            expected, options = original(graph, roots, phases, lengths)
            assert model.value(model.full, lengths) == expected
            actual = model.literal_options(model.full, lengths)
            assert {v for _, v in actual} == set(options.values())
            count += 1
    return count


def literal_candidate_options(checker, mask, n, m):
    # Expand the actual tail positions; compare to the verifier's residue
    # interval compression even on tables whose values are not valid SGs.
    options = set()
    for v, ns in enumerate(checker.graph):
        if not mask >> v & 1:
            continue
        rest, lengths, extra = mask & ~(ns | 1 << v), [n, m], 0
        for side, (root, p, length) in enumerate(zip(checker.roots, checker.phases, (n, m))):
            if length and v == root:
                lengths[side] = 0
                extra ^= ((0, 1, 2, 1), (0, 1, 2, 3))[p ^ 1][(length - 1) % 4] ^ p
        options.add(checker.at(rest, *lengths) ^ extra)
    for side, (root, p, length, other) in enumerate(zip(checker.roots, checker.phases, (n, m), (m, n))):
        if not length:
            continue
        new = [n, m]
        new[side] = 0
        c = lambda phase, size: ((0, 1, 2, 1), (0, 1, 2, 3))[phase][size % 4]
        options.add(checker.terminal(mask & ~(1 << root), tuple(new))
                    ^ c(p, max(length - 2, 0)) ^ ((1 - p) if length >= 2 else 0))
        for i in range(length):
            tooth = p ^ (i % 2)
            if tooth:
                a, b = (i, other) if side == 0 else (other, i)
                options.add(checker.at(mask, a, b) ^ c(0, length - i - 1))
            if i:
                a, b = (i - 1, other) if side == 0 else (other, i - 1)
                options.add(checker.at(mask, a, b) ^ c(p ^ ((i + 2) % 2), max(length - i - 2, 0))
                            ^ ((1 + int(i + 1 < length)) * (1 - tooth) % 2))
    return options


def main():
    original_count = original_checks()
    rng = random.Random(2342)
    coefficient_checks = refined_checks = reachable = 0
    query_original_checks = query_spectrum_checks = 0
    for name in ("square_opposite", "strip_middle2", "strip_middle3"):
        proof = load(name)
        query = CertifiedTwoComb(proof)
        checker = Checker(proof, check_unary=False)
        for n in range(7):
            for m in range(7):
                lengths = n, m
                expected, raw_options = original(checker.graph, checker.roots, checker.phases, lengths)
                options = query.representative_options(lengths)
                assert query.value(lengths) == expected
                assert {value for _, value in options} == set(raw_options.values())
                assert all(raw_options[move] == value for move, value in options)
                winning = query.winning_move(lengths)
                assert (winning is None) == (expected == 0)
                if winning is not None:
                    assert raw_options[winning] == 0
                query_original_checks += 1
        for _ in range(80):
            lengths = rng.randrange(10 ** 18), rng.randrange(10 ** 18)
            # Include the precise diagonal and either short axis, not just far cells.
            if query_spectrum_checks % 4 == 0:
                lengths = lengths[0], lengths[0] + rng.randrange(-3, 4)
            elif query_spectrum_checks % 4 == 1:
                lengths = rng.randrange(8), lengths[1]
            elif query_spectrum_checks % 4 == 2:
                lengths = lengths[0], rng.randrange(8)
            options = query.representative_options(lengths)
            assert query.value(lengths) == checker.at(checker.full, *lengths)
            assert {value for _, value in options} == checker.options(checker.full, *lengths)
            for move, _ in options:
                if move[0] != "core":
                    kind, side, index = move
                    assert 0 <= index < lengths[side]
                    assert kind == "spine" or checker.phases[side] ^ (index % 2)
            winning = query.winning_move(lengths)
            assert (winning is None) == (query.value(lengths) == 0)
            if winning is not None:
                assert (winning, 0) in options
            query_spectrum_checks += 1
    for name in ("square_opposite", "strip_middle2", "strip_middle3", "strip_middle4"):
        checker = Checker(load(name), check_unary=False)
        for mask in checker.rows:
            for _ in range(12):
                n, m = rng.randrange(300), rng.randrange(300)
                assert literal_candidate_options(checker, mask, n, m) == checker.options(mask, n, m)
                coefficient_checks += 1
                n, m = rng.randrange(10 ** 12), rng.randrange(10 ** 12)
                a, b = checker.refined(n, m)
                assert checker.candidate(mask, n, m) == checker.candidate(mask, a, b)
                assert checker.options(mask, n, m) == checker.options(mask, a, b)
                refined_checks += 1
    # Random semilinear tables verify interval compression independently of
    # the actual SG values; it must not rely on accidental spectrum overlap.
    fake = Checker(load("square_opposite"), check_unary=False)
    for mask, rows in fake.rows.items():
        coordinates = {fake.canonical(n, m) for n in range(len(rows)) for m in range(len(rows))}
        values = {position: rng.randrange(64) for position in coordinates}
        fake.rows[mask] = [[values[fake.canonical(n, m)] for m in range(len(rows))]
                          for n in range(len(rows))]
    for _ in range(1000):
        n, m = rng.randrange(500), rng.randrange(500)
        assert literal_candidate_options(fake, fake.full, n, m) == fake.options(fake.full, n, m)
        a, b = fake.refined(n, m)
        assert fake.options(fake.full, n, m) == fake.options(fake.full, a, b)
        coefficient_checks += 1
        refined_checks += 1

    # Explicit real histories on EMPTY three-row boards. Arbitrary positive
    # arm lengths retain all three rows in the complete central thick block.
    families = []
    for width in (2, 3, 4, 5):
        checker = Checker(load(f"strip_middle{width}"), check_unary=False)
        for n in range(1, 16):
            for m in range(1, 16):
                parent_width = n + width + m + 2
                available = {(r, c) for c in range(parent_width) for r in range(3)}
                history = [(0, c) for c in range(n - 1, -1, -2)]
                history += [(0, c) for c in range(n + width + 2, parent_width, 2)]
                for r, c in history:
                    assert (r, c) in available
                    available -= {(r, c), (r - 1, c), (r + 1, c), (r, c - 1), (r, c + 1)}
                arm_left = tuple(4 if i % 2 == 0 else 6 for i in range(n))[::-1]
                arm_right = tuple(4 if i % 2 == 0 else 6 for i in range(m))
                cols = arm_left + (6,) + (7,) * width + (6,) + arm_right
                actual = tuple(sum(1 << r for r in range(3) if (r, c) in available)
                               for c in range(parent_width))
                assert actual == cols
                core_columns = (6,) + (7,) * width + (6,)
                graph, cells = column_graph(core_columns)
                roots = (cells.index((2, 0)), cells.index((2, width + 1)))
                phases = checker.phases
                assert phases == (0, 0)
                # Certificates must use phase zero for this geometric history.
                assert graph == checker.graph and roots == checker.roots
                if n + m <= 10:
                    expected, _ = original(graph, roots, phases, (n, m))
                    assert checker.at(checker.full, n, m) == expected
                reachable += 1
        matrix = [[checker.at(checker.full, 8 + i, 8 + j) for j in range(4)] for i in range(4)]
        if width == 3:
            residue_value = (0, 3, 2, 1)
            assert all(matrix[i][j] == (residue_value[i] ^ residue_value[j] ^ 1)
                       for i in range(4) for j in range(4))
            # Closed-form winning play uses only the last two spine vertices
            # (or last leaf), so no search of a long arm or core is required.
            query = CertifiedTwoComb(load("strip_middle3"))
            for n in range(12, 28):
                for m in range(12, 28):
                    value = query.value((n, m))
                    if value == 0:
                        continue
                    if value == 2:
                        side = 0 if n % 2 == 0 else 1
                        length = (n, m)[side]
                        move = ("spine", side, length - 1)
                    elif value == 1:
                        move = ("leaf", 0, n - 1) if n % 2 == 0 else ("spine", 0, n - 2)
                    else:
                        assert value == 3
                        move = ("spine", 0, n - 2) if n % 2 == 0 else ("spine", 0, n - 1)
                    kind, side, index = move
                    length, other = (n, m)[side], (m, n)[side]
                    if kind == "leaf":
                        prefix, extra = index, ((0, 1, 2, 1)[(length - index - 1) % 4])
                    else:
                        prefix = index - 1
                        pure_tail = ((0, 1, 2, 1), (0, 1, 2, 3))[(index + 2) % 2]
                        extra = pure_tail[max(length - index - 2, 0) % 4]
                        extra ^= (1 + int(index + 1 < length)) * (1 - index % 2) % 2
                    pair = (prefix, other) if side == 0 else (other, prefix)
                    assert checker.at(checker.full, *pair) ^ extra == 0
        assert all(checker.at(checker.full, n, m) == matrix[n % 4][m % 4]
                   for n in range(checker.b, checker.R + 1)
                   for m in range(checker.b, checker.R + 1))
        # Phase-zero reachable examples have independent periods even in every
        # retained descendant. The generic certificate also permits a diagonal band.
        assert all(checker.at(mask, n, m) == checker.at(mask, 8 + n % 4, 8 + m % 4)
                   for mask in checker.rows for n in range(8, checker.R + 1)
                   for m in range(8, checker.R + 1))
        families.append({"complete_middle_columns": width, "minimum": checker.b,
                         "matrix_mod4": matrix,
                         "zero_pairs_mod4": [[i, j] for i in range(4) for j in range(4)
                                             if matrix[i][j] == 0]})

    # Phase-one C4 with opposite roots has exact length discrimination forever.
    square = Checker(load("square_opposite"), check_unary=False)
    n = 10 ** 18 + 1
    assert square.at(square.full, n, n) == 0
    assert square.at(square.full, n + 4, n) == 2
    assert square.at(square.full, n, n + 2) == 5
    # Absorb the first phase-one vertex and its leaf into a bipartite core.
    # This gives the same persistent obstruction with two phase-zero arms.
    absorbed_graph = (26, 5, 74, 5, 33, 16, 132, 64)
    absorbed = TwoCombCore(absorbed_graph, (4, 6), (0, 0)).build(20)
    assert all(absorbed.value(absorbed.full, (i, j))
               == square.at(square.full, i + 1, j + 1)
               for i in range(21) for j in range(21))
    tampered = []
    proof = load("square_opposite")
    bad = deepcopy(proof); bad["states"][0]["rows"][8][8] ^= 1; tampered.append(bad)
    bad = deepcopy(proof); bad["states"].pop(); tampered.append(bad)
    bad = deepcopy(proof); bad["unary"].pop(); tampered.append(bad)
    bad = deepcopy(proof); bad["unary"][0]["states"][0]["values"][0] ^= 1; tampered.append(bad)
    for bad in tampered:
        try:
            Checker(bad).verify()
        except (ValueError, StopIteration):
            continue
        raise AssertionError("altered certificate accepted")
    audit = {"format": "node-kayles-two-comb-audit-v1", "original_graph_checks": original_count,
             "literal_coefficient_checks": coefficient_checks, "refined_spectrum_checks": refined_checks,
             "reachable_histories": reachable, "families": families,
             "query_original_checks": query_original_checks,
             "query_spectrum_checks": query_spectrum_checks,
             "tamper_rejections": len(tampered),
             "persistent_difference_obstruction": {"n": n, "equal": 0,
                                                   "difference4": 2, "difference2": 5}}
    (ROOT / "data/proofs/two_comb_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps(audit))


if __name__ == "__main__":
    main()
