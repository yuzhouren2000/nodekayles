"""Independent finite checks of cube rigidity, proof rules and strip parity DP."""
import itertools
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from history_cube_geometry import (cube_witness, independent_history_count,
                                   maximal_history_parities, neighbourhoods, parity_sg)
from test_strip3 import subset_table
from verify_strip3_cube_proof import check_text


def all_independent(adjacency):
    return [s for s in range(1 << len(adjacency))
            if all(not (s & adjacency[v]) for v in range(len(adjacency)) if s >> v & 1)]


def remainder(adjacency, history):
    gone = history
    for v, neighbours in enumerate(adjacency):
        if history >> v & 1:
            gone |= neighbours
    return ((1 << len(adjacency))-1) & ~gone


def subset_sg(adjacency):
    values = [0] * (1 << len(adjacency))
    for mask in range(1, len(values)):
        options = {values[mask & ~(adjacency[v] | (1 << v))]
                   for v in range(len(adjacency)) if mask >> v & 1}
        while values[mask] in options:
            values[mask] += 1
    return values


def sparse_sg(mask):
    vertices = [v for v in range(mask.bit_length()) if mask >> v & 1]
    adjacency = [sum(1 << j for j, u in enumerate(vertices) if
                     abs(v//3-u//3)+abs(v%3-u%3) == 1) for v in vertices]
    return subset_sg(adjacency)[-1]


def check_geometry():
    pairs = list(itertools.combinations(range(4), 2))
    checked = nontrivial = 0
    for edge_set in range(1 << len(pairs)):
        adjacency = [0]*4
        for i, (v, u) in enumerate(pairs):
            if edge_set >> i & 1:
                adjacency[v] |= 1 << u
                adjacency[u] |= 1 << v
        histories = all_independent(adjacency)
        sg = subset_sg(adjacency)
        for letters in itertools.product(range(3), repeat=4):
            yes = sum(1 << v for v, a in enumerate(letters) if a == 1)
            no = sum(1 << v for v, a in enumerate(letters) if a == 0)
            members = [s for s in histories if s & yes == yes and not s & no]
            if not members:
                try:
                    cube_witness(adjacency, yes, no)
                except ValueError:
                    continue
                raise AssertionError("accepted empty cube")
            checked += 1
            witness = cube_witness(adjacency, yes, no)
            if witness is None:
                assert members == [yes]
            else:
                nontrivial += 1
                a, b = witness
                assert a in members and b in members and (b ^ a).bit_count() == 1
                assert sg[remainder(adjacency, a)] != sg[remainder(adjacency, b)]
    adjacency = neighbourhoods(4)
    sg = subset_table(4)
    histories = all_independent(adjacency)
    assert independent_history_count(4) == len(histories)
    for mask in range(1 << 12):
        expected = 0
        for s in histories:
            if s & ~mask:
                continue
            blocked = s
            for v in range(12):
                if s >> v & 1:
                    blocked |= adjacency[v]
            if not mask & ~blocked:
                expected |= 1 << (s.bit_count() % 2)
        assert maximal_history_parities(mask) == expected
        if expected in (1, 2):
            assert parity_sg(mask) == sg[mask]
    # A single affine parity system cannot cover the entire SG=0 class either.
    assert [sg[remainder(adjacency, s)] for s in (0, 10, 130, 10 ^ 130)] == [0, 0, 0, 2]
    # There ARE succinct nonconvex exact classes: a disjoint union of edges.
    for m in range(1, 7):
        matching = [1 << (v ^ 1) for v in range(2*m)]
        game = subset_sg(matching)
        hs = all_independent(matching)
        assert len(hs) == 3**m
        for s in hs:
            assert game[remainder(matching, s)] == (m-s.bit_count()) % 2
        assert sum(game[remainder(matching, s)] == 0 for s in hs) == (3**m+1)//2
    print(f"{checked} nonempty cubes on all 64 four-vertex graphs; {nontrivial} explicit one-move obstructions")
    print("all 4096 induced 3x4 maximal-history parity transfers agree with independent enumeration/mex")


def main():
    check_geometry()
    rng = random.Random(37)
    adjacency = neighbourhoods(5)
    sg = subset_table(5)
    histories = all_independent(adjacency)
    rng.shuffle(histories)
    values = {h: sg[remainder(adjacency, h)] for h in histories}
    with tempfile.TemporaryDirectory(prefix="strip3-cube-test-") as folder:
        binary = Path(folder) / "cube"
        subprocess.run(["gcc", "-O3", "-Wall", "-Wextra", "-o", str(binary),
                        str(ROOT / "src/cfast/strip3_cube.c")], check=True, timeout=5)
        env = {k: v for k, v in os.environ.items() if not k.startswith("STRIP3_")}

        def run(mode, n, text=None, **extra):
            return subprocess.run([str(binary), mode, str(n)], input=text, text=True,
                                  capture_output=True, env={**env, **extra}, timeout=3)

        for mode in (0, 1, 2):
            queries = [(h, k) for h in histories for k in rng.sample(range(9), 9)]
            r = run("queries", 5, "".join(f"{h:x} {k}\n" for h, k in queries),
                    STRIP3_CUBE_LOOKUP=str(mode))
            assert r.returncode == 0, r.stderr
            assert list(map(int, r.stdout.split())) == [int(values[h] == k) for h, k in queries]
            assert "region_hits=0" in r.stderr, r.stderr
        r = run("histories", 5, "".join(f"{h:x}\n" for h in histories))
        assert r.returncode == 0 and list(map(int, r.stdout.split())) == [values[h] for h in histories]

        for n, expected in ((1, 2), (2, 1), (3, 1), (4, 0), (5, 3), (6, 3)):
            path = Path(folder) / f"cube-{n}.txt"
            r = subprocess.run([str(binary), "board", str(n), "--certificate", str(path)],
                               capture_output=True, text=True, env=env, timeout=3)
            assert r.returncode == 0, r.stderr
            result = check_text(path.read_text())
            assert result["n"] == n and result["SG"] == expected
            for line in path.read_text().splitlines()[1:]:
                row = line.split()
                yes, no = int(row[4], 16), int(row[5], 16)
                # In this inference system even the negative facts remain
                # singleton histories. The checker itself permits broad N.
                assert cube_witness(neighbourhoods(n), yes, no) is None
        text = path.read_text()

        def rejected(data):
            try:
                check_text(data)
            except (AssertionError, ValueError, IndexError):
                return
            raise AssertionError("accepted corrupted proof")

        rejected(text.rsplit("\n", 2)[0]+"\n")
        rejected(text+"extra\n")
        rows = text.splitlines()
        for operation in ("not_legal", "future", "missing", "wrong_target"):
            changed = list(rows)
            for i in range(1, len(rows)):
                row = rows[i].split()
                if row[3] != "M":
                    continue
                if operation == "not_legal":
                    row[5] = f"{int(row[5], 16) & ~(1 << int(row[6])):x}"
                elif operation == "future":
                    row[7] = row[0]
                elif operation == "missing":
                    row.pop()
                else:
                    row[1] = str(18 if int(row[1]) != 18 else 17)
                changed[i] = " ".join(row)
                break
            rejected("\n".join(changed)+"\n")

        # Near-terminal histories spanning 300 bits; these are positions, not
        # claims to have solved any complete 100-column board.
        for n in (22, 42, 100):
            adjacency = neighbourhoods(n)
            maximal = sum(1 << v for v in range(3*n) if (v//3+v%3) % 2 == 0)
            selected = [v for v in range(3*n) if maximal >> v & 1]
            samples = [maximal]
            for _ in range(12):
                h = maximal
                for v in rng.sample(selected, 5):
                    h &= ~(1 << v)
                samples.append(h)
            expected = [sparse_sg(remainder(adjacency, h)) for h in samples]
            r = run("histories", n, "".join(f"{h:x}\n" for h in samples))
            assert r.returncode == 0 and list(map(int, r.stdout.split())) == expected, r.stderr
        for extra in ({"STRIP3_CUBE_QUERY_LIMIT": "1"}, {"STRIP3_CUBE_BITS": "8"}):
            r = run("board", 8, **extra)
            assert r.returncode == 3 and not r.stdout and "incomplete" in r.stderr, r
        for mode, n, text, extra in (("board", 101, None, {}), ("histories", 3, "3\n", {}),
                                    ("histories", 3, "200\n", {}), ("queries", 3, "0 -1\n", {}),
                                    ("histories", 3, "0 extra\n", {}),
                                    ("board", 3, None, {"STRIP3_CUBE_LOOKUP": "3"})):
            r = run(mode, n, text, **extra)
            assert r.returncode == 2 and not r.stdout, r
    print(f"{3*len(histories)*9} shuffled target queries; six independent cube proofs and six corruptions checked")
    print("39 sparse histories across 64/126/300 bits agree; work/memory/input failures emit no false SG")


if __name__ == "__main__":
    main()
