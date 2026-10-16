#!/usr/bin/env python3
"""Check chronological history-cube inferences without SG or BDD evaluation.

Each cube is interpreted only on independent histories of the coordinate grid.
The checker verifies set implications under independence and legal substitution,
then applies mex to already checked premises. No source from the solver is used.
"""
import argparse
import gzip
from pathlib import Path


def check_text(text):
    lines = text.splitlines()
    assert lines, "empty proof"
    header = lines[0].split()
    assert len(header) == 4 and header[0] == "node-kayles-history-cubes-v1", "format"
    n, count, root = map(int, header[1:])
    assert 1 <= n <= 100 and 1 <= count <= 1_048_575 and 1 <= root <= count, "header bounds"
    assert len(lines) == count + 1, "record count"
    vertices = [(r, c) for c in range(n) for r in range(3)]
    size, full = len(vertices), (1 << len(vertices)) - 1
    closed = [sum(1 << u for u, (a, b) in enumerate(vertices)
                  if abs(a-r)+abs(b-c) <= 1) for r, c in vertices]
    open_neighbours = [mask ^ (1 << v) for v, mask in enumerate(closed)]
    facts = [None]

    def implication(parent, child, move=None):
        yes, _, zeros, _, _ = parent
        cyes, cno, _, _, _ = child
        if move is not None:
            assert not cno >> move & 1, "successor forbids chosen move"
            cyes &= ~(1 << move)
            # Legal moves force their open neighbourhood to zero.
            zeros |= open_neighbours[move]
        assert not cyes & ~yes and not cno & ~zeros, "condition does not imply premise"

    def premise(token, target, truth):
        index = int(token)
        assert 1 <= index < len(facts), "nonchronological reference"
        fact = facts[index]
        assert fact[3:] == (target, truth), "premise target/truth"
        return fact

    p_count = edges = literals = 0
    for expected, line in enumerate(lines[1:], 1):
        row = line.split()
        assert len(row) >= 6, "short record"
        index, k, truth = map(int, row[:3])
        kind = row[3]
        yes, no = map(lambda value: int(value, 16), row[4:6])
        assert index == expected and 0 <= k <= size and truth in (0, 1), "record bounds"
        assert yes >= 0 and no >= 0 and not (yes | no) & ~full and not yes & no, "condition bounds"
        zeros = no
        for v in range(size):
            if yes >> v & 1:
                assert not yes & open_neighbours[v], "nonindependent positive literals"
                zeros |= open_neighbours[v]
        parent = (yes, no, zeros, k, truth)
        refs = row[6:]
        literals += yes.bit_count() + no.bit_count()
        if kind == "P":
            assert truth == 1 and len(refs) == size + k, "universal premise count"
            p_count += 1
            for v in range(size):
                ref = int(refs[v])
                if ref < 0:
                    blocker = -ref - 1
                    assert 0 <= blocker < size and yes >> blocker & 1, "missing blocker"
                    assert closed[v] >> blocker & 1, "blocker not adjacent"
                else:
                    child = premise(refs[v], k, 0)
                    assert not yes & closed[v], "unnecessary move premise"
                    implication(parent, child, v)
                    edges += 1
            for j in range(k):
                implication(parent, premise(refs[size+j], j, 0))
                edges += 1
        elif kind == "M":
            assert truth == 0 and len(refs) == 2, "move premise count"
            v = int(refs[0])
            assert 0 <= v < size and not yes & closed[v], "illegal witness"
            assert not closed[v] & ~zeros, "move not forced legal"
            implication(parent, premise(refs[1], k, 1), v)
            edges += 1
        elif kind == "L":
            assert truth == 0 and len(refs) == 1, "lower-target premise count"
            index = int(refs[0])
            assert 1 <= index < len(facts), "nonchronological reference"
            child = facts[index]
            assert child[4] == 1 and child[3] < k, "not a smaller exact value"
            implication(parent, child)
            edges += 1
        else:
            raise AssertionError("unknown inference")
        facts.append(parent)
    result = facts[root]
    assert result[4] == 1 and result[0] == 0, "root does not prove an exact full-board value"
    return dict(n=n, SG=result[3], facts=count, p_facts=p_count,
                premise_edges=edges, literals=literals)


def verify(path):
    path = Path(path)
    text = gzip.decompress(path.read_bytes()).decode() if path.suffix == ".gz" else path.read_text()
    return check_text(text)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("proof", nargs="+", type=Path)
    args = parser.parse_args()
    for path in args.proof:
        print(path, verify(path))


if __name__ == "__main__":
    main()
