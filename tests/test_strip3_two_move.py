#!/usr/bin/env python3
"""Independent coordinate checks of the parameterised two-move identities.

The expected pieces below are written as coordinate sets, without importing
the builder. Small cases also use the independent recursive SG implementation.
Certificate mutation checks ensure that bad proofs cannot pass the verifier.
"""

import argparse
import copy
import json
from itertools import combinations
from pathlib import Path

from test_strip3_ports import rec_value
from verify_strip3_two_move import check
from verify_strip3_two_move import encode, family as premise_family

ROOT = Path(__file__).resolve().parents[1]


def full(n):
    return frozenset((r, c) for c in range(1, n + 1) for r in (1, 2, 3))


def model(name, n, j=None):
    cells = full(n)
    if name == "A":
        return cells - {(2, n)}
    if name == "D":
        return cells - {(2, 1), (2, n)}
    if name == "E":
        return cells - {(1, j - 1), (1, j), (1, j + 1), (2, j), (2, 1)}
    if name == "B":
        return cells - {(1, j - 1), (1, j), (1, j + 1), (2, j)}
    if name == "T":
        return cells - {(1, n)}
    if name == "C":
        return cells - {(1, n), (3, n)}
    raise ValueError(name)


def mapped(cells, n, offset=0, horizontal=False, vertical=False):
    return frozenset((4 - r if vertical else r,
                      offset + (n + 1 - c if horizontal else c))
                     for r, c in cells)


def distance(a, b):
    return abs(a[0] - b[0]) + abs(a[1] - b[1])


def tail(n, a, b):
    return frozenset(v for v in full(n)
                     if distance(v, a) > 1 and distance(v, b) > 1)


def assert_sum(n, a, b, pieces, sg):
    if distance(a, b) <= 1:
        raise AssertionError(("illegal move pair", n, a, b))
    successor = tail(n, a, b)
    used = frozenset()
    for piece in pieces:
        if used & piece or any(distance(v, w) == 1 for v in piece for w in used):
            raise AssertionError(("overlap or crossing edge", n, a, b))
        used |= piece
    if used != successor:
        raise AssertionError(("false set identity", n, a, b, used ^ successor))
    if sg:
        value = 0
        for piece in pieces:
            value ^= rec_value(piece)
        if rec_value(successor) != value:
            raise AssertionError(("false SG identity", n, a, b))


def identities(max_n, max_sg_n):
    counts = {"middle_middle": 0, "mixed": 0, "opposite": 0}
    for n in range(1, max_n + 1):
        sg = n <= max_sg_n
        for u in range(1, n + 1):
            for v in range(u + 2, n + 1):
                pieces = [model("A", u - 1),
                          mapped(model("D", v - u - 1), v - u - 1, u),
                          mapped(model("A", n - v), n - v, v, horizontal=True)]
                assert_sum(n, (2, u), (2, v), pieces, sg)
                counts["middle_middle"] += 1
        for j in range(1, n + 1):
            for i in range(1, n + 1):
                if i != j:
                    if i < j:
                        pieces = [model("A", i - 1),
                                  mapped(model("E", n - i, j - i), n - i, i)]
                    else:
                        pieces = [mapped(model("E", i - 1, i - j), i - 1,
                                         horizontal=True),
                                  mapped(model("A", n - i), n - i, i,
                                         horizontal=True)]
                    for row in (1, 3):
                        actual = [mapped(p, n, vertical=row == 3) for p in pieces]
                        assert_sum(n, (row, j), (2, i), actual, sg)
                        counts["mixed"] += 1
                delta = abs(i - j)
                if delta > 2:
                    continue
                lo, hi = min(i, j), max(i, j)
                if delta == 0:
                    pieces = [model("C", lo - 1),
                              mapped(model("C", n - hi), n - hi, hi,
                                     horizontal=True)]
                elif delta == 1:
                    pieces = [mapped(model("T", lo - 1), lo - 1, vertical=i < j),
                              mapped(model("T", n - hi), n - hi, hi,
                                     horizontal=True, vertical=j < i)]
                else:
                    pieces = [mapped(model("B", lo, 1), lo, horizontal=True,
                                     vertical=i < j),
                              frozenset({(2, lo + 1)}),
                              mapped(model("B", n - hi + 1, 1), n - hi + 1,
                                     hi - 1, vertical=j < i)]
                assert_sum(n, (1, j), (3, i), pieces, sg)
                counts["opposite"] += 1
    return counts


def connected(cells):
    if not cells:
        return True
    seen = {min(cells)}
    todo = list(seen)
    for r, c in todo:
        for w in ((r - 1, c), (r + 1, c), (r, c - 1), (r, c + 1)):
            if w in cells and w not in seen:
                seen.add(w)
                todo.append(w)
    return seen == cells


def nonseparating(max_n):
    """The missing pair types always leave one connected piece."""
    count = 0
    for n in range(1, max_n + 1):
        for a, b in combinations(sorted(full(n)), 2):
            same_outer = a[0] == b[0] and a[0] in (1, 3)
            distant_opposite = {a[0], b[0]} == {1, 3} and abs(a[1] - b[1]) >= 3
            if distance(a, b) <= 1 or not (same_outer or distant_opposite):
                continue
            successor = tail(n, a, b)
            if not successor or not connected(successor):
                raise AssertionError(("false connectedness claim", n, a, b))
            count += 1
    return count


def auxiliary_records():
    checked = 0
    for name in ("strip3_ports_extension.json", "strip3_two_move_auxiliary.json"):
        data = json.loads((ROOT / "data" / name).read_text())
        keys = set()
        for row in data["rows"]:
            key = row["family"], row["n"], row.get("j")
            cells = premise_family(*key)
            if key in keys or int(row["mask"], 16) != encode(cells):
                raise AssertionError("invalid auxiliary record geometry")
            keys.add(key)
            if row["value"] is not None and (type(row["value"]) is not int
                                              or row["value"] < 0):
                raise AssertionError("invalid auxiliary SG value")
            if "options" in row:
                moves = [tuple(v) for v in row["moves"]]
                values = row["options"]
                if len(moves) != len(values) or len(moves) != len(set(moves)) \
                        or set(moves) != cells:
                    raise AssertionError("auxiliary mex omitted a legal move")
                if all(type(v) is int and v >= 0 for v in values):
                    mex = 0
                    while mex in set(values):
                        mex += 1
                    if row["value"] != mex:
                        raise AssertionError("false auxiliary mex")
                elif row["value"] is not None:
                    raise AssertionError("partial options claimed an exact nimber")
            checked += 1
    values = {json.loads(line)["n"]: json.loads(line)["options"][0]
              for line in (ROOT / "data/strip3_survey.jsonl").read_text().splitlines()}
    for n in values:
        if n + 2 in values and values[n + 2] == values[n]:
            raise AssertionError("endpoint shift-two inequality violated")
    return checked


def width_bound(n, j, limit):
    """No zero separator response using a port wider than limit is considered.

    Within the limit every needed SG entry must exist; all geometric pieces
    are reconstructed independently before asserting the xor is nonzero.
    The width bound is conditional on these recorded exact SG values.
    """
    values = {("A", 0, None): 0, ("T", 0, None): 0, ("C", 0, None): 0,
              ("point", 1, None): 1}
    for line in (ROOT / "data/strip3_survey.jsonl").read_text().splitlines():
        row = json.loads(line)
        values["A", row["n"], None] = row["end_gap"]
        values["B", row["n"], 1] = row["options"][0]
    for line in (ROOT / "data/strip3_ports_n12.jsonl").read_text().splitlines():
        row = json.loads(line)
        if row.get("family") == "E":
            values["E", row["n"], row["j"]] = row["value"]
    for name in ("strip3_ports_extension.json", "strip3_two_move_auxiliary.json"):
        for row in json.loads((ROOT / "data" / name).read_text())["rows"]:
            if type(row["value"]) is int:
                values[row["family"], row["n"], row.get("j")] = row["value"]
    candidates = []
    for i in range(1, n + 1):
        if i == j:
            continue
        if i < j:
            keys = [("A", i - 1, None), ("E", n - i, j - i)]
            pieces = [model("A", i - 1), mapped(model("E", n - i, j - i), n - i, i)]
        else:
            keys = [("E", i - 1, i - j), ("A", n - i, None)]
            pieces = [mapped(model("E", i - 1, i - j), i - 1, horizontal=True),
                      mapped(model("A", n - i), n - i, i, horizontal=True)]
        candidates.append(((2, i), keys, pieces))
    for i in range(max(1, j - 2), min(n, j + 2) + 1):
        lo, hi = min(i, j), max(i, j)
        if i == j:
            keys = [("C", lo - 1, None), ("C", n - hi, None)]
            pieces = [model("C", lo - 1),
                      mapped(model("C", n - hi), n - hi, hi, horizontal=True)]
        elif abs(i - j) == 1:
            keys = [("T", lo - 1, None), ("T", n - hi, None)]
            pieces = [mapped(model("T", lo - 1), lo - 1, vertical=i < j),
                      mapped(model("T", n - hi), n - hi, hi,
                             horizontal=True, vertical=j < i)]
        else:
            keys = [("B", lo, 1), ("point", 1, None), ("B", n - hi + 1, 1)]
            pieces = [mapped(model("B", lo, 1), lo, horizontal=True, vertical=i < j),
                      frozenset({(2, lo + 1)}),
                      mapped(model("B", n - hi + 1, 1), n - hi + 1,
                             hi - 1, vertical=j < i)]
        candidates.append(((3, i), keys, pieces))
    checked = 0
    for reply, keys, pieces in candidates:
        if max(key[1] for key in keys) > limit:
            continue
        assert_sum(n, (1, j), reply, pieces, False)
        total = 0
        for key in keys:
            if key not in values:
                raise AssertionError(("missing width-bound premise", key))
            total ^= values[key]
        if total == 0:
            raise AssertionError(("false lower width bound", n, j, reply))
        checked += 1
    return checked


def rejects(report):
    try:
        check(report)
    except (ValueError, KeyError, TypeError):
        return
    raise AssertionError("invalid certificate accepted")


def mutations():
    original = json.loads((ROOT / "data/proofs/strip3_two_move_30.json").read_text())
    check(original)
    bad = copy.deepcopy(original)
    bad["cases"][0]["reply"] = bad["cases"][0]["first"]
    rejects(bad)
    bad = copy.deepcopy(original)
    bad["cases"][0]["tail_mask"] = "0"
    rejects(bad)
    bad = copy.deepcopy(original)
    bad["cases"][0]["pieces"][0]["isomorphism"]["horizontal"] = False
    rejects(bad)
    bad = copy.deepcopy(original)
    bad["cases"].pop()
    rejects(bad)
    bad = copy.deepcopy(original)
    bad["cases"][0]["pieces"].pop()
    rejects(bad)
    bad = copy.deepcopy(original)
    bad["cases"][0]["pieces"][0]["value"] ^= 1
    rejects(bad)
    bad = copy.deepcopy(original)
    bad["premises"][0]["value"] = True
    rejects(bad)
    bad = copy.deepcopy(original)
    bad["max_premise_width"] -= 1
    rejects(bad)
    return 8


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-n", type=int, default=30)
    parser.add_argument("--max-sg-n", type=int, default=7)
    args = parser.parse_args()
    print("parameterised geometry checks:", identities(args.max_n, args.max_sg_n))
    print(f"complete separator classification: {nonseparating(args.max_n)} "
          "other legal move pairs leave a connected successor")
    print(f"independent SG recursion checked all identities for n <= {args.max_sg_n}")
    print(f"auxiliary data: {auxiliary_records()} descriptors/recorded mexes checked; "
          "endpoint shift-two inequality holds on the saved prefix")
    print("zero separator width bounds (conditional on saved SG values): "
          f"3x30@(1,7) needs >=22 ({width_bound(30, 7, 21)} replies excluded); "
          f"3x36@(1,17) needs >=31 ({width_bound(36, 17, 30)} replies excluded)")
    print(f"certificate validation: {mutations()} false proofs rejected")


if __name__ == "__main__":
    main()
