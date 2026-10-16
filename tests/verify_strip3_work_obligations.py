"""Replay the unfinished query chain using independent coordinate geometry.

This verifies a diagnostic, not the SG value of the unfinished root query.
Only the small detached pieces are evaluated by independent subset mex.
"""
from pathlib import Path
import json

from test_strip3 import sparse_value

ROOT = Path(__file__).resolve().parents[1]


def cells(mask):
    return {(v//3, v%3) for v in range(mask.bit_length()) if mask >> v & 1}


def canonical(points):
    c0 = min(c for c, _ in points); r0 = min(r for _, r in points)
    pts = {(c-c0, r-r0) for c, r in points}
    n = 1+max(c for c, _ in pts); h = 1+max(r for _, r in pts)
    bases = [pts]
    if n == h == 3:
        bases.append({(r, c) for c, r in pts})
    return min(sum(1 << (3*(n-1-c if fc else c)+(h-1-r if fr else r)) for c, r in base)
               for base in bases for fc in (False, True) for fr in (False, True))


def components(points):
    unseen = set(points)
    while unseen:
        start = unseen.pop(); comp = {start}; todo = [start]
        while todo:
            c, r = todo.pop()
            for v in ((c-1,r), (c+1,r), (c,r-1), (c,r+1)):
                if v in unseen:
                    unseen.remove(v); comp.add(v); todo.append(v)
        yield comp


def main():
    report = json.loads((ROOT / "data/strip3_work_obligations.json").read_text())
    row, = report["samples"]
    assert row["status"] == "failed" and row["returncode"] == 3 and not row["stdout"]
    assert row["counts"]["moves"] == 750000
    chain = row["obligations"]
    root = chain[0]
    assert (root["n"], root["nodes"], root["k"]) == (28, 84, 4)
    assert int(root["mask"], 16) == (1 << 84)-1
    detached = 0
    for parent, child in zip(chain, chain[1:]):
        pm, cm = int(parent["mask"], 16), int(child["mask"], 16)
        assert pm.bit_count() == parent["nodes"] and cm.bit_count() == child["nodes"]
        if pm == cm:
            assert child["k"] < parent["k"]
            continue
        v = parent["move"]; assert v >= 0 and pm >> v & 1
        pc, pr = divmod(v, 3)
        remaining = {(c,r) for c,r in cells(pm) if abs(c-pc)+abs(r-pr) > 1}
        parts = [canonical(comp) for comp in components(remaining)]
        # Equal components cancel, independent of their placement.
        parts = [part for part in set(parts) if parts.count(part)%2]
        assert cm in parts, (parent, child)
        parts.remove(cm)
        offset = 0
        for part in parts:
            assert part.bit_count() <= 12  # keep the independent replay tiny
            offset ^= sparse_value(part); detached += 1
        assert child["k"] == parent["k"] ^ offset
    a24 = chain[1]
    assert int(a24["mask"], 16) == ((1 << 72)-1) ^ (1 << 70)
    old = {d["n"]: d for d in map(json.loads, (ROOT / "data/strip3_survey.jsonl").read_text().splitlines())}
    assert old[3]["end_gap"] == 1 and old[24]["end_gap"] == 9
    print(f"{len(chain)-1} query-dependency steps and {detached} small detached components replayed independently")
    print("the unfinished root target 4 reaches E(A(24),5); stored SG(A(24))=9, so this obligation is false")


if __name__ == "__main__":
    main()
