"""Check universal response certificates without importing their constructor.

First verify every interface option from the original vertex-set graphs, using
the established independent checker. Then verify each response and all first
moves, including independent Nim moves, against earlier proven pair facts.
"""
import gzip
import json
from pathlib import Path
import sys

from verify_kayles_boundary import freeze, verify as verify_interfaces


def verify(report):
    assert report["format"] == "node-kayles-response-equivalence-v1"
    base = report["interfaces"]
    checked = verify_interfaces(base, require_equal=False)
    nodes = tuple(map(freeze, base["nodes"]))
    facts, active, used = {}, set(), set()
    for record in report["pairs"]:
        a, b, k = record["a"], record["b"], record["k"]
        assert all(type(v) is int for v in (a, b, k))
        assert 0 <= a < len(nodes) and 0 <= b < len(nodes) and a != b and k >= 0
        assert (a, b, k) not in facts
        facts[a, b, k] = record

    def known(left, right, k, parent=None):
        a, x = left; b, y = right; k ^= x ^ y
        if parent is not None:
            assert (a+b, k) < (parent[0]+parent[1], parent[2]), "nondecreasing proof dependency"
        if a == b:
            assert k == 0, "invalid reflexive leaf"
            return
        key = a, b, k
        assert key in facts, "missing response fact"
        if key in used:
            return
        assert key not in active, "cyclic response proof"
        active.add(key)
        record = facts[key]
        na, nb = nodes[a], nodes[b]
        assert na[0] == nb[0]
        da, db, ma, mb = map(dict, (na[1], nb[1], na[2], nb[2]))
        assert record["deletions"] == record["ports"] == list(na[0])
        for label in na[0]:
            known(da[label], db[label], k, key)
            assert nodes[ma[label][0]][0] == nodes[mb[label][0]][0]
            known(ma[label], mb[label], k, key)
        rows = record["internal"]
        assert [(r["first_side"], r["first_move"]) for r in rows] == [
            (side, i) for side, node in enumerate((na, nb)) for i in range(len(node[3]))]
        for row in rows:
            side, index, response = row["first_side"], row["first_move"], row["response"]
            source, other = (na, nb) if side == 0 else (nb, na)
            other_form = b if side == 0 else a
            removed, (child, extra) = source[3][index]
            post = k ^ extra
            if set(response) == {"nim"}:
                j = response["nim"]
                assert not removed and type(j) is int and 0 <= j < post
                known((child, 0), (other_form, 0), j, key)
            else:
                assert set(response) == {"side", "move"}
                response_side, move = response["side"], response["move"]
                assert type(response_side) is int and response_side in (0, 1)
                assert type(move) is int
                if response_side != side:
                    assert 0 <= move < len(other[3])
                    attack, after = other[3][move]
                    assert attack == removed, "different attacks leave unmatched outside ports"
                    known((child, 0), after, post, key)
                else:
                    assert not removed, "same-side reply cannot restore a deleted port"
                    assert 0 <= move < len(nodes[child][3])
                    attack, after = nodes[child][3][move]
                    assert not attack
                    known(after, (other_form, 0), post, key)
        assert [r["first"] for r in record["nim"]] == list(range(k))
        for row in record["nim"]:
            j, response = row["first"], row["response"]
            assert set(response) == {"side", "move"}
            side, move = response["side"], response["move"]
            assert type(side) is int and side in (0, 1) and type(move) is int
            source, other_form = (na, b) if side == 0 else (nb, a)
            assert 0 <= move < len(source[3])
            attack, after = source[3][move]
            assert not attack, "Nim response must preserve every port"
            known(after, (other_form, 0), j, key)
        active.remove(key); used.add(key)

    roots = report["roots"]
    assert len(roots) == 2 and all(type(i) is int and 0 <= i < len(base["states"]) for i in roots)
    k = report["xor"]
    assert type(k) is int and k >= 0
    known(tuple(base["states"][roots[0]]), tuple(base["states"][roots[1]]), k)
    assert used == set(facts), "unreferenced response facts"
    return dict(interface_nodes=len(nodes), induced_states_checked=checked["induced_states_checked"],
                pair_facts=len(used), xor=k)


if __name__ == "__main__":
    path = Path(sys.argv[1]); raw = path.read_bytes()
    print(json.dumps(verify(json.loads(gzip.decompress(raw) if path.suffix == ".gz" else raw))))
