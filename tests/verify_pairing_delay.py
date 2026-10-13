"""Check the two-round obstruction, including legal moves and full coverage.

No search or SG solver is imported. Neighbourhoods use the independent set
implementation in verify_fixed_pairing; all referenced refutations are checked.

Run: python3 tests/verify_pairing_delay.py data/proofs/pairing_delay_5x8.json.gz
"""

import gzip
import json
import sys
from pathlib import Path

from verify_fixed_pairing import neighbourhoods, verify_result


def verify_obstruction(data):
    if data.get("format") != "nodekayles-pairing-delay-v1" or data.get("status") != "OBSTRUCTION":
        raise ValueError("not a pairing-delay obstruction")
    rows, cols = data["rows"], data["cols"]
    if type(rows) is not int or type(cols) is not int or rows <= 0 or cols <= 0:
        raise ValueError("invalid board dimensions")
    full = frozenset(range(rows * cols))
    nb = neighbourhoods(rows, cols, full)
    first = data["first"]
    if type(first) is not int or first not in full:
        raise ValueError("illegal first attack")
    checked = set()

    def check_refutation(key, vertices):
        result = data["refutations"][key]
        if result["rows"] != rows or result["cols"] != cols or result["status"] != "UNSAT":
            raise ValueError("wrong graph or status in a subproblem")
        if frozenset(result["playable"]) != vertices:
            raise ValueError("refutation concerns a different playable set")
        if key not in checked:
            verify_result(result)
            checked.add(key)

    check_refutation(data["initial_refutation"], full)
    after_first = full - nb[first]
    branches = data["branches"]
    replies = [branch["reply"] for branch in branches]
    if len(set(replies)) != len(replies) or set(replies) != after_first:
        raise ValueError("not all legal first replies are covered")
    for branch in branches:
        reply = branch["reply"]
        after_reply = after_first - nb[reply]
        check_refutation(branch["after_reply_refutation"], after_reply)
        attack = branch["blocking_attack"]
        if type(attack) is not int or attack not in after_reply:
            raise ValueError("illegal blocking attack")
        after_attack = after_reply - nb[attack]
        tails = branch["second_replies"]
        second_replies = [tail["reply"] for tail in tails]
        if len(set(second_replies)) != len(second_replies) or set(second_replies) != after_attack:
            raise ValueError("not all legal second replies are covered")
        for tail in tails:
            remaining = after_attack - nb[tail["reply"]]
            check_refutation(tail["refutation"], remaining)
    return True


def read_obstruction(path):
    path = Path(path)
    encoded = path.read_bytes()
    if path.suffix == ".gz":
        encoded = gzip.decompress(encoded)
    return json.loads(encoded)


def main():
    if not sys.argv[1:]:
        print(__doc__)
        return 2
    for name in sys.argv[1:]:
        data = read_obstruction(name)
        verify_obstruction(data)
        print(f"verified {name}: all {len(data['branches'])} first replies covered")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
