"""Compare compatibility search against raw matching enumeration.

Run: python3 tests/test_fixed_pairing.py
"""

import copy
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from fixed_pairing import rotation_pairing, solve_pairing
from verify_fixed_pairing import neighbourhoods, verify_pairs, verify_result
from pairing_delay import build_obstruction
from verify_pairing_delay import verify_obstruction


def brute_exists(rows, cols, vertices):
    def matchings(remaining):
        if not remaining:
            yield []
            return
        a, *rest = remaining
        for i, b in enumerate(rest):
            for tail in matchings(rest[:i] + rest[i + 1 :]):
                yield [(a, b), *tail]

    for pairs in matchings(vertices):
        try:
            verify_pairs(rows, cols, vertices, pairs)
        except ValueError:
            continue
        return True
    return False


def must_reject(data):
    try:
        verify_result(data)
    except ValueError:
        return
    raise AssertionError("invalid proof accepted")


def main():
    checked = 0
    for rows, cols in [(2, 4), (3, 3)]:
        for mask in range(1 << (rows * cols)):
            result = solve_pairing(rows, cols, mask)
            vertices = [v for v in range(rows * cols) if mask >> v & 1]
            assert (result["status"] == "SAT") == brute_exists(rows, cols, vertices)
            assert verify_result(result)
            checked += 1
    print(f"all {checked} induced positions agree with raw matching enumeration")

    refutation = solve_pairing(5, 8)
    assert refutation["status"] == "UNSAT" and verify_result(refutation)
    truncated = copy.deepcopy(refutation)
    root = truncated["proof"]["nodes"][truncated["proof"]["root"]]
    assert root["branches"]
    root["branches"].pop()
    must_reject(truncated)
    false_parity = copy.deepcopy(refutation)
    false_parity["proof"] = {"kind": "odd-order"}
    must_reject(false_parity)
    print("5x8 refutation verified; omitted branches and false parity rejected")

    limited = solve_pairing(5, 8, node_limit=0)
    assert limited["status"] == "UNKNOWN"
    must_reject(limited)
    timed = solve_pairing(5, 8, timeout=0)
    assert timed["status"] == "UNKNOWN"
    print("search limits return UNKNOWN, never a false nonexistence claim")

    # For odd x even boards the only adjacent rotational pair is in the
    # central row. The six neighbours of this pair, when present, are exactly
    # the first moves whose rotational response gives a pairing certificate.
    cases = 0
    for rows in (1, 3, 5, 7):
        for cols in (2, 4, 6, 8):
            vertices = set(range(rows * cols))
            nb = neighbourhoods(rows, cols, vertices)
            center = (rows // 2) * cols + cols // 2 - 1
            dangerous = {center, center + 1}
            safe_first = (nb[center] | nb[center + 1]) - dangerous
            for first in vertices:
                reply = rows * cols - 1 - first
                legal = reply not in nb[first]
                tail = vertices - nb[first] - nb[reply]
                mask = sum(1 << v for v in tail)
                certificate = rotation_pairing(rows, cols, mask) if legal else None
                assert (certificate is not None) == (first in safe_first)
                if certificate is not None:
                    verify_pairs(rows, cols, tail, certificate)
                cases += 1
    print(f"rotational-response theorem checked on {cases} first moves")

    obstruction = build_obstruction(1, 14, 0)
    assert obstruction["status"] == "OBSTRUCTION" and verify_obstruction(obstruction)
    dropped = copy.deepcopy(obstruction)
    dropped["branches"].pop()
    try:
        verify_obstruction(dropped)
    except ValueError:
        pass
    else:
        raise AssertionError("missing first reply in obstruction accepted")
    assert build_obstruction(2, 4, 0)["status"] == "NO_OBSTRUCTION"
    print("delay-obstruction verification rejects incomplete attack/reply coverage")


if __name__ == "__main__":
    main()
