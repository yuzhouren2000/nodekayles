"""Check saved three-row strategy and middle-exclusion certificates.

Uses only the saved auxiliary values and the proved split/zero-tail lemmas;
does not import the certificate generator or evaluate an unknown long strip.
It checks every middle-move orbit is covered by a valid premise or explicitly
left unresolved. This does not independently recompute the input nimbers.

Run: python3 tests/verify_strip3_analysis.py
"""

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def verify(report, survey):
    assert [x["n"] for x in survey] == list(range(1, len(survey) + 1))
    a = [0] + [x["end_gap"] for x in survey]
    d = [0] + [x["both_gaps"] for x in survey]
    assert report["n_max"] == len(survey)
    assert report["zero_tail_forbidden_gaps"] == [i + 1 for i in range(1, len(d)) if d[i] == 0]
    gaps = report["zero_tail_forbidden_gaps"]
    assert report["impossible_eventual_periods_of_A"] == sorted(
        {p for gap in gaps for p in range(1, gap + 1) if gap % p == 0})

    def check_tail(premise, delta):
        t = premise["zero_tail_length"]
        assert 1 <= t < len(d) and d[t] == 0
        assert premise["gap"] == delta == t + 1

    for certificate in report["middle_strategy_certificates"]:
        p, q = certificate["split_lengths"]
        n = certificate["n"]
        assert 0 <= p <= q < len(a) and n == p + q + 1
        assert certificate["row"] == 2 and certificate["columns"] == [p + 1, q + 1]
        assert a[p] == a[q] == certificate["equal_nimber"]
        assert certificate["beyond_full_survey"] == (n > len(survey))

    no_middle = []
    screens = report["middle_screens"]
    assert [s["n"] for s in screens] == list(range(2, 2 * len(a) + 1, 2))
    for screen in screens:
        n = screen["n"]
        covered = set()

        def cover(p, q):
            assert 0 <= p <= q and p + q == n - 1 and p not in covered
            covered.add(p)

        for proof in screen["exclusions"]:
            p, q = proof["left"], proof["right"]
            cover(p, q)
            if proof["reason"] == "known_unequal":
                assert q < len(a) and a[p] != a[q]
                assert proof["values"] == [a[p], a[q]]
            elif proof["reason"] == "zero_tail_gap":
                check_tail(proof, q - p)
            elif proof["reason"] == "known_equal_witness":
                w = proof["witness"]
                assert 0 <= w < len(a) and p < len(a)
                assert a[p] == a[w] == proof["equal_nimber"]
                check_tail(proof, abs(q - w))
            else:
                raise AssertionError("unknown proof rule")
        for j in screen["winning_columns_in_left_half"]:
            p, q = j - 1, n - j
            cover(p, q)
            assert q < len(a) and a[p] == a[q]
        for p, q in screen["unresolved_split_lengths"]:
            cover(p, q)
        assert covered == set(range(n // 2)), "middle-move orbits incompletely covered"
        expected = "middle_win" if screen["winning_columns_in_left_half"] else (
            "incomplete" if screen["unresolved_split_lengths"] else "no_middle_win")
        assert screen["status"] == expected
        if expected == "no_middle_win" and n > len(survey):
            no_middle.append(n)
    return no_middle


def main():
    survey = [json.loads(s) for s in (ROOT / "data/strip3_survey.jsonl").read_text().splitlines()]
    report = json.loads((ROOT / "data/strip3_analysis.json").read_text())
    no_middle = verify(report, survey)
    print(f"{len(report['middle_strategy_certificates'])} strategy certificates and "
          f"{len(report['middle_screens'])} complete orbit partitions verified")
    print(f"Beyond the full survey, all middle moves excluded for n={no_middle}; outer moves remain unresolved")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
