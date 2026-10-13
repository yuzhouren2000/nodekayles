"""Verify a response survey with independent proofs and the C exact solver.

Run: python3 tests/verify_response_structure.py data/responses_5x8.json.gz
This builds a temporary C wrapper for induced positions up to 64 cells. It
uses the repository's C solver without modifying its source or binary.
"""

import gzip
import json
import subprocess
import sys
import tempfile
from pathlib import Path

from verify_fixed_pairing import neighbourhoods, verify_pairs, verify_result

ROOT = Path(__file__).resolve().parents[1]


def verify_survey(data):
    if data.get("format") != "nodekayles-response-survey-v1":
        raise ValueError("unknown response survey format")
    rows, cols = data["rows"], data["cols"]
    if rows <= 0 or cols <= 0 or rows * cols > 64:
        raise ValueError("this C wrapper supports positive boards up to 64 cells")
    full = set(range(rows * cols))
    nb = neighbourhoods(rows, cols, full)

    def mask(vertices):
        return sum(1 << v for v in vertices)

    expected_values = {mask(full): data["grundy"]}
    covered_first = set()
    checks = 0
    pairing_checks = set()
    for case in data["first_move_orbits"]:
        first = case["first"]
        r, c = divmod(first, cols)
        orbit = {
            rr * cols + cc
            for rr, cc in ((r, c), (rows - 1 - r, c), (r, cols - 1 - c), (rows - 1 - r, cols - 1 - c))
        }
        if rows == cols:
            orbit |= {cc * cols + rr for rr, cc in [divmod(v, cols) for v in orbit]}
        if first not in full or covered_first & orbit:
            raise ValueError("invalid or duplicate first-move orbit")
        covered_first |= orbit
        rest = full - nb[first]
        expected_values[mask(rest)] = case["after_first_grundy"]
        replies = case["replies"]
        if len({t["reply"] for t in replies}) != len(replies) or {t["reply"] for t in replies} != rest:
            raise ValueError("not every legal response is present")
        fixed = {t["reply"]: t["certificate"] for t in case["fixed_pairings"]}
        rotation = {t["reply"]: t["pairs"] for t in case["rotation_replies"]}
        if len(fixed) != len(case["fixed_pairings"]) or len(rotation) != len(case["rotation_replies"]):
            raise ValueError("duplicate prescribed response")
        if not set(fixed) <= rest or not set(rotation) <= rest:
            raise ValueError("illegal prescribed response")
        for entry in replies:
            reply = entry["reply"]
            tail = rest - nb[reply]
            m = mask(tail)
            g = entry["grundy"]
            if m in expected_values and expected_values[m] != g:
                raise ValueError("inconsistent values for the same position")
            expected_values[m] = g
            result = data["pairing_results"][str(m)]
            if result["rows"] != rows or result["cols"] != cols or set(result["playable"]) != tail:
                raise ValueError("pairing result concerns the wrong induced graph")
            if m not in pairing_checks:
                verify_result(result)
                pairing_checks.add(m)
            if (result["status"] == "SAT") != (reply in fixed):
                raise ValueError("fixed-pairing response list is incomplete")
            if reply in fixed:
                verify_result(fixed[reply])
                if fixed[reply]["status"] != "SAT" or set(fixed[reply]["playable"]) != tail:
                    raise ValueError("response certificate uses the wrong graph")
            rotation_valid = all(
                rows * cols - 1 - v in tail
                and rows * cols - 1 - v not in nb[v]
                for v in tail
            )
            if rotation_valid != (reply in rotation):
                raise ValueError("rotational response list is incomplete")
            if reply in rotation:
                verify_pairs(rows, cols, tail, rotation[reply])
                if any(a + b != rows * cols - 1 for a, b in rotation[reply]):
                    raise ValueError("a claimed rotation pair is not rotational")
            checks += 1
    if covered_first != full:
        raise ValueError("first-move orbits do not cover the board")

    # The wrapper calls the independently implemented C SG recursion on
    # arbitrary masks. Temporary output avoids disturbing running campaigns.
    include_path = str(ROOT / "src" / "cfast" / "kayles.c")
    source = (
        '#define main kayles_cli_main\n'
        f'#include "{include_path}"\n'
        '#undef main\n'
        'int main(void) {\n'
        '    int h, w; unsigned long long mask; ht_init(20);\n'
        '    while (scanf("%d %d %llx", &h, &w, &mask) == 3) {\n'
        '        uint64_t position[WORDS] = {0}; position[0] = (uint64_t)mask;\n'
        '        printf("%d\\n", xor_components(h, w, position));\n'
        '    }\n'
        '    return 0;\n'
        '}\n'
    )
    with tempfile.TemporaryDirectory(prefix="kayles-response-check-") as temp:
        temp = Path(temp)
        (temp / "positions.c").write_text(source)
        subprocess.run(["gcc", "-O2", "-o", str(temp / "positions"), str(temp / "positions.c")], check=True)
        positions = sorted(expected_values)
        input_text = "".join(f"{rows} {cols} {m:x}\n" for m in positions)
        output = subprocess.run([str(temp / "positions")], input=input_text, text=True, capture_output=True, check=True, timeout=60)
        values = [int(x) for x in output.stdout.splitlines()]
        if len(values) != len(positions):
            raise ValueError("C solver returned incomplete results")
        for m, g in zip(positions, values):
            if g != expected_values[m]:
                raise ValueError(f"SG mismatch for mask {m}: C={g}, report={expected_values[m]}")
    print(f"{checks} replies covered; {len(pairing_checks)} pairing proofs verified; {len(positions)} C/Python SG values agree")
    return True


def main():
    if not sys.argv[1:]:
        print(__doc__)
        return 2
    for name in sys.argv[1:]:
        raw = Path(name).read_bytes()
        if name.endswith(".gz"):
            raw = gzip.decompress(raw)
        verify_survey(json.loads(raw))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
