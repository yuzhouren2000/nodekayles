#!/usr/bin/env python3
"""Independently check two-move certificates and recompute their SG premises.

Uses coordinate sets, explicit isomorphisms, and the original multiword C
solver.  Imports neither the certificate builder nor the three-row solver.
--geometry-only checks the proof conditional on the recorded SG premises.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]


def full(n):
    return {(r, c) for c in range(1, n + 1) for r in (1, 2, 3)}


def targets(v):
    r, c = v
    return {(r, c), (r - 1, c), (r + 1, c), (r, c - 1), (r, c + 1)}


def family(f, n, j):
    if type(n) is not int or not 0 <= n <= 128:
        raise ValueError("invalid family width")
    if f in ("A", "D", "T", "C") and j is not None:
        raise ValueError("one-parameter family has a port index")
    cells = full(n)
    if f == "A":
        return cells - {(2, n)}
    if f == "D":
        return cells - {(2, 1), (2, n)}
    if f == "T":
        return cells - {(1, n)}
    if f == "C":
        return cells - {(1, n), (3, n)}
    if f == "point" and n == 1 and j is None:
        return {(2, 1)}
    if f in ("B", "E", "H") and type(j) is int and 1 <= j <= n:
        cells -= {(1, j - 1), (1, j), (1, j + 1), (2, j)}
        if f in ("E", "H"):
            cells.discard((2, 1))
        if f == "H":
            cells.discard((2, n))
        return cells
    raise ValueError(f"invalid family descriptor {f},{n},{j}")


def encode(cells):
    return sum(1 << (3 * (c - 1) + r - 1) for r, c in cells)


def decode(s):
    mask = int(s, 16)
    if mask < 0:
        raise ValueError("negative mask")
    return {(v % 3 + 1, v // 3 + 1) for v in range(mask.bit_length())
            if mask >> v & 1}


def transform(cells, n, horizontal, vertical, offset=0):
    return {(4 - r if vertical else r,
             offset + (n + 1 - c if horizontal else c)) for r, c in cells}


def check(report):
    if report.get("format") != "nodekayles-strip3-two-move-v1":
        raise ValueError("unrecognised certificate")
    n = report["n"]
    if type(n) is not int or n < 2 or n % 2:
        raise ValueError("expected a positive even strip width")
    board = full(n)
    premises = {}
    for p in report["premises"]:
        key = p["family"], p["n"], p["j"]
        if key in premises or type(p["value"]) is not int or p["value"] < 0:
            raise ValueError("duplicate or invalid premise")
        cells = family(*key)
        if decode(p["mask"]) != cells:
            raise ValueError("premise geometry mismatch")
        if p["family"] == "point" and p["value"] != 1:
            raise ValueError("an isolated point has SG value 1")
        premises[key] = p
    representatives = {(r, c) for r in (1, 2) for c in range(1, n // 2 + 1)}
    checked, covered, used_keys = set(), set(), set()
    for case in report["cases"]:
        first, reply = tuple(case["first"]), tuple(case["reply"])
        if first not in representatives or first in checked:
            raise ValueError("invalid or duplicate orbit representative")
        checked.add(first)
        child = board - targets(first)
        if reply not in child:
            raise ValueError("illegal response")
        tail = child - targets(reply)
        if decode(case["tail_mask"]) != tail:
            raise ValueError("wrong two-move successor")
        if case["rule"] == "safe_rotation":
            if case["pieces"]:
                raise ValueError("rotation proof must not contain unchecked pieces")
            for v in tail:
                partner = (4 - v[0], n + 1 - v[1])
                if partner not in tail or partner in targets(v):
                    raise ValueError("unsafe rotational pairing")
        elif case["rule"] in ("middle_middle", "middle_outer", "outer_middle",
                              "opposite_outer"):
            used, pieces, value = set(), [], 0
            for p in case["pieces"]:
                key = p["family"], p["n"], p["j"]
                if key not in premises or p["value"] != premises[key]["value"]:
                    raise ValueError("piece lacks a matching SG premise")
                used_keys.add(key)
                mapping = p["isomorphism"]
                if type(mapping["horizontal"]) is not bool or \
                        type(mapping["vertical"]) is not bool or \
                        type(mapping["column_offset"]) is not int:
                    raise ValueError("invalid isomorphism")
                cells = decode(p["mask"])
                mapped = transform(family(*key), p["n"], mapping["horizontal"],
                                   mapping["vertical"], mapping["column_offset"])
                if cells != mapped or not cells or cells & used or not cells <= tail:
                    raise ValueError("false piece isomorphism or overlapping pieces")
                used |= cells
                pieces.append(cells)
                value ^= p["value"]
            if used != tail or value != 0:
                raise ValueError("pieces fail to cover a zero-valued successor")
            for k, piece in enumerate(pieces):
                for other in pieces[k + 1:]:
                    if any(targets(v) & other for v in piece):
                        raise ValueError("edges cross the claimed disjoint sum")
        else:
            raise ValueError("unknown proof rule")
        for horizontal in (False, True):
            for vertical in (False, True):
                covered |= transform({first}, n, horizontal, vertical)
    unresolved = [tuple(v) for v in report["unresolved"]]
    if len(unresolved) != len(set(unresolved)) or set(unresolved) != representatives - checked:
        raise ValueError("unresolved list does not match orbit coverage")
    status = "PARTIAL" if unresolved else "P_CERTIFICATE"
    if report["status"] != status:
        raise ValueError("a partial cover cannot prove a P-position")
    if not unresolved and covered != board:
        raise ValueError("not every legal first move is covered")
    if used_keys != set(premises):
        raise ValueError("unused or missing premises")
    if report["max_premise_width"] != max((p["n"] for p in premises.values()), default=0):
        raise ValueError("incorrect width bound")
    return sorted(premises.values(), key=lambda p: (p["n"], p["family"], p["j"] or 0)), len(covered)


def recompute(premises, timeout):
    source = '''#define main original_cli_main
#include "%s"
#undef main
int main(void) {
    const char *limit = getenv("KAYLES_HT_MAXBITS");
    if (limit) ht_maxbits = (size_t)atoi(limit);
    int n; unsigned long long words[WORDS]; ht_init(18);
    while (scanf("%%d", &n) == 1) {
        if (n < 1 || 3*n > MAXBITS) return 2;
        for (int w = 0; w < WORDS; w++)
            if (scanf("%%llx", &words[w]) != 1) return 2;
        uint64_t mask[WORDS] = {0};
        for (int c = 0; c < n; c++) for (int r = 0; r < 3; r++) {
            int v = 3*c+r, u = r*n+c;
            if ((words[v/64] >> (v %% 64)) & 1) mask[u/64] |= 1ULL << (u %% 64);
        }
        printf("%%d\\n", xor_components(3, n, mask)); fflush(stdout);
    }
    return 0;
}
''' % (ROOT / "src/cfast/kayles.c")
    selected = [p for p in premises if p["family"] != "point"]
    inputs = "".join(str(p["n"]) + " " + " ".join(
        f"{(int(p['mask'], 16) >> (64*w)) & ((1 << 64) - 1):x}"
        for w in range(6)) + "\n" for p in selected)
    with tempfile.TemporaryDirectory(prefix="strip3-two-move-check-") as folder:
        folder = Path(folder)
        (folder / "original.c").write_text(source)
        subprocess.run(["gcc", "-O3", "-march=native", "-o", str(folder / "original"),
                        str(folder / "original.c")], check=True)
        started = time.monotonic()
        result = subprocess.run([str(folder / "original")], input=inputs, text=True,
                                capture_output=True, timeout=timeout,
                                env={**os.environ, "KAYLES_HT_MAXBITS": "26"})
        if result.returncode:
            raise ValueError(result.stderr)
        values = list(map(int, result.stdout.split()))
        if values != [p["value"] for p in selected]:
            raise ValueError(f"independent SG recomputation failed: {values}")
        return {"solver": "original multiword C", "positions": len(selected),
                "source_sha256": hashlib.sha256(
                    (ROOT / "src/cfast/kayles.c").read_bytes()).hexdigest(),
                "elapsed_s": round(time.monotonic() - started, 3)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--geometry-only", action="store_true")
    parser.add_argument("--timeout", type=float, default=180)
    args = parser.parse_args()
    report = json.loads(args.certificate.read_text())
    premises, covered = check(report)
    print(f"geometry: {len(report['cases'])} orbit responses, {covered}/{3 * report['n']} "
          f"first moves, {len(premises)} SG premises, width <= {report['max_premise_width']}")
    if args.geometry_only:
        print("SG premises not recomputed; this is a conditional proof check")
    else:
        print("independent recomputation:", recompute(premises, args.timeout))
    if report["status"] == "P_CERTIFICATE":
        prefix = "conditional conclusion" if args.geometry_only else "conclusion"
        print(f"{prefix}: G(3x{report['n']}) = 0")
    else:
        print("PARTIAL: board outcome unresolved; remaining representatives:", report["unresolved"])


if __name__ == "__main__":
    main()
