#!/usr/bin/env python3
"""Port families B, E, H on the 3 x n strip: geometry, survey, closure, target.

Bit encoding is column-major, vertex (row r in 1..3, column c >= 1) is bit
3*(c-1)+(r-1).  All geometry below is built from one-based coordinates and then
encoded, so the encoding itself can be cross-checked against a second
coordinate-based implementation (tests/test_strip3_ports.py).

Families, with F_n the full strip:
  A(n)   = F_n minus the right end middle cell (2,n)
  D(n)   = F_n minus (2,1) and (2,n);  D(1) removes one cell twice
  S(n,j) = {(1,j-1),(1,j),(1,j+1),(2,j)} intersected with F_n, 1 <= j <= n
  B(n,j) = F_n minus S(n,j)                       (top-row first move at j)
  E(n,j) = B(n,j) minus (2,1)                     (left end middle gap)
  H(n,j) = B(n,j) minus (2,1) and (2,n)           (both end middle gaps)

Masks are passed to src/cfast/strip3 via its `masks` mode (decimal n and hex
column-major mask per line).  The Python solver in this file is a slow
reference used for small positions and structural exploration only.

Subcommands:
  survey   --max-n N --ht-bits B --timeout S --output FILE
  closure  --max-n N --output FILE
  target   --n N --columns J,... --ht-bits B --timeout S --output FILE
  recompute --max-n N            (reference solver values, for tests/spot checks)
"""

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BINARY = ROOT / "src/cfast/strip3"
BASE_SURVEY = ROOT / "data/strip3_survey.jsonl"
RUN_LOG = ROOT / "data/strip3_ports_run.json"
VERSION = 1

# ------------------------------------------------------------------ geometry


def cells_F(n):
    return {(r, c) for r in (1, 2, 3) for c in range(1, n + 1)}


def cells_S(n, j):
    return {(r, c) for r, c in ((1, j - 1), (1, j), (1, j + 1), (2, j))
            if 1 <= c <= n}


def family_cells(name, n, j=None):
    """One-based cell sets.  n=0 is the empty position for F/A/D."""
    F = cells_F(n)
    if name == "F":
        return F
    if name == "A":
        return F - {(2, n)} if n >= 1 else set()
    if name == "D":
        return F - {(2, 1), (2, n)} if n >= 1 else set()
    if name in ("B", "E", "H"):
        base = F - cells_S(n, j)
        if name == "B":
            return base
        base = base - {(2, 1)}
        if name == "E":
            return base
        return base - {(2, n)}
    raise ValueError(name)


def to_mask(cells):
    mask = 0
    for r, c in cells:
        mask |= 1 << (3 * (c - 1) + (r - 1))
    return mask


def mask_to_cells(mask):
    cells = set()
    while mask:
        low = mask & -mask
        v = low.bit_length() - 1
        mask ^= low
        cells.add((v % 3 + 1, v // 3 + 1))
    return cells


def family_mask(name, n, j=None):
    return to_mask(family_cells(name, n, j))


def closed_cells(cell):
    r, c = cell
    ur, uc = r, c
    return {(rr, cc) for rr, cc in
            ((r, c), (r - 1, c), (r + 1, c), (r, c - 1), (r, c + 1))
            if 1 <= rr <= 3 and cc >= 1}


def legal_moves(cells):
    """All one-based cells at which a move is legal, sorted column-major."""
    return sorted(cells, key=lambda cell: (cell[1], cell[0]))


def successor_cells(cells, cell):
    return cells - closed_cells(cell)


def cell_to_display(cell):
    return [cell[0], cell[1]]


ROWMASK = []
for _r in range(3):
    _m = 0
    for _c in range(1, 43):
        _m |= 1 << (3 * (_c - 1) + _r)
    ROWMASK.append(_m)


def neighbours_mask(mask):
    """Orthogonal neighbours of a column-major mask; same shifts as strip3.c.
    Out-of-board bits cannot reappear inside the position mask."""
    return ((mask << 3) | (mask >> 3)
            | ((mask & (ROWMASK[0] | ROWMASK[1])) << 1)
            | ((mask & (ROWMASK[1] | ROWMASK[2])) >> 1))


def closed_mask(mask, bit):
    return mask & ~(bit | neighbours_mask(bit))


def split_components(mask):
    """Connected components of a mask, largest width last, deterministic."""
    out = []
    while mask:
        low = mask & -mask
        comp = low
        frontier = low
        while frontier:
            frontier = neighbours_mask(frontier) & mask & ~comp
            comp |= frontier
        mask &= ~comp
        out.append(comp)
    return sorted(out, key=lambda m: (m.bit_length(), m))


# --------------------------------------------------------- reference solver

_reference_memo = {}


def ref_value(mask):
    """Independent slow SG: coordinates -> closed neighbourhoods -> mex."""
    memo = _reference_memo
    if mask == 0:
        return 0
    if mask in memo:
        return memo[mask]
    components = split_components(mask)
    if len(components) > 1:
        g = 0
        for comp in components:
            g ^= ref_value(comp)
        memo[mask] = g
        return g
    seen = set()
    x = mask
    while x:
        bit = x & -x
        x ^= bit
        seen.add(ref_value(closed_mask(mask, bit)))
    g = 0
    while g in seen:
        g += 1
    memo[mask] = g
    return g


def ref_value_coordinates(cells):
    return ref_value(to_mask(cells))


# ------------------------------------------------------------- C solver I/O


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def run_masks(positions, binary=DEFAULT_BINARY, ht_bits=26, timeout=300,
              reductions=None, precision="exact"):
    """Run one strip3 `masks` batch.  Returns (results, info).

    results[i] is None (UNKNOWN: budget/timeout before its turn) or the exact
    nimber.  info records return code, timing and stderr statistics.
    """
    env = dict(os.environ)
    env["STRIP3_HT_BITS"] = str(ht_bits)
    if reductions is not None:
        env["STRIP3_REDUCTIONS"] = str(reductions)
    text = "".join(f"{n} {mask:x}\n" for n, mask in positions)
    started = time.monotonic()
    status = "ok"
    stderr = ""
    try:
        proc = subprocess.run([str(binary), "masks"], input=text, text=True,
                              capture_output=True, timeout=timeout, env=env)
        code = proc.returncode
        stderr = proc.stderr
        lines = proc.stdout.splitlines()
    except subprocess.TimeoutExpired as exc:
        code = None
        status = "timeout"
        stderr = (exc.stderr or b"").decode() if isinstance(exc.stderr, bytes) else (exc.stderr or "")
        lines = (exc.stdout or b"").decode().splitlines() if isinstance(exc.stdout, bytes) else (exc.stdout or "").splitlines()
    elapsed = time.monotonic() - started
    if code not in (0, None) and status == "ok":
        status = "exit_%s" % code
    results = [None] * len(positions)
    for i, line in enumerate(lines[:len(positions)]):
        results[i] = int(line)
    info = {"returncode": code, "status": status, "elapsed_s": round(elapsed, 3),
            "emitted": len(lines), "requested": len(positions),
            "ht_bits": ht_bits, "reductions": reductions,
            "binary_sha256": sha256_file(binary) if binary.exists() else None,
            "stderr_tail": stderr.strip().splitlines()[-3:] if stderr else []}
    return results, info


# ------------------------------------------------------- template matching

def _shift_min_column(cells):
    if not cells:
        return cells
    shift = min(c for _r, c in cells) - 1
    return {(r, c - shift) for r, c in cells}


def _flip_vertical(cells):
    return {(4 - r, c) for r, c in cells}


def _flip_horizontal(cells, width):
    return {(r, width + 1 - c) for r, c in cells}


def _is_path_or_cycle(cells):
    if not cells:
        return None
    degree = {}
    for r, c in cells:
        degree[(r, c)] = sum(1 for nb in closed_cells((r, c)) - {(r, c)}
                             if nb in cells)
    if max(degree.values()) <= 2:
        total = sum(degree.values())
        if any(d < 2 for d in degree.values()):
            return "path(%d)" % len(cells)
        return "cycle(%d)" % len(cells)
    return None


def identify_piece(cells):
    """Label a piece by declared geometric templates up to translation and the
    two board reflections.  This vocabulary claim says nothing about game
    equivalence beyond exact geometry."""
    cells = set(cells)
    if not cells:
        return {"label": "empty", "width": 0, "reflection": "identity"}
    if not (1 <= min(r for r, _ in cells) and max(r for r, _ in cells) <= 3):
        return {"label": "off_board", "width": None}
    normal = _shift_min_column(cells)
    candidates = []
    for name, transform in (("identity", lambda s: s),
                            ("vertical", _flip_vertical),
                            ("horizontal", lambda s: _flip_horizontal(s, max(c for _r, c in s))),
                            ("both", lambda s: _flip_horizontal(_flip_vertical(s), max(c for _r, c in s)))):
        t = _shift_min_column(transform(normal))
        width = max(c for _r, c in t)
        if t == cells_F(width):
            return {"label": "F(%d)" % width, "width": width, "reflection": name}
        for small in ("A", "D"):
            if t == family_cells(small, width):
                return {"label": "%s(%d)" % (small, width), "width": width,
                        "reflection": name}
        for j in range(1, width + 1):
            for small in ("B", "E", "H"):
                if t == family_cells(small, width, j):
                    return {"label": "%s(%d,%d)" % (small, width, j),
                            "width": width, "reflection": name}
        candidates.append(t)
    shape = _is_path_or_cycle(cells)
    if shape:
        return {"label": shape, "width": max(c for _r, c in cells), "reflection": "none"}
    if len(cells) == 4 and len(split_components(to_mask(cells))) == 1:
        degree = max(sum(1 for nb in closed_cells(v) - {v} if nb in cells)
                     for v in cells)
        if degree == 3:
            return {"label": "T4(claw)", "width": max(c for _r, c in cells),
                    "reflection": "none"}
    return {"label": "unknown", "width": max(c for _r, c in cells),
            "cells": sorted(cells),
            "normalised": sorted(candidates[0]) if candidates else []}


def describe_position(mask):
    """Components with labels, for reporting only."""
    pieces = []
    for comp in split_components(mask):
        cells = mask_to_cells(comp)
        info = identify_piece(cells)
        info["mask"] = "%x" % comp
        info["cells"] = sorted(cells)
        pieces.append(info)
    return pieces


# ------------------------------------------------------------------ survey

def collect_family_positions(max_n, families=("B", "E", "H")):
    rows = []
    for n in range(1, max_n + 1):
        for j in range(1, n + 1):
            for fam in families:
                rows.append({"family": fam, "n": n, "j": j,
                             "mask": family_mask(fam, n, j)})
    return rows


def survey(args):
    base = collect_family_positions(args.max_n)
    solve = {}  # mask -> n

    def add(mask, n):
        solve[mask] = n

    for row in base:
        add(row["mask"], row["n"])
    level1 = {}
    for row in base:
        cells = mask_to_cells(row["mask"])
        level1[row["mask"]] = []
        for cell in legal_moves(cells):
            child = to_mask(successor_cells(cells, cell))
            level1[row["mask"]].append((cell, child))
            add(child, row["n"])
    # B(10,4) response certificates: solve the options of every option.
    responses = {}
    focus = family_mask("B", 10, 4)
    if focus in {row["mask"] for row in base}:
        for _cell, child in level1[focus]:
            sub = []
            for cell2 in legal_moves(mask_to_cells(child)):
                grand = to_mask(successor_cells(mask_to_cells(child), cell2))
                sub.append((cell2, grand))
                add(grand, 10)
            responses[child] = sub
    positions = [(solve[m], m) for m in solve]
    results, info = run_masks(positions, binary=args.binary, ht_bits=args.ht_bits,
                              timeout=args.timeout, reductions=args.reductions)
    value = {mask: results[i] for i, (n, mask) in enumerate(positions)}

    old_survey = {}
    if BASE_SURVEY.exists():
        for line in BASE_SURVEY.read_text().splitlines():
            if line.strip():
                rec = json.loads(line)
                old_survey[rec["n"]] = rec

    with open(args.output, "w") as handle:
        handle.write(json.dumps({
            "version": VERSION, "kind": "strip3_ports_survey",
            "max_n": args.max_n, "families": ["B", "E", "H"],
            "ht_bits": args.ht_bits, "reductions": args.reductions,
            "strip3_sha256": sha256_file(ROOT / "src/cfast/strip3.c"),
            "focus": "B(10,4) options and one response level"}) + "\n")
        cross_checked = 0
        cross_mismatch = []
        for row in base:
            n, fam, j, mask = row["n"], row["family"], row["j"], row["mask"]
            v = value.get(mask)
            options = []
            for cell, child in level1[mask]:
                entry = {"cell": cell_to_display(cell), "mask": "%x" % child,
                         "value": value.get(child),
                         "pieces": describe_position(child)}
                if fam == "B" and (n, j) == (10, 4) and child in responses:
                    entry["zero_replies"] = [
                        {"cell": cell_to_display(cell2), "mask": "%x" % grand,
                         "value": value.get(grand),
                         "pieces": describe_position(grand)}
                        for cell2, grand in responses[child] if value.get(grand) == 0]
                options.append(entry)
            record = {"family": fam, "n": n, "j": j, "mask": "%x" % mask,
                      "value": v, "options": options}
            if fam == "B" and n in old_survey and 3 * (j - 1) < len(old_survey[n]["options"]):
                known = old_survey[n]["options"][3 * (j - 1)]
                record["known_b_from_strip3_survey"] = known
                if known == v:
                    cross_checked += 1
                else:
                    cross_mismatch.append([n, j, known, v])
            handle.write(json.dumps(record) + "\n")
    run_info = {"kind": "survey", "output": str(args.output), "max_n": args.max_n,
                "positions": len(positions), "cross_checked_b": cross_checked,
                "cross_mismatch": cross_mismatch, "c_info": info}
    append_run(run_info)
    if cross_mismatch:
        print("MISMATCH against strip3_survey options:", cross_mismatch[:5])
        return 1
    print(f"survey: {len(base)} positions, {len(positions)} solved, "
          f"{cross_checked} b-values cross-checked against the full survey, "
          f"elapsed {info['elapsed_s']}s, status {info['status']}")
    return 0


# ----------------------------------------------------------------- closure

def closure(args):
    """Every legal move of every B/E/H position up to max_n, plus B(10,4):
    split the successor at empty columns (the only legal cuts), then check for
    every run of consecutive non-empty columns whether it matches a declared
    family position (F/A/D/B/E/H under the board symmetries) or is a path or
    cycle.  Minimal unrepresented runs are recorded with full context."""
    report = {"version": VERSION, "kind": "strip3_ports_closure",
              "max_n": args.max_n, "vocabulary":
              ["F", "A", "D", "B", "E", "H", "path", "cycle"],
              "vocabulary_extended": ["F", "A", "D", "B", "E", "H", "path",
                                      "cycle", "T4(claw)"],
              "symmetries": ["column translation", "horizontal reflection",
                             "vertical reflection (top-bottom)"],
              "note": "template representation of whole column runs; equal "
                      "nimbers are not treated as matches; T4(claw) counts as "
                      "unrepresented in the base vocabulary and represented in "
                      "the extended one",
              "unrepresented": [], "extended_unrepresented": [],
              "counts": {}, "positions": 0, "moves": 0,
              "moves_by_type": {"middle": 0, "outer": 0}}
    unrepresented = {}
    extended = {}

    def runs_of(cells):
        if not cells:
            return []
        columns = sorted({c for _r, c in cells})
        runs, start, prev = [], columns[0], columns[0]
        for c in columns[1:]:
            if c != prev + 1:
                runs.append((start, prev))
                start = c
            prev = c
        runs.append((start, prev))
        return runs

    unrep_count = 0
    extended_count = 0

    def audit_one(cells, origin):
        nonlocal unrep_count, extended_count
        for cell in legal_moves(cells):
            kind = "middle" if cell[0] == 2 else "outer"
            report["moves"] += 1
            report["moves_by_type"][kind] += 1
            child = successor_cells(cells, cell)
            if not child:
                report["counts"]["empty"] = report["counts"].get("empty", 0) + 1
                continue
            ok = True
            ok_extended = True
            pieces = []
            for start, end in runs_of(child):
                run_cells = {(r, c) for r, c in child if start <= c <= end}
                label = identify_piece(run_cells)
                pieces.append({"columns": [start, end], "label": label["label"],
                               "reflection": label.get("reflection"),
                               "size": len(run_cells),
                               "cells": sorted(run_cells)})
                key = label["label"] if label["label"] != "unknown" else "unknown"
                report["counts"][key] = report["counts"].get(key, 0) + 1
                if key == "unknown":
                    ok = False
                    ok_extended = False
                elif key == "T4(claw)":
                    ok = False
            if not ok:
                unrep_count += 1
                signature = (len(child), "%x" % to_mask(child))
                if signature not in unrepresented:
                    unrepresented[signature] = {
                        "successor_mask": "%x" % to_mask(child),
                        "size": len(child),
                        "occurrences": 0,
                        "origin": dict(origin), "move": cell_to_display(cell),
                        "move_type": kind, "pieces": pieces}
                unrepresented[signature]["occurrences"] += 1
            if not ok_extended:
                extended_count += 1
                signature = (len(child), "%x" % to_mask(child))
                if signature not in extended:
                    extended[signature] = {
                        "successor_mask": "%x" % to_mask(child),
                        "size": len(child),
                        "occurrences": 0,
                        "origin": dict(origin), "move": cell_to_display(cell),
                        "move_type": kind, "pieces": pieces}
                extended[signature]["occurrences"] += 1

    for n in range(1, args.max_n + 1):
        for j in range(1, n + 1):
            for fam in ("B", "E", "H"):
                report["positions"] += 1
                audit_one(family_cells(fam, n, j), {"family": fam, "n": n, "j": j})
    report["positions"] += 1
    audit_one(family_cells("B", 10, 4), {"family": "B", "n": 10, "j": 4})
    report["unrepresented"] = [unrepresented[k] for k in sorted(unrepresented)][:100]
    report["unrepresented_count"] = unrep_count
    report["extended_unrepresented"] = [extended[k] for k in sorted(extended)][:100]
    report["extended_unrepresented_count"] = extended_count
    Path(args.output).write_text(json.dumps(report, indent=2) + "\n")
    append_run({"kind": "closure", "output": str(args.output),
                "max_n": args.max_n, "positions": report["positions"],
                "moves": report["moves"],
                "unrepresented_successors": unrep_count,
                "distinct_unrepresented": len(unrepresented),
                "extended_unrepresented_successors": extended_count,
                "distinct_extended_unrepresented": len(extended)})
    print(f"closure: {report['positions']} positions, {report['moves']} moves "
          f"({report['moves_by_type']['middle']} middle, "
          f"{report['moves_by_type']['outer']} outer), "
          f"{unrep_count} unrepresented successors "
          f"({len(unrepresented)} distinct); with T4(claw) added, "
          f"{extended_count} remain ({len(extended)} distinct)")
    return 0


# ------------------------------------------------------------------ target

def target(args):
    columns = [int(c) for c in args.columns.split(",")]
    n = args.n
    middle_reference = None
    analysis = ROOT / "data/strip3_analysis.json"
    if analysis.exists():
        screens = {s["n"]: s for s in json.loads(analysis.read_text())["middle_screens"]}
        if n in screens:
            middle_reference = {"source": str(analysis.relative_to(ROOT)),
                                "status": screens[n]["status"],
                                "exclusions": len(screens[n]["exclusions"]),
                                "unresolved": screens[n]["unresolved_split_lengths"],
                                "note": "all middle first moves were excluded by "
                                        "zero-tail / known-value proofs; the full "
                                        "board outcome still needs every edge-row "
                                        "orbit"}
    out = {"version": VERSION, "kind": "strip3_target",
           "n": n, "columns": columns, "ht_bits": args.ht_bits,
           "timeout_s": args.timeout, "reductions": args.reductions,
           "middle_screens_reference": middle_reference,
           "note": "each column is one strip3 masks process; UNKNOWN means the "
                   "bounded memo or the timeout hit before a value was emitted; "
                   "a zero value proves n=N only if every edge orbit is covered, "
                   "and is only a P-certificate with all option values nonzero",
           "results": []}
    for j in columns:
        mask = family_mask("B", n, j)
        positions = [(n, mask)]
        results, info = run_masks(positions, binary=args.binary,
                                  ht_bits=args.ht_bits, timeout=args.timeout,
                                  reductions=args.reductions)
        entry = {"column": j, "mask": "%x" % mask, "value": results[0],
                 "children": [], "c_info": info}
        if results[0] == 0:
            # P-position certificate: every legal option must be non-zero.
            cells = mask_to_cells(mask)
            child_positions = [(n, to_mask(successor_cells(cells, cell)))
                               for cell in legal_moves(cells)]
            child_values, child_info = run_masks(
                child_positions, binary=args.binary, ht_bits=args.ht_bits,
                timeout=args.timeout, reductions=args.reductions)
            entry["children"] = [{"cell": cell_to_display(cell),
                                  "mask": "%x" % child[1],
                                  "value": child_values[i]}
                                 for i, (cell, child) in enumerate(
                                     zip(legal_moves(cells), child_positions))]
            entry["children_info"] = child_info
            entry["p_certificate"] = all(v is not None and v != 0
                                         for v in child_values)
            entry["untested_children"] = sum(v is None for v in child_values)
        out["results"].append(entry)
        Path(args.output).write_text(json.dumps(out, indent=2) + "\n")
        v = entry["value"]
        print(f"B({n},{j}) = {'UNKNOWN' if v is None else v} "
              f"[{info['status']}, {info['elapsed_s']}s]", flush=True)
    append_run({"kind": "target", "output": str(args.output), "n": n,
                "columns": columns, "results": [
                    {"column": e["column"], "value": e["value"],
                     "status": e["c_info"]["status"],
                     "elapsed_s": e["c_info"]["elapsed_s"]}
                    for e in out["results"]]})
    return 0


# --------------------------------------------------------------- recompute

def recompute(args):
    values = {}
    for n in range(0, args.max_n + 1):
        for fam in ("F", "A", "D"):
            if n == 0 or fam == "F":
                values["%s(%d)" % (fam, n)] = ref_value(family_mask(fam, n))
        for j in range(1, n + 1):
            for fam in ("B", "E", "H"):
                values["%s(%d,%d)" % (fam, n, j)] = ref_value(
                    family_mask(fam, n, j))
    print(json.dumps(values))
    return 0


# ----------------------------------------------------------------- collect

def collect(args):
    """Assemble a target result file from raw per-column run outputs produced
    elsewhere (value/log files plus summary_*.txt lines).  Column identity is
    recovered by matching the recorded masks against freshly built ones."""
    import re
    n = args.n
    want = {("%x" % family_mask("B", n, j)): j for j in range(1, n + 1)}
    out = {"version": VERSION, "kind": "strip3_target", "n": n,
           "ht_bits": args.ht_bits, "timeout_s": args.timeout,
           "machine": args.machine, "strip3_sha256": args.strip3_sha256,
           "middle_screens_reference": None,
           "note": "assembled from raw per-column strip3 masks runs; UNKNOWN "
                   "means budget or timeout before a value was emitted; a zero "
                   "value proves the full board is an N-position only if every "
                   "edge orbit is covered and all options are nonzero",
           "results": []}
    analysis = ROOT / "data/strip3_analysis.json"
    if analysis.exists():
        screens = {s["n"]: s for s in json.loads(analysis.read_text())["middle_screens"]}
        if n in screens:
            out["middle_screens_reference"] = {
                "source": str(analysis.relative_to(ROOT)),
                "status": screens[n]["status"],
                "exclusions": len(screens[n]["exclusions"])}
    pattern = re.compile(r"seq=(\d+) mask=([0-9a-f]+) rc=(\d+) value=(\S+) seconds=(\d+)")
    seen = {}
    for summary in sorted(Path(args.source_dir).glob("summary_*.txt")):
        for line in summary.read_text().splitlines():
            m = pattern.match(line)
            if not m:
                continue
            seq, mask, rc, value, seconds = m.groups()
            j = want.get(mask)
            if j is None:
                continue
            stem = summary.stem.replace("summary_", "")
            log = Path(args.source_dir) / f"{stem}_{seq}.log"
            value_file = Path(args.source_dir) / f"{stem}_{seq}.value"
            raw = value_file.read_text().strip() if value_file.exists() else ""
            entry = {"column": j, "mask": mask,
                     "value": int(value) if value.isdigit() else None,
                     "returncode": int(rc), "seconds": int(seconds),
                     "status": "ok" if rc == "0" else
                               "budget_exhausted" if rc == "3" else
                               "timeout" if rc == "124" else f"exit_{rc}",
                     "raw_output": raw,
                     "stderr_tail": log.read_text().strip().splitlines()[-2:]
                     if log.exists() else []}
            seen[j] = entry
    out["results"] = [seen[j] for j in sorted(seen)]
    out["covers_columns"] = sorted(seen)
    orbits = range(1, (n + 1) // 2 + 1)
    out["complete_edge_orbit_cover"] = all(
        (j in seen) or ((n + 1 - j) in seen) for j in orbits)
    Path(args.output).write_text(json.dumps(out, indent=2) + "\n")
    append_run({"kind": "target_collect", "output": str(args.output), "n": n,
                "source_dir": str(args.source_dir), "columns": sorted(seen),
                "machine": args.machine, "strip3_sha256": args.strip3_sha256,
                "ht_bits": args.ht_bits, "timeout_s": args.timeout,
                "statuses": {j: seen[j]["status"] for j in sorted(seen)},
                "values": {j: seen[j]["value"] for j in sorted(seen)}})
    for j in sorted(seen):
        v = seen[j]["value"]
        print(f"B({n},{j}) = {'UNKNOWN' if v is None else v} "
              f"[{seen[j]['status']}, {seen[j]['seconds']}s]")
    return 0


# ------------------------------------------------------------------ value

def position_value(args):
    """Exact value of a single position, parallelised over its move children.

    All children are evaluated in one shared-memo process (the memo is shared
    across the children, so overlapping subtrees are computed once); the mex of
    the children's values is the position's Grundy value.  This gives a large
    job to many threads even for one board, which the coarse per-line
    parallelism cannot do.
    """
    n = args.n
    mask = int(args.mask, 16)
    if mask & ~((1 << (3 * n)) - 1):
        raise ValueError("mask does not fit the board")
    cells = mask_to_cells(mask)
    children = []
    for cell in legal_moves(cells):
        children.append((n, to_mask(successor_cells(cells, cell))))
    if not children:
        raise ValueError("position has no legal move")
    env = dict(os.environ)
    env["STRIP3_HT_BITS"] = str(args.ht_bits)
    if args.load is not None:
        env["STRIP3_HT_LOAD"] = str(args.load)
    if args.prealloc:
        env["STRIP3_HT_PREALLOC"] = "1"
    if args.threads:
        env["STRIP3_THREADS"] = str(args.threads)
    text = "".join(f"{c[0]} {c[1]:x}\n" for c in children)
    started = time.monotonic()
    proc = subprocess.run([str(args.binary), "masks"], input=text, text=True,
                          capture_output=True, timeout=args.timeout, env=env)
    elapsed = time.monotonic() - started
    if proc.returncode != 0:
        raise RuntimeError(f"solver exited {proc.returncode}: {proc.stderr[-400:]}")
    values = [int(x) for x in proc.stdout.split()]
    if len(values) != len(children):
        raise RuntimeError(f"expected {len(children)} child values, got {len(values)}")
    present = set(values)
    value = 0
    while value in present:
        value += 1
    out = {"version": VERSION, "kind": "strip3_value", "n": n,
           "mask": args.mask.lower().lstrip("0") or "0",
           "family": args.family, "value": value,
           "children": values, "child_cells": [cell_to_display(c) for c in legal_moves(cells)],
           "ht_bits": args.ht_bits, "threads": args.threads,
           "load": args.load, "prealloc": args.prealloc,
           "elapsed_s": round(elapsed, 1),
           "binary_sha256": sha256_file(args.binary) if args.binary.exists() else None,
           "note": "all legal moves evaluated; the value is their mex"}
    Path(args.output).write_text(json.dumps(out, indent=1) + "\n")
    append_run({"kind": "value", "output": str(args.output), "n": n,
                "family": args.family, "value": value, "elapsed_s": round(elapsed, 1),
                "threads": args.threads, "ht_bits": args.ht_bits})
    print(f"{args.family}({n}) = {value}   [{len(children)} children, "
          f"{elapsed:.0f}s, threads={args.threads}]")
    return 0


# ------------------------------------------------------------------- sweep

def sweep(args):
    """One shared-memo process over many B(n,j) columns.

    The memo is shared across all columns in the same process, which removes
    almost all repeated work between columns (measured at n = 24: the first
    column alone needs 3.44M states, each further column adds 0.17-0.57M).
    Completed lines survive both the budget exit (code 3) and timeouts thanks
    to the per-line flush, so a partial run still yields the leading columns.
    """
    n = args.n
    if args.columns:
        columns = [int(c) for c in args.columns.split(",")]
    else:
        # orbit representatives: 1..n/2 for even n, 1..(n+1)/2 for odd n
        # (the centre column of an odd strip is its own mirror image)
        columns = list(range(1, (n + 1) // 2 + 1))
    env_bits = args.ht_bits
    binary = args.binary
    text = "".join(f"{n} {family_mask('B', n, j):x}\n" for j in columns)
    env = dict(os.environ)
    env["STRIP3_HT_BITS"] = str(env_bits)
    if args.load is not None:
        env["STRIP3_HT_LOAD"] = str(args.load)
    if args.prealloc:
        env["STRIP3_HT_PREALLOC"] = "1"
    if args.progress:
        env["STRIP3_PROGRESS"] = str(args.progress)
    if args.threads:
        env["STRIP3_THREADS"] = str(args.threads)
    started = time.monotonic()
    status = "ok"
    stderr = ""
    try:
        proc = subprocess.run([str(binary), "masks"], input=text, text=True,
                              capture_output=True, timeout=args.timeout, env=env)
        code, stderr, lines = proc.returncode, proc.stderr, proc.stdout.splitlines()
    except subprocess.TimeoutExpired as exc:
        code, status = None, "timeout"
        raw_out = exc.stdout or b""
        lines = [l for l in (raw_out.decode() if isinstance(raw_out, bytes)
                             else raw_out).splitlines() if l.strip()]
        stderr = (exc.stderr or b"").decode() if isinstance(exc.stderr, bytes) \
            else (exc.stderr or "")
    elapsed = time.monotonic() - started
    if code not in (0, None) and status == "ok":
        status = f"exit_{code}"
    results = []
    union_memo = None
    memo_at_exit = None
    for line in stderr.splitlines():
        if line.startswith("done ") and "memo=" in line:
            union_memo = int(line.split("memo=")[1].split()[0])
        if "memo budget reached at" in line:
            memo_at_exit = int(line.split("memo budget reached at")[1]
                               .split()[0])
    for idx, j in enumerate(columns):
        value = int(lines[idx]) if idx < len(lines) else None
        results.append({"column": j, "mask": "%x" % family_mask("B", n, j),
                        "value": value})
    out = {"version": VERSION, "kind": "strip3_sweep", "n": n,
           "columns": columns, "ht_bits": env_bits, "load": args.load,
           "prealloc": args.prealloc, "threads": args.threads,
           "timeout_s": args.timeout,
           "status": status, "returncode": code,
           "elapsed_s": round(elapsed, 1), "union_memo": union_memo,
           "memo_at_exit": memo_at_exit,
           "binary_sha256": sha256_file(binary) if binary.exists() else None,
           "results": results,
           "note": "one process, shared memo; UNKNOWN means the process ended "
                   "before that line was emitted; fully covered columns are "
                   "exact values",
           "stderr_tail": stderr.strip().splitlines()[-5:] if stderr else []}
    Path(args.output).write_text(json.dumps(out, indent=2) + "\n")
    append_run({"kind": "sweep", "output": str(args.output), "n": n,
                "columns": columns, "status": status,
                "values": {r["column"]: r["value"] for r in results},
                "union_memo": union_memo, "elapsed_s": round(elapsed, 1),
                "ht_bits": env_bits, "load": args.load,
                "prealloc": bool(args.prealloc)})
    solved = sum(r["value"] is not None for r in results)
    print(f"sweep n={n}: {solved}/{len(columns)} columns exact "
          f"[{status}, {elapsed:.0f}s, union memo={union_memo}]")
    for r in results:
        print(f"  B({n},{r['column']}) = "
              f"{'UNKNOWN' if r['value'] is None else r['value']}")
    return 0


def append_run(entry):
    entries = []
    if RUN_LOG.exists():
        try:
            entries = json.loads(RUN_LOG.read_text())
        except json.JSONDecodeError:
            entries = []
    if not isinstance(entries, list):
        entries = []
    entry = {"timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"), **entry}
    entries.append(entry)
    RUN_LOG.write_text(json.dumps(entries, indent=2) + "\n")


# --------------------------------------------------------------------- CLI

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--binary", type=Path, default=DEFAULT_BINARY)
    parser.add_argument("--reductions", type=int, default=None)
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("survey")
    p.add_argument("--max-n", type=int, default=12)
    p.add_argument("--ht-bits", type=int, default=24)
    p.add_argument("--timeout", type=float, default=300)
    p.add_argument("--output", type=Path,
                   default=ROOT / "data/strip3_ports_n12.jsonl")
    p.set_defaults(func=survey)

    p = sub.add_parser("closure")
    p.add_argument("--max-n", type=int, default=8)
    p.add_argument("--output", type=Path,
                   default=ROOT / "data/strip3_ports_closure.json")
    p.set_defaults(func=closure)

    p = sub.add_parser("target")
    p.add_argument("--n", type=int, required=True)
    p.add_argument("--columns", required=True)
    p.add_argument("--ht-bits", type=int, default=27)
    p.add_argument("--timeout", type=float, default=1200)
    p.add_argument("--output", type=Path, required=True)
    p.set_defaults(func=target)

    p = sub.add_parser("recompute")
    p.add_argument("--max-n", type=int, default=6)
    p.set_defaults(func=recompute)

    p = sub.add_parser("collect")
    p.add_argument("--n", type=int, required=True)
    p.add_argument("--source-dir", type=Path, required=True,
                   help="directory with summary_*.txt, <grp>_<seq>.value/.log")
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--ht-bits", type=int, default=30)
    p.add_argument("--timeout", type=float, default=3600)
    p.add_argument("--machine", default="")
    p.add_argument("--strip3-sha256", default="")
    p.set_defaults(func=collect)

    p = sub.add_parser("value")
    p.add_argument("--n", type=int, required=True)
    p.add_argument("--family", required=True,
                   help="family label for the record (A, D, B, E, H, F, ...)")
    p.add_argument("--mask", required=True, help="hex column-major mask")
    p.add_argument("--ht-bits", type=int, default=28)
    p.add_argument("--load", type=int, default=None)
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--prealloc", action="store_true")
    p.add_argument("--timeout", type=float, default=7200)
    p.add_argument("--output", type=Path, required=True)
    p.set_defaults(func=position_value)

    p = sub.add_parser("sweep")
    p.add_argument("--n", type=int, required=True)
    p.add_argument("--columns", help="comma list; default 1..n//2")
    p.add_argument("--ht-bits", type=int, default=28)
    p.add_argument("--load", type=int, default=None,
                   help="hash load percent (10..95); default keeps the binary default")
    p.add_argument("--prealloc", action="store_true",
                   help="allocate the full table upfront (no doubling)")
    p.add_argument("--progress", type=int, default=0,
                   help="print a progress line to stderr every N calls")
    p.add_argument("--threads", type=int, default=1,
                   help="worker threads in one shared-memo process (needs --prealloc)")
    p.add_argument("--timeout", type=float, default=7200)
    p.add_argument("--output", type=Path, required=True)
    p.set_defaults(func=sweep)

    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
