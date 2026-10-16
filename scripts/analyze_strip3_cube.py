#!/usr/bin/env python3
"""Audit representation barriers without rerunning expensive game searches.

--generate optionally constructs the small six-column cube proof. Each build,
construction and independent check is single-core, 512 MiB and at most 2 s.
The default reads the saved proof and saved equality DAG; no SG search runs.
"""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time

from history_cube_geometry import (cube_witness, independent_history_count,
                                   neighbourhoods, maximal_history_parities, parity_sg)
from screen_strip3 import atomic_write_json, limits, source_provenance, BUILD_FLAGS

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests"))
from verify_strip3_cube_proof import check_text


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generate", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cpu = min(os.sched_getaffinity(0))
    os.sched_setaffinity(0, {cpu})
    started = time.monotonic()
    proof = ROOT / "data/proofs/strip3_history_cubes_6.txt.gz"
    provenance = source_provenance(ROOT / "src/cfast/strip3_cube.c", True)
    provenance["source_sha256"].update({str(path.relative_to(ROOT)): sha(path) for path in
        (Path(__file__).resolve(), ROOT / "scripts/history_cube_geometry.py",
         ROOT / "tests/verify_strip3_cube_proof.py", ROOT / "tests/test_strip3_cube.py")})
    report = dict(format="strip3-history-cube-audit-v1", provenance=provenance,
                  limits=dict(cpu=cpu, threads=1, process_wall_s=2, address_space_mib=512),
                  work=dict(builds=0, solver_processes=0))

    def bounded(command, env=None):
        remaining = 20 - (time.monotonic()-started)
        assert remaining > 0, "batch budget"
        before = time.monotonic()
        run = subprocess.run(command, text=True, capture_output=True, env=env,
                             timeout=min(2, remaining), preexec_fn=lambda: limits(512, 2))
        assert run.returncode == 0, run.stderr
        return dict(command=list(map(str, command)), wall_s=round(time.monotonic()-before, 6),
                    returncode=run.returncode, stdout=run.stdout.strip(), stderr=run.stderr.strip())

    if args.generate:
        with tempfile.TemporaryDirectory(prefix="strip3-cube-proof-") as folder:
            binary, raw = Path(folder) / "cube", Path(folder) / "proof.txt"
            report["build"] = bounded([provenance["build"]["compiler"], *BUILD_FLAGS,
                                       "-o", str(binary), str(ROOT / "src/cfast/strip3_cube.c")])
            report["work"]["builds"] += 1
            env = {k: v for k, v in os.environ.items() if not k.startswith("STRIP3_")}
            env.update(STRIP3_CUBE_LOOKUP="1", STRIP3_CUBE_BITS="18")
            report["generation"] = bounded([str(binary), "board", "6", "--certificate", str(raw)], env)
            report["work"]["solver_processes"] += 1
            assert report["generation"]["stdout"].split() == ["3", "6", "3"]
            raw_data = raw.read_bytes()
            check_text(raw_data.decode())
            proof.write_bytes(gzip.compress(raw_data, mtime=0))
    report["independent_check"] = bounded([sys.executable, str(ROOT / "tests/verify_strip3_cube_proof.py"), str(proof)])
    raw = gzip.decompress(proof.read_bytes())
    stats = check_text(raw.decode())
    adjacency = neighbourhoods(stats["n"])
    non_singletons = 0
    for line in raw.decode().splitlines()[1:]:
        row = line.split()
        non_singletons += cube_witness(adjacency, int(row[4], 16), int(row[5], 16)) is not None
    assert non_singletons == 0
    report["proof"] = dict(path=str(proof.relative_to(ROOT)), sha256=sha(proof),
                           raw_sha256=hashlib.sha256(raw).hexdigest(), raw_bytes=len(raw),
                           gzip_bytes=proof.stat().st_size, non_singleton_facts=non_singletons, **stats)
    report["full_strip_parity_test"] = [
        dict(n=n, independent_histories=independent_history_count(n),
             maximal_history_parity_support=maximal_history_parities((1 << (3*n))-1),
             SG_from_parity=parity_sg((1 << (3*n))-1)) for n in range(1, 101)]
    assert all(row["SG_from_parity"] is None for row in report["full_strip_parity_test"])
    corpus = ROOT / "data/proofs/strip3_equality/strip3_equality_target_tree_min_26.json.gz"
    records = json.loads(gzip.decompress(corpus.read_bytes()))["nodes"]
    geometry = {int(row["mask"], 16) for row in records}
    hits = {mask: parity_sg(mask) for mask in geometry if parity_sg(mask) is not None}
    for row in records:
        mask = int(row["mask"], 16)
        if mask in hits:
            assert (row["k"] == hits[mask]) == row["equal"]
    report["saved_geometry_parity_test"] = dict(
        corpus=str(corpus.relative_to(ROOT)), corpus_sha256=sha(corpus),
        geometries=len(geometry), hits=len(hits), maximum_hit_vertices=max(map(int.bit_count, hits)),
        hit_masks=[dict(mask=f"{mask:x}", SG=g) for mask, g in sorted(hits.items())])
    report["single_affine_class_counterexample"] = dict(
        n=4, target=0, histories=["0", "a", "82", "88"], values=[0, 0, 0, 2],
        validation="tests/test_strip3_cube.py independent subset mex; 0xa xor 0x82 = 0x88")
    report["conclusion"] = (
        "Equality cubes on independent histories are singletons. The tested local rules also "
        "leave all negative cubes singleton. Tries cannot generalise these facts. Uniform "
        "maximal-history parity is sufficient for SG 0/1 but solves no complete strip through "
        "100 and only tiny geometries in the saved DAG. No new full-board SG value is proved.")
    report["wall_s"] = round(time.monotonic()-started, 6)
    atomic_write_json(args.output, report)
    print(json.dumps(dict(proof=report["proof"], corpus=report["saved_geometry_parity_test"],
                          wall_s=report["wall_s"], work=report["work"])))


if __name__ == "__main__":
    main()
