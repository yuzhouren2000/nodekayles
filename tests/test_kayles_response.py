"""Small raw-graph contexts, transfer SG and certificate damage checks."""
from copy import deepcopy
import json
import os
from pathlib import Path
import random
import resource
import signal
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from kayles_response_equivalence import ResponseAlgebra, ResponseTransfer, column_example
from kayles_boundary_transfer import ThreeRowTransfer, column_graph
from test_kayles_boundary import graphs, labels, glue, raw_sg
from test_kayles_boundary_transfer import random_graph, random_labels, named_join
from test_strip3 import subset_table
from verify_kayles_response import verify


def main():
    os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
    resource.setrlimit(resource.RLIMIT_AS, (512*1024**2,)*2)
    resource.setrlimit(resource.RLIMIT_CPU, (3,)*2)
    signal.setitimer(signal.ITIMER_REAL, 2)
    started = time.monotonic()
    reduced = ResponseAlgebra()
    checks = 0
    for n in range(1, 5):
        for k in range(1, min(n, 3)+1):
            for graph in graphs(n):
                la = labels(n, k)
                state = reduced.state(graph, la)
                assert reduced.form_value(state[0]) ^ state[1] == raw_sg(glue(graph, la, (), ()))
                for context in graphs(k+1):
                    lb = labels(k+1, k)
                    other = reduced.state(context, lb)
                    assert reduced.glued_value(state, other) == raw_sg(glue(graph, la, context, lb))
                    checks += 1
    assert reduced.metrics["merged_forms"] > 0
    print(f"{checks} actual small-graph contexts agree after response reductions", flush=True)

    # Exercise all transfer operations, including deleted/forgotten ports, on
    # arbitrary input masks, not just the reachable-history sublanguage.
    transfer = ThreeRowTransfer(ResponseTransfer())
    expected = subset_table(3)
    for mask, g in enumerate(expected):
        word = tuple(mask >> (3*c) & 7 for c in range(3))
        assert transfer.value(word) == g
    rng = random.Random(3900)
    random_contexts = 0
    for _ in range(60):
        word = tuple(rng.randrange(8) for _ in range(rng.randrange(1, 6)))
        graph, cells = column_graph(word)
        state = transfer.state(word)
        k = max(1, word[-1].bit_count())
        if not word[-1]:
            assert transfer.algebra.value(state) == raw_sg(glue(graph, (0,)*len(graph), (), ()))
            continue
        active_labels = [r+1 for r in range(3) if word[-1] >> r & 1]
        la = tuple(r+1 if c == len(word)-1 else 0 for r, c in cells)
        for _ in range(3):
            context = rng.choice(tuple(graphs(k+1)))
            lb = tuple(active_labels)+(0,)
            other = transfer.algebra.state(context, lb)
            assert transfer.algebra.glued_value(state, other) == raw_sg(glue(graph, la, context, lb))
            random_contexts += 1
    print(f"all 512 induced 3x3 transfer values and {random_contexts} random-context joins agree", flush=True)

    # Semantic equivalence must survive composition when the other side has
    # additional labels, and when some former ports become internal vertices.
    composed = transfer.algebra
    for _ in range(50):
        items = [(random_graph(rng, n), random_labels(rng, n)) for n in (5, 4, 3)]
        states = [composed.state(g, la) for g, la in items]
        joined = composed.join(states[0], states[1])
        assert composed.value(composed.join(joined, states[2])) == raw_sg(named_join(*items))
        forgotten = composed.forget(states[0], (1, 3))
        a, la = items[0]; c, lc = items[2]
        la = tuple(0 if i in (1, 3) else i for i in la)
        lc = tuple(0 if i in (1, 3) else i for i in lc)
        assert composed.value(composed.join(forgotten, composed.state(c, lc))) == raw_sg(glue(a, la, c, lc))
        mapping = {i: i+10 for i in range(1, 5)}
        renamed = composed.rename(joined, mapping)
        assert composed.value(renamed) == composed.value(joined)
    print("50 partial-overlap compositions, port forgetting and renamings agree with actual graphs", flush=True)

    models = tuple(column_example(word) for word in ("573137", "463137"))
    proof = ResponseAlgebra(reduce=False).certificate(models, 1)
    proof = json.loads(json.dumps(proof))
    checked = verify(proof)
    assert checked["pair_facts"] == 32
    assert proof["interfaces"]["states"][0] != proof["interfaces"]["states"][1]
    # Independently test the universal XOR claim on actual graphs. The proof,
    # not this bounded test, is what supplies arbitrary-context validity.
    examples = 0
    for context in graphs(4):
        lb = (1, 2, 3, 0)
        values = [raw_sg(glue(g, la, context, lb)) for g, la in models]
        assert values[0] ^ values[1] == 1
        examples += 1
    invalid = []
    item = deepcopy(proof); item["xor"] = 0; invalid.append(item)
    item = deepcopy(proof); item["pairs"][-1]["internal"].pop(); invalid.append(item)
    item = deepcopy(proof); item["pairs"][-1]["ports"].pop(); invalid.append(item)
    item = deepcopy(proof); item["pairs"][-1]["deletions"].pop(); invalid.append(item)
    item = deepcopy(proof)
    changed = False
    for row in item["pairs"]:
        for attack in row["internal"]:
            if "nim" in attack["response"]:
                attack["response"]["nim"] = 127; changed = True; break
        if changed: break
    assert changed; invalid.append(item)
    item = deepcopy(proof); item["pairs"].pop(0); invalid.append(item)
    item = deepcopy(proof); item["pairs"][-1]["nim"].pop(); invalid.append(item)
    item = deepcopy(proof); item["interfaces"]["nodes"][-1][3].pop(); invalid.append(item)
    for item in invalid:
        try:
            verify(item)
        except (AssertionError, KeyError, IndexError):
            continue
        raise AssertionError("damaged response proof accepted")
    try:
        tiny = ResponseAlgebra(reduce=False, pair_limit=1)
        tiny.certificate(models, 1)
    except RuntimeError as error:
        assert "incomplete" in str(error)
    else:
        raise AssertionError("pair budget ignored")
    # Both bare games have SG=1, but attaching a leaf distinguishes them. A
    # failed response search is UNKNOWN, never an inequivalence certificate.
    false_pair = ResponseAlgebra(reduce=False)
    a = false_pair.state((2, 1), (1, 0))
    b = false_pair.state((0,), (1,))
    assert not false_pair.pair(a, b, 0)
    context = ((2, 1), (1, 0))
    assert raw_sg(glue((2, 1), (1, 0), *context)) != raw_sg(glue((0,), (1,), *context))
    print(f"64-node/32-fact universal proof verified; {examples} actual contexts and {len(invalid)} corruptions checked", flush=True)
    print(f"single-core batch wall_s={time.monotonic()-started:.3f}")
    signal.setitimer(signal.ITIMER_REAL, 0)


if __name__ == "__main__":
    main()
