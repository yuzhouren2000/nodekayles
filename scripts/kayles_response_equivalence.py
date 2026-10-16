#!/usr/bin/env python3
"""Sufficient universal interface equivalence by paired two-move responses.

Forms retain the complete interface semantics of kayles_boundary_algebra.
Two forms may also be merged when a finite response proof shows that their
contextual SG values always differ by a fixed nimber. Failure is UNKNOWN,
not a distinguishing context. No claims of completeness or bounded size.
"""
from collections import Counter
from functools import cache
import argparse
import gzip
import json
from pathlib import Path

from kayles_boundary_algebra import BoundaryAlgebra, mex
from kayles_boundary_transfer import BoundaryTransfer, column_graph


class ResponseAlgebra(BoundaryAlgebra):
    def __init__(self, reduce=True, pair_limit=200000):
        self.enabled = reduce
        self.pair_limit = pair_limit
        self.buckets = {}
        self.merges = []
        self.metrics = Counter()
        self.proofs = {}
        self.raw_signatures = {}
        super().__init__()
        self.form_value = cache(self.form_value)
        self.certify = cache(self.certify)

    def _form(self, graph, labels):
        state = super()._form(graph, labels)
        return self.reduce_state(state)

    def reduce_state(self, state):
        if not self.enabled or not state[0]:
            return state
        form, offset = state
        if form in self.raw_signatures:
            representative, delta = self.raw_signatures[form]
            return representative, offset ^ delta
        g = self.form_value(form)
        node = self.nodes[form]
        # All proper descendant states have already been reduced. Match the
        # direct-deletion and port-move observations up to the candidate xor.
        signature = (node[0],
                     tuple((label, child, extra ^ g) for label, (child, extra) in node[1]),
                     tuple((label, child, extra ^ g) for label, (child, extra) in node[2]))
        bucket = self.buckets.setdefault(signature, [])
        for other in bucket:
            delta = g ^ self.form_value(other)
            self.metrics["candidate_pairs"] += 1
            if self.certify(form, other, delta):
                self.raw_signatures[form] = other, delta
                self.merges.append((form, other, delta))
                self.metrics["merged_forms"] += 1
                return other, offset ^ delta
        bucket.append(form)
        self.raw_signatures[form] = form, 0
        return state

    def form_value(self, form):
        self.metrics["value_evaluations"] += 1
        if not form:
            return 0
        node = self.nodes[form]
        self.metrics["value_options"] += len(node[2])+len(node[3])
        return mex(self.form_value(child) ^ extra for _, (child, extra) in node[2]+node[3])

    def pair(self, left, right, target):
        a, x = left
        b, y = right
        return self.certify(a, b, target ^ x ^ y)

    def certify(self, a, b, k):
        """Prove SG(a glued C) xor SG(b glued C) == k for every C.

        Recursion decreases the sum of form ids, then the independent Nim
        heap. Input forms are immutable topologically ordered DAG nodes.
        A success records every deletion, port move and two-move response.
        """
        self.metrics["pair_queries"] += 1
        if self.metrics["pair_queries"] > self.pair_limit:
            raise RuntimeError("response-equivalence pair budget reached; incomplete")
        if a == b:
            return k == 0
        na, nb = self.nodes[a], self.nodes[b]
        if na[0] != nb[0] or self.form_value(a) ^ self.form_value(b) != k:
            return False
        key = a, b, k
        proof = dict(a=a, b=b, k=k, deletions=[], ports=[], internal=[], nim=[])
        da, db, ma, mb = map(dict, (na[1], nb[1], na[2], nb[2]))
        # Outside moves can directly remove any subset of ports. Single-port
        # checks recursively include every further deletion, hence all subsets.
        for label in na[0]:
            if not self.pair(da[label], db[label], k):
                return False
            proof["deletions"].append(label)
        for label in na[0]:
            if self.nodes[ma[label][0]][0] != self.nodes[mb[label][0]][0]:
                return False
            if not self.pair(ma[label], mb[label], k):
                return False
            proof["ports"].append(label)

        def reply(first_side, first_index):
            source, other = (na, nb) if first_side == 0 else (nb, na)
            removed, (child, extra) = source[3][first_index]
            other_form = b if first_side == 0 else a
            post = k ^ extra
            # Match an attack on ports by an internal move in the other core.
            # Identical attacks restore the two copies of the outside context.
            for index, (attack, after) in enumerate(other[3]):
                if attack == removed and self.pair((child, 0), after, post):
                    return dict(side=1-first_side, move=index)
            if removed:
                return None
            # With no port touched, a remaining Nim move or another safe move
            # in the same core is a legal response as well.
            target = self.form_value(child) ^ self.form_value(other_form)
            if target < post and self.certify(child, other_form, target):
                return dict(nim=target)
            for index, (attack, after) in enumerate(self.nodes[child][3]):
                if not attack and self.pair(after, (other_form, 0), post):
                    return dict(side=first_side, move=index)
            return None

        for side, node in enumerate((na, nb)):
            for index in range(len(node[3])):
                response = reply(side, index)
                if response is None:
                    return False
                proof["internal"].append(dict(first_side=side, first_move=index, response=response))
        # An opponent can reduce the independent Nim heap. A safe core move
        # must restore a certified zero sum; no port can disappear on one side.
        for j in range(k):
            response = None
            for side, (node, other_form) in enumerate(((na, b), (nb, a))):
                for index, (attack, after) in enumerate(node[3]):
                    if not attack and self.pair(after, (other_form, 0), j):
                        response = dict(side=side, move=index)
                        break
                if response is not None:
                    break
            if response is None:
                return False
            proof["nim"].append(dict(first=j, response=response))
        self.proofs[key] = proof
        self.metrics["certified_pairs"] += 1
        return True

    def report(self):
        return dict(forms=len(self.nodes)-1, representatives=len({p for p, _ in self.raw_signatures.values()}),
                    pairs_cached=self.certify.cache_info().currsize, **self.metrics)

    def certificate(self, graphs, xor):
        """Export unreduced graph interfaces and only the required pair facts."""
        if self.enabled or len(graphs) != 2:
            raise ValueError("certificate extraction needs two graphs and reduction disabled")
        base = self.proof(graphs)
        left, right = base["states"]
        if not self.pair(left, right, xor):
            raise ValueError("no response equivalence proved")
        used = set()

        def visit(left, right, k):
            a, x = left; b, y = right; k ^= x ^ y
            if a == b:
                assert k == 0
                return
            key = a, b, k
            if key in used:
                return
            used.add(key)
            record = self.proofs[key]
            na, nb = self.nodes[a], self.nodes[b]
            for section in (1, 2):
                for (_, sa), (_, sb) in zip(na[section], nb[section]):
                    visit(sa, sb, k)
            for row in record["internal"]:
                side, index, response = row["first_side"], row["first_move"], row["response"]
                source, other = (na, nb) if side == 0 else (nb, na)
                other_form = b if side == 0 else a
                _, (child, extra) = source[3][index]
                if "nim" in response:
                    visit((child, 0), (other_form, 0), response["nim"])
                elif response["side"] != side:
                    visit((child, 0), other[3][response["move"]][1], k ^ extra)
                else:
                    visit(self.nodes[child][3][response["move"]][1], (other_form, 0), k ^ extra)
            for row in record["nim"]:
                side, move = row["response"]["side"], row["response"]["move"]
                node, other_form = (na, b) if side == 0 else (nb, a)
                visit(node[3][move][1], (other_form, 0), row["first"])

        visit(left, right, xor)
        return dict(format="node-kayles-response-equivalence-v1", interfaces=base, roots=[0, 1], xor=xor,
                    pairs=[self.proofs[key] for key in sorted(used, key=lambda key: (key[0]+key[1], key[2], key))])


class ResponseTransfer(ResponseAlgebra, BoundaryTransfer):
    """Compose the reduced states; do not enumerate the original graph masks."""
    def _finish_node(self, node):
        return self.reduce_state((self._intern_node(node), 0))


def column_example(word):
    graph, cells = column_graph(tuple(map(int, word)))
    return graph, tuple(r+1 if c == len(word)-1 else 0 for r, c in cells)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--certificate", type=Path, required=True)
    args = parser.parse_args()
    # A genuine full-end-column relation. All ports remain exposed to an
    # arbitrary graph, including arbitrarily long complete thick strips.
    graphs = tuple(column_example(word) for word in ("573137", "463137"))
    report = ResponseAlgebra(reduce=False).certificate(graphs, 1)
    report["column_words"] = ["573137", "463137"]
    raw = (json.dumps(report, separators=(",", ":"))+"\n").encode()
    args.certificate.write_bytes(gzip.compress(raw, mtime=0) if args.certificate.suffix == ".gz" else raw)
