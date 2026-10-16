# convenience targets for the Node Kayles research workspace
PY ?= python3
CC ?= gcc
CFLAGS ?= -O2 -march=native

.PHONY: build build-strip3 build-strip3-language export-strip3-comb test-strip3-comb benchmark-strip3-comb test-protected-strip-search test test-strip3 test-strip3-ports test-strip3-target test-strip3-exact test-strip3-two-move test-strip3-interface test-strip3-defect test-comb-stabilization test-clique-strips test-strip3-language test-strip3-language-33 test-kayles-boundary test-kayles-transfer test-guarded-symmetry test-neutral-core test-comb-bridge test-two-comb-core test-parametric-bridge test-connected-comb-core test-three-comb-chain test-interface-classes analyze analyze-strip3 analyze-strip3-ports data wmoves clean help

help:
	@echo "make test-strip3-proof-cost - verify equality certificates and offline minimum-tree recurrence"
	@echo "make test-strip3-absence - check absent-option mex bounds against independent subset SG"
	@echo "make test-strip3-connection-rank - independently check finite ranks, nonsingular minors and reachable joins"
	@echo "make test-strip3-search-experiments - check coupled, symbolic and proof-number prototypes on independent small graphs"
	@echo "make test-strip3-target-search - independently verify the target-nimber prototype (seconds)"
	@echo "make screen-strip3 - bounded local screening: 2 s/process, 20 s total, 512 MiB"
	@echo "make build    - compile the C solver"
	@echo "make build-strip3 - compile the specialised three-row solver"
	@echo "make export-strip3-comb - independently verify and export bounded comb tables"
	@echo "make test-strip3-comb - check comb recognition and SG against recursive search"
	@echo "make benchmark-strip3-comb - locally measure full-board search with/without comb rules"
	@echo "make test-protected-strip-search - verify the bounded-thickness search obstruction and cache witnesses"
	@echo "make build-strip3-language - compile the perfect-index one-byte solver"
	@echo "make test-strip3  - independent induced-position and survey checks"
	@echo "make test-strip3-ports - B/E/H port-family coordinate and data checks"
	@echo "make test-strip3-target - verify the 3x30 edge-orbit conclusion"
	@echo "make test-strip3-exact - re-derive G(3xn) for n=29..32 from sweep data"
	@echo "make test-strip3-two-move - check separator identities and recompute proof premises"
	@echo "make test-strip3-interface - independently check the boundary-summary obstruction"
	@echo "make test-strip3-defect - check move-closed defect algebra and infinite period certificates"
	@echo "make test-comb-stabilization - exhaust the universal block lemma and check tail algebra"
	@echo "make test-clique-strips - check three-row defect algebra and infinite period certificates"
	@echo "make test-strip3-language - check the complete reachable-component language and byte table"
	@echo "make test-strip3-language-33 - verify f(33)=1: a-sequence middle options and mex over the survey certificate"
	@echo "make test-interface-classes - check reachable-family transfer classes and boundary-fixing isomorphisms"
	@echo "make test-kayles-boundary - verify arbitrary-context square/tree replacement and boundary products"
	@echo "make test-kayles-transfer - check complete column generation and infinite-interface premises"
	@echo "make test-guarded-symmetry - check exact guarded-symmetry SG and dense three-row families"
	@echo "make test-neutral-core - verify paired-exterior absorption and arbitrary-prefix tail rewrites"
	@echo "make test-comb-bridge - verify the two-core bridge period and reachable families"
	@echo "make test-two-comb-core - verify the two-arm finite formula certificates"
	@echo "make test-parametric-bridge - verify bridge periods uniform in both outer-arm lengths"
	@echo "make test-connected-comb-core - verify component-normalized comb certificates"
	@echo "make test-three-comb-chain - verify all three lengths and the explicit nim formula"
	@echo "make analyze-strip3 - structural analysis and constructive first moves"
	@echo "make analyze-strip3-ports - regenerate the port-family survey/closure data"
	@echo "make test     - Python reference checks vs literature (~1 min)"
	@echo "make analyze  - sequences, periodicity, conjecture monitors, all checks"
	@echo "make data     - run/continue the full data campaign (background, resumable)"
	@echo "make wmoves   - winning first moves for small boards"

build:
	$(CC) $(CFLAGS) -o src/cfast/kayles src/cfast/kayles.c

build-strip3: src/cfast/strip3_comb_tables.h
	$(CC) $(CFLAGS) -Wall -Wextra -o src/cfast/strip3 src/cfast/strip3.c

export-strip3-comb:
	$(PY) scripts/export_strip3_comb.py --verify

test-strip3-comb:
	$(PY) scripts/export_strip3_comb.py --check
	$(PY) tests/test_strip3_comb.py

benchmark-strip3-comb:
	$(PY) scripts/benchmark_strip3_comb.py

test-protected-strip-search:
	$(PY) tests/verify_protected_strip_search.py

build-strip3-language:
	$(CC) $(CFLAGS) -Wall -Wextra -pthread -o src/cfast/strip3_language src/cfast/strip3_language.c

build-strip3-target:
	$(CC) $(CFLAGS) -Wall -Wextra -pthread -o src/cfast/strip3_target src/cfast/strip3_target.c

test-strip3-target-search:
	$(PY) tests/test_strip3_target.py

test-strip3-absence:
	$(PY) tests/test_strip3_absence.py
	$(PY) tests/test_strip3_sequence_screen.py

test-strip3-proof-cost:
	$(PY) tests/test_strip3_proof_cost.py
	$(PY) tests/test_strip3_target_proof.py
	$(PY) tests/verify_strip3_equality_proof.py data/proofs/strip3_equality/strip3_equality_target_tree_min_26.json.gz

test-strip3-search-experiments:
	$(PY) tests/test_strip3_coupled.py
	$(PY) tests/test_strip3_search_experiments.py
	$(PY) tests/verify_strip3_work_obligations.py

screen-strip3:
	$(PY) scripts/screen_strip3.py --output data/strip3_target_screen.json

test-strip3:
	$(PY) tests/test_strip3.py --survey data/strip3_survey.jsonl
	$(PY) tests/verify_strip3_analysis.py

test-strip3-ports:
	$(PY) tests/test_strip3_ports.py --max-n 12 --max-sg-n 8
	$(PY) tests/verify_strip3_ports.py data/strip3_ports_n12.jsonl
	$(PY) tests/verify_fixed_pairing.py data/proofs/strip3_pairing_[A-Z]_*.json

test-strip3-target:
	$(PY) tests/verify_strip3_target.py

test-strip3-exact:
	$(PY) tests/verify_strip3_exact.py --n 29 --sweep data/strip3_sweep_29.json
	$(PY) tests/verify_strip3_exact.py --n 30 --sweep data/strip3_sweep_30.json
	$(PY) tests/verify_strip3_exact.py --n 31 --sweep data/strip3_sweep_31.json
	$(PY) tests/verify_strip3_exact.py --n 32 --sweep data/strip3_sweep_32.json
	$(PY) tests/verify_strip3_exact.py --n 33 --sweep data/proofs/strip3_language_sweep_33.json
	$(PY) tests/verify_strip3_exact.py --n 34 --sweep data/proofs/strip3_language_sweep_34.json

test-strip3-two-move:
	$(PY) tests/test_strip3_two_move.py
	$(PY) tests/verify_strip3_two_move.py data/proofs/strip3_two_move_4.json
	$(PY) tests/verify_strip3_two_move.py data/proofs/strip3_two_move_22.json
	$(PY) tests/verify_strip3_two_move.py data/proofs/strip3_two_move_26.json
	$(PY) tests/verify_strip3_two_move.py data/proofs/strip3_two_move_30.json

test-strip3-interface:
	$(PY) tests/verify_strip3_interface.py data/proofs/strip3_interface_obstruction.json

test-strip3-defect:
	$(PY) tests/test_strip3_defect_algebra.py
	$(PY) tests/verify_strip3_defect_period.py data/proofs/strip3_defect_period.json

test-comb-stabilization: test-strip3-defect
	$(PY) tests/verify_comb_stabilization.py

test-clique-strips:
	$(PY) tests/test_clique_column_strips.py
	$(PY) tests/verify_clique_strip_period.py data/proofs/clique_strip_period.json

test-strip3-language:
	$(PY) tests/test_strip3_history_language.py

test-strip3-language-33:
	$(PY) tests/verify_strip3_language_33.py

test-interface-classes:
	$(PY) tests/verify_interface_classes.py

test-strip3-connection-rank:
	$(PY) tests/verify_strip3_connection_rank.py data/strip3_connection_rank.json.gz data/strip3_guarded_connection_rank.json

test-kayles-boundary:
	$(PY) tests/verify_kayles_boundary.py data/proofs/kayles_boundary_replacement.json
	$(PY) tests/test_kayles_boundary.py

test-kayles-transfer:
	$(PY) tests/test_kayles_boundary_transfer.py

test-guarded-symmetry:
	$(PY) tests/test_guarded_symmetry.py

test-neutral-core:
	$(PY) tests/test_neutral_core.py
	$(PY) tests/verify_neutral_core.py data/proofs/neutral_core_sg5.json
	$(PY) tests/verify_neutral_core.py data/proofs/paired_neutral_core_9.json
	$(PY) tests/verify_neutral_core.py data/proofs/bipartite_neutral_obstruction_12.json

test-comb-bridge:
	$(PY) tests/test_comb_bridge.py
	$(PY) tests/verify_comb_bridge.py data/proofs/comb_bridge_2_4.json

test-reachable-boundary-search:
	$(PY) tests/test_reachable_boundary_search.py

test-two-comb-core:
	$(PY) tests/test_two_comb_core.py
	$(PY) tests/verify_two_comb_core.py data/proofs/two_comb_square_opposite.json.gz
	$(PY) tests/verify_two_comb_core.py data/proofs/two_comb_strip_middle2.json.gz
	$(PY) tests/verify_two_comb_core.py data/proofs/two_comb_strip_middle3.json.gz
	$(PY) tests/verify_two_comb_core.py data/proofs/two_comb_strip_middle4.json.gz
	$(PY) tests/verify_two_comb_core.py data/proofs/two_comb_strip_middle5.json.gz

test-parametric-bridge:
	$(PY) tests/test_parametric_comb_bridge.py
	$(PY) tests/verify_parametric_comb_bridge.py

test-connected-comb-core:
	$(PY) tests/test_connected_comb_core.py

test-three-comb-chain: test-connected-comb-core
	$(PY) tests/test_three_comb_chain.py
	$(PY) tests/verify_three_comb_chain.py

analyze-strip3:
	$(PY) scripts/analyze_strip3.py

analyze-strip3-ports:
	$(PY) scripts/strip3_ports.py survey --max-n 12 --ht-bits 24 --timeout 300 \
	  --output data/strip3_ports_n12.jsonl
	$(PY) scripts/strip3_ports.py closure --max-n 8 \
	  --output data/strip3_ports_closure.json

test:
	$(PY) tests/test_literature.py
	$(PY) tests/dumb_solver.py 6 7
	$(PY) tests/verify_certificates.py
	$(PY) tests/test_pairing_classification.py
	$(PY) tests/test_fixed_pairing.py

analyze:
	$(PY) scripts/analyze.py

data:
	bash scripts/run_campaign.sh all

wmoves:
	$(PY) scripts/winning_moves.py 6

clean:
	rm -f src/cfast/kayles src/cfast/strip3 src/cfast/strip3_language src/cfast/strip3_target

.PHONY: build-strip3-target test-strip3-target-search screen-strip3
.PHONY: test-strip3-search-experiments
.PHONY: test-strip3-connection-rank
.PHONY: test-strip3-absence
.PHONY: test-strip3-proof-cost
