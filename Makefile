# convenience targets for the Node Kayles research workspace
PY ?= python3
CC ?= gcc
CFLAGS ?= -O2 -march=native

.PHONY: build test analyze data wmoves clean help

help:
	@echo "make build    - compile the C solver"
	@echo "make test     - Python reference checks vs literature (~1 min)"
	@echo "make analyze  - sequences, periodicity, conjecture monitors, all checks"
	@echo "make data     - run/continue the full data campaign (background, resumable)"
	@echo "make wmoves   - winning first moves for small boards"

build:
	$(CC) $(CFLAGS) -o src/cfast/kayles src/cfast/kayles.c

test:
	$(PY) tests/test_literature.py
	$(PY) tests/dumb_solver.py 6 7

analyze:
	$(PY) scripts/analyze.py

data:
	bash scripts/run_campaign.sh all

wmoves:
	$(PY) scripts/winning_moves.py 6

clean:
	rm -f src/cfast/kayles
