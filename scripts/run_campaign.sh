#!/usr/bin/env bash
# Data campaign: compute Grundy values for grid boards and strips.
# Re-runnable: a finished job (marker file .done) is skipped. Partial results
# are kept because the solver flushes one TSV line per finished board.
#
# Usage:  scripts/run_campaign.sh strips|boards|all
set -u
cd "$(dirname "$0")/.."
BIN=./src/cfast/kayles
mkdir -p data/raw logs

run() { # run <name> <cap_bits> <outer_timeout_secs> <args...>
  local name=$1 cap=$2 tmo=$3; shift 3
  if [ -e "data/raw/$name.done" ]; then echo "skip  $name (done)"; return; fi
  if [ -e "data/raw/$name.pid" ] && kill -0 "$(cat "data/raw/$name.pid")" 2>/dev/null; then
    echo "skip  $name (already running, pid $(cat "data/raw/$name.pid"))"; return
  fi
  echo "start $name : kayles $* (cap 2^$cap, timeout ${tmo}s)"
  nohup bash -c "timeout $tmo env KAYLES_HT_MAXBITS=$cap '$BIN' $* \
      > 'data/raw/$name.tsv' 2> 'logs/$name.log'; \
      rc=\$?; if [ \$rc -eq 0 ] || [ \$rc -eq 124 ]; then touch 'data/raw/$name.done'; fi; exit \$rc" \
      > /dev/null 2>&1 &
  echo $! > "data/raw/$name.pid"
}

strips() {
  # 1xn is Dawson's chess: computed by a dedicated O(n^2) recurrence (octal 0.137)
  if [ ! -e "data/raw/s1_1xn.done" ]; then
    echo "start s1_1xn : python3 scripts/dawsons_chess.py"
    python3 scripts/dawsons_chess.py > logs/s1_1xn.log 2>&1 && touch data/raw/s1_1xn.done
  else echo "skip  s1_1xn (done)"; fi
  run s2_2xn   26 900    strip 2 150
  run s3_3xn   25 7200   strip 3 32      # published only to n=17 (A316632)
  run s4_4xn   25 7200   strip 4 26      # theorem: all zero; independent check deep
  run s5_5xn   24 5400   strip 5 18
  run s6_6xn   24 5400   strip 6 15
  run s7_7xn   24 3600   strip 7 13
  run s8_8xn   24 3600   strip 8 12
}

boards() {
  run b9_9xn   25 7200   strip 9 11
  run b10_10xn 25 7200   strip 10 11
  run b11_11xn 25 7200   strip 11 11
  run b12_12xn 25 7200   strip 12 12
}

case "${1:-all}" in
  strips) strips ;;
  boards) boards ;;
  all) strips; boards ;;
  *) echo "usage: $0 strips|boards|all"; exit 2 ;;
esac
sleep 2
ps -eo pid,pcpu,rss,args | grep "[c]fast/kayles" | head -30
