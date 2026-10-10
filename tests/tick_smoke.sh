#!/usr/bin/env bash
# End-to-end test of the real-time path: tick/tp.q -> tick/rdb.q, fed by feed.q.
# Usage (from the repository root, after `q hawkes_quotes.q`):  bash tests/tick_smoke.sh
# Uses ports 15010/15011 and a temporary journal directory; cleans up on exit.
set -euo pipefail
cd "$(dirname "$0")/.."
[ -d hdbq ] || { echo "hdbq/ missing: run 'q hawkes_quotes.q' first"; exit 1; }
DATE=${DATE:-2024.01.02}
TP=${TP_PORT:-15010}; RDB=${RDB_PORT:-15011}
LOG=$(mktemp -d)
pids=()
cleanup() { for p in "${pids[@]}"; do kill "$p" 2>/dev/null || true; done; rm -rf "$LOG"; }
trap cleanup EXIT

q tick/tp.q "$LOG" -p "$TP" > "$LOG/tp.out" 2>&1 < /dev/null & pids+=($!)
sleep 1
q tick/rdb.q ":$TP" -p "$RDB" > "$LOG/rdb.out" 2>&1 < /dev/null & pids+=($!)
sleep 1
# x100000: the 6h30 session is replayed in about a quarter of a second
q feed.q 100000 "$DATE" "$TP" > "$LOG/feed.out" 2>&1 < /dev/null & pids+=($!)

if q tests/tick_check.q "$RDB" "$DATE" 120 < /dev/null; then
  echo "tick smoke test: PASS"
else
  echo "tick smoke test: FAIL"; for f in tp rdb feed; do echo "--- $f"; cat "$LOG/$f.out"; done; exit 1
fi
