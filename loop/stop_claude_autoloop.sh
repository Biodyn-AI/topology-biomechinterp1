#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUNTIME_DIR="$ROOT_DIR/runtime"
PID_FILE="$RUNTIME_DIR/claude_autoloop.pid"
STOP_FILE="$RUNTIME_DIR/claude_STOP"
PY_SCRIPT="$ROOT_DIR/loop/run_claude_topology_autoloop.py"
STATUS_FILE="$RUNTIME_DIR/loop_status.json"

mkdir -p "$RUNTIME_DIR"
touch "$STOP_FILE"

find_loop_pid() {
  pgrep -f "$PY_SCRIPT" 2>/dev/null | tail -n 1 || true
}

PID=""
if [[ -f "$PID_FILE" ]]; then
  PID="$(cat "$PID_FILE" || true)"
fi

if [[ -z "$PID" ]] || ! kill -0 "$PID" 2>/dev/null; then
  PID="$(find_loop_pid)"
fi

if [[ -n "$PID" ]] && kill -0 "$PID" 2>/dev/null; then
  kill "${PID}" || true
  sleep 1
  if kill -0 "${PID}" 2>/dev/null; then
    kill -9 "${PID}" || true
  fi

  # Clean up any orphan Claude workers tied to this project root.
  ORPHANS="$(pgrep -f "claude -p .*--add-dir $ROOT_DIR" || true)"
  if [[ -n "$ORPHANS" ]]; then
    echo "$ORPHANS" | xargs kill >/dev/null 2>&1 || true
    sleep 1
    STILL_RUNNING="$(echo "$ORPHANS" | xargs -I{} sh -c 'kill -0 {} 2>/dev/null && echo {}' || true)"
    if [[ -n "$STILL_RUNNING" ]]; then
      echo "$STILL_RUNNING" | xargs kill -9 >/dev/null 2>&1 || true
    fi
  fi

  echo "Stopped autoloop process ${PID}"
  rm -f "$PID_FILE"
else
  echo "No running autoloop process found. STOP flag created."
fi

if [[ -f "$STATUS_FILE" ]]; then
  python3 - "$STATUS_FILE" <<'PY'
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

status_path = Path(sys.argv[1])
payload = json.loads(status_path.read_text(encoding="utf-8"))
payload["running"] = False
payload["phase"] = "stopped"
payload["ts"] = datetime.now(timezone.utc).isoformat()
status_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
PY
fi
