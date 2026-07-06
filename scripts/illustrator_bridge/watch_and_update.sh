#!/usr/bin/env bash
# watch_and_update.sh  (OPTIONAL — the "fully automatic" loop)
# RUN ON THE MAC. Watches the local figure mirror; whenever a PDF changes it
# tells Illustrator to reload all linked files in the front-most document.
# Uses plain osascript, so it does NOT require the MCP to be running.
#
# Requires: fswatch  ->  brew install fswatch
#
# Usage:
#   ./watch_and_update.sh
# Pair it with a periodic sync (cron/launchd running sync_figures_from_hpc.sh),
# or just run the sync manually when you regenerate figures on the HPC.
set -euo pipefail

LOCAL_DIR="${LOCAL_DIR:-$HOME/masld_figures/}"
JSX="${JSX:-$HOME/masld_illustrator/update_links.jsx}"

command -v fswatch >/dev/null 2>&1 || { echo "Install fswatch first:  brew install fswatch"; exit 1; }
[ -f "$JSX" ] || { echo "update_links.jsx not found at: $JSX"; exit 1; }

echo "Watching $LOCAL_DIR for PDF changes; reloading Illustrator links on change."
echo "Ctrl-C to stop."

# -o = batch changes into a single event count; --latency debounces rapid rsync writes.
# -e '.*' excludes everything, then -i '\.pdf$' re-includes only PDFs.
fswatch -o --latency 3 -e '.*' -i '\.pdf$' "$LOCAL_DIR" | while read -r _; do
  echo "[$(date '+%H:%M:%S')] change detected -> reloading Illustrator links"
  osascript <<OSA || echo "  (Illustrator not reachable / no document open — skipped)"
tell application "Adobe Illustrator"
  do javascript POSIX file "$JSX"
end tell
OSA
done
