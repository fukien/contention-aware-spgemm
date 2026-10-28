#!/usr/bin/env bash
# Drop Linux filesystem caches before a benchmark (requires root).
# Mode: 1 = page cache, 2 = dentries/inodes, 3 = both (default).
set -euo pipefail
if [ "$#" -gt 1 ]; then
    echo "Usage: sudo bash $0 [1|2|3]" >&2
    exit 1
fi
MODE=${1:-3}
case "$MODE" in
    1|2|3) ;;
    *) echo "Usage: sudo bash $0 [1|2|3]" >&2; exit 1 ;;
esac
if [ "$(uname -s)" != Linux ]; then
    echo "clear.sh requires Linux." >&2
    exit 1
fi
if [ "$EUID" -ne 0 ]; then
    echo "Run with sudo: sudo bash $0 $MODE" >&2
    exit 1
fi
sync
printf '%s\n' "$MODE" > /proc/sys/vm/drop_caches
