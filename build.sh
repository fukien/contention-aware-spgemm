#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cmake -S "$ROOT" -B "$ROOT/build" "$@"
cmake --build "$ROOT/build" --parallel "${BUILD_JOBS:-4}"
