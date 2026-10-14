#!/usr/bin/env bash
set -euo pipefail

if [[ -n ${CI_ENV_SCRIPT:-} ]]; then
  set +u # Vendor environment scripts may reference unset variables.
  # shellcheck disable=SC1090
  source "$CI_ENV_SCRIPT"
  set -u
fi

build_dir="$PWD/$CI_BUILD_DIR"
mkdir -p "$build_dir/logs"

make -j"$(nproc)" COMPILER="$CI_COMPILER" "CC=$CC" "CXX=$CXX" "BIN_DIR=$build_dir/bin"

for binary in "$build_dir"/bin/tsvc_*; do
  log="$build_dir/logs/$(basename "$binary").log"
  printf 'Running %s\n' "$(basename "$binary")"
  if [[ -z ${CI_RUN_TIMEOUT:-} ]]; then
    "$binary" | tee "$log"
  elif ! timeout "$CI_RUN_TIMEOUT" "$binary" | tee "$log"; then
    status=${PIPESTATUS[0]}
    if [[ $status != 124 ]]; then exit "$status"; fi
    printf 'Timed out after %s seconds\n' "$CI_RUN_TIMEOUT"
  fi
done
