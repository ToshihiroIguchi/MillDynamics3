#!/usr/bin/env bash
# Builds crates/mill-wasm and emits the wasm-bindgen "web" target bundle into web/src/wasm/
# (gitignored; regenerate with this script after any change to crates/mill-core or crates/mill-wasm).
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"

wasm-pack build "${repo_root}/crates/mill-wasm" \
  --target web \
  --release \
  --out-dir "${repo_root}/web/src/wasm" \
  --out-name mill_wasm
