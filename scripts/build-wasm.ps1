# Builds crates/mill-wasm and emits the wasm-bindgen "web" target bundle into web/src/wasm/
# (gitignored; regenerate with this script after any change to crates/mill-core or crates/mill-wasm).
$ErrorActionPreference = "Stop"

$repoRoot = Split-Path -Parent $PSScriptRoot

wasm-pack build (Join-Path $repoRoot "crates/mill-wasm") `
  --target web `
  --release `
  --out-dir (Join-Path $repoRoot "web/src/wasm") `
  --out-name mill_wasm
