"""Rebuild the WASM core + web UI and copy the result into the Python package.

Maintainer-only (needs Rust, wasm-pack, Node.js). End users never run this.
"""

import shutil
import subprocess
import sys
from pathlib import Path

root = Path(__file__).resolve().parent.parent
shell = sys.platform == "win32"

subprocess.run(
    ["wasm-pack", "build", str(root / "crates" / "mill-wasm"), "--target", "web", "--release",
     "--out-dir", str(root / "web" / "src" / "wasm"), "--out-name", "mill_wasm"],
    check=True, shell=shell,
)
subprocess.run(["npm", "run", "build"], cwd=root / "web", check=True, shell=shell)

dest = root / "python" / "milldynamics" / "static"
shutil.rmtree(dest, ignore_errors=True)
shutil.copytree(root / "web" / "dist", dest)
print(f"bundled web/dist -> {dest}")
