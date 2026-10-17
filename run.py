"""Run the simulator straight from a checkout: `python run.py` (no install, no Rust/Node)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / "python"))

from milldynamics.__main__ import main  # noqa: E402

main()
