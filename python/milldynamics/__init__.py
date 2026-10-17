"""MillDynamics3: 2D tumbling ball mill simulator, runnable from Python alone.

The simulator core (Rust) and UI (TypeScript) are shipped prebuilt as static files inside this
package, so no Rust or Node.js toolchain is needed. ``milldynamics.serve()`` hosts them on a
local HTTP server and opens the browser.
"""

from .server import serve

__all__ = ["serve"]
__version__ = "0.1.0"
