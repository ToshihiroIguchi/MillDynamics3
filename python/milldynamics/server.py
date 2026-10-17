"""Minimal static file server for the prebuilt web app (standard library only)."""

import functools
import http.server
import threading
import webbrowser
from pathlib import Path
from typing import Optional

from .fetch import STATIC_DIR, ensure_static

# The bundle is built with Vite `base: "/MillDynamics3/"`, so every asset URL is rooted here.
BASE_PATH = "/MillDynamics3/"
DEFAULT_PORT = 8765

_TYPES = {
    ".wasm": "application/wasm",
    ".js": "text/javascript",
    ".mjs": "text/javascript",
    ".css": "text/css",
    ".svg": "image/svg+xml",
}


class _Handler(http.server.SimpleHTTPRequestHandler):
    extensions_map = {**http.server.SimpleHTTPRequestHandler.extensions_map, **_TYPES}

    def do_GET(self) -> None:
        if self.path in ("", "/"):
            self._redirect(302, BASE_PATH)
            return
        if self.path == BASE_PATH.rstrip("/"):
            self._redirect(301, BASE_PATH)
            return
        super().do_GET()

    def _redirect(self, code: int, location: str) -> None:
        self.send_response(code)
        self.send_header("Location", location)
        self.send_header("Content-Length", "0")
        self.end_headers()

    def translate_path(self, path: str) -> str:
        if path.startswith(BASE_PATH):
            path = "/" + path[len(BASE_PATH):]
        else:
            path = "/__not_found__"
        return super().translate_path(path)

    def log_message(self, format: str, *args) -> None:  # noqa: A002
        pass


def make_server(
    host: str = "127.0.0.1", port: int = 8765, static_dir: Path = STATIC_DIR
) -> http.server.ThreadingHTTPServer:
    """Create (but do not start) the HTTP server serving ``static_dir``."""
    handler = functools.partial(_Handler, directory=str(static_dir))
    return http.server.ThreadingHTTPServer((host, port), handler)


def serve(
    host: str = "127.0.0.1", port: Optional[int] = None, open_browser: bool = True
) -> None:
    """Serve the simulator until interrupted (Ctrl+C).

    ``port=None`` means the default (8765), falling back to a free port if it is busy; an
    explicit port that is busy is an error.
    """
    ensure_static(STATIC_DIR)
    try:
        httpd = make_server(host, DEFAULT_PORT if port is None else port)
    except OSError as exc:
        if port is not None:
            raise RuntimeError(
                f"Port {port} is not available ({exc}); choose another with --port, or 0 for any free port."
            ) from exc
        httpd = make_server(host, 0)
    with httpd:
        url = f"http://{host}:{httpd.server_address[1]}{BASE_PATH}"
        print(f"MillDynamics3 running at {url}  (Ctrl+C to stop)", flush=True)
        if open_browser:
            threading.Timer(0.3, webbrowser.open, args=(url,)).start()
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            print("stopped")
