"""Minimal static file server for the prebuilt web app (standard library only)."""

import functools
import http.server
import threading
import webbrowser
from pathlib import Path

STATIC_DIR = Path(__file__).parent / "static"
# The bundle is built with Vite `base: "/MillDynamics3/"`, so every asset URL is rooted here.
BASE_PATH = "/MillDynamics3/"

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
            self.send_response(302)
            self.send_header("Location", BASE_PATH)
            self.end_headers()
            return
        super().do_GET()

    def translate_path(self, path: str) -> str:
        if path.startswith(BASE_PATH):
            path = "/" + path[len(BASE_PATH):]
        else:
            path = "/__not_found__"
        return super().translate_path(path)

    def log_message(self, format: str, *args) -> None:  # noqa: A002
        pass


def serve(host: str = "127.0.0.1", port: int = 8765, open_browser: bool = True) -> None:
    """Serve the simulator until interrupted (Ctrl+C)."""
    if not (STATIC_DIR / "index.html").exists():
        raise RuntimeError(
            f"{STATIC_DIR} has no index.html; run `python scripts/bundle_python.py` first."
        )
    handler = functools.partial(_Handler, directory=str(STATIC_DIR))
    with http.server.ThreadingHTTPServer((host, port), handler) as httpd:
        url = f"http://{host}:{httpd.server_address[1]}{BASE_PATH}"
        print(f"MillDynamics3 running at {url}  (Ctrl+C to stop)")
        if open_browser:
            threading.Timer(0.3, webbrowser.open, args=(url,)).start()
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            print("stopped")
