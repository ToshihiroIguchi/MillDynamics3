import argparse

from .server import serve


def main() -> None:
    p = argparse.ArgumentParser(prog="milldynamics", description="Run the MillDynamics3 simulator.")
    p.add_argument("--host", default="127.0.0.1")
    p.add_argument("--port", type=int, default=None,
                   help="port (default 8765, or a free port if busy); 0 picks a free port")
    p.add_argument("--no-browser", action="store_true", help="do not open the browser")
    a = p.parse_args()
    try:
        serve(host=a.host, port=a.port, open_browser=not a.no_browser)
    except RuntimeError as exc:
        raise SystemExit(f"error: {exc}")


if __name__ == "__main__":
    main()
