"""Download the prebuilt web bundle on first use (standard library only)."""

import shutil
import tempfile
import urllib.request
import zipfile
from pathlib import Path

STATIC_DIR = Path(__file__).parent / "static"
DEFAULT_URL = (
    "https://github.com/ToshihiroIguchi/MillDynamics3/releases/download/"
    "python-latest/milldynamics-static.zip"
)


def _safe_extract(zf: zipfile.ZipFile, dest: Path) -> None:
    root = dest.resolve()
    for name in zf.namelist():
        target = (root / name).resolve()
        if target != root and root not in target.parents:
            raise RuntimeError(f"unsafe path in archive: {name!r}")
    zf.extractall(str(root))


def ensure_static(static_dir: Path = STATIC_DIR, url: str = DEFAULT_URL) -> Path:
    """Make sure ``static_dir/index.html`` exists, downloading the bundle if needed."""
    static_dir = Path(static_dir)
    if (static_dir / "index.html").exists():
        return static_dir
    print(f"Web bundle not found; downloading {url} ...")
    try:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            zip_path = tmp_path / "bundle.zip"
            with urllib.request.urlopen(url, timeout=60) as resp, open(zip_path, "wb") as out:
                shutil.copyfileobj(resp, out)
            print("Downloaded; extracting ...")
            extract_dir = tmp_path / "extract"
            extract_dir.mkdir()
            with zipfile.ZipFile(zip_path) as zf:
                _safe_extract(zf, extract_dir)
            if not (extract_dir / "index.html").exists():
                raise RuntimeError("archive has no index.html at its root")
            if static_dir.exists():
                shutil.rmtree(static_dir)
            static_dir.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(extract_dir), str(static_dir))
    except Exception as exc:
        raise RuntimeError(
            f"Could not obtain the web bundle ({exc}). Alternatives: "
            "pip install the release wheel "
            "(https://github.com/ToshihiroIguchi/MillDynamics3/releases/download/python-latest/"
            "milldynamics-0.1.0-py3-none-any.whl), or build it locally with "
            "`python scripts/bundle_python.py`."
        ) from exc
    print("Web bundle ready.")
    return static_dir
