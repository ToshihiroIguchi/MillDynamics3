import http.client
import tempfile
import threading
import unittest
import zipfile
from pathlib import Path

from milldynamics.fetch import ensure_static
from milldynamics.server import make_server


class ServerTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory()
        root = Path(cls._tmp.name)
        (root / "index.html").write_text("<html></html>", encoding="utf-8")
        (root / "a.wasm").write_bytes(b"\x00asm")
        (root / "a.js").write_text("1", encoding="utf-8")
        cls.httpd = make_server("127.0.0.1", 0, root)
        cls.port = cls.httpd.server_address[1]
        cls.thread = threading.Thread(target=cls.httpd.serve_forever, daemon=True)
        cls.thread.start()

    @classmethod
    def tearDownClass(cls):
        cls.httpd.shutdown()
        cls.httpd.server_close()
        cls._tmp.cleanup()

    def get(self, path):
        conn = http.client.HTTPConnection("127.0.0.1", self.port, timeout=5)
        conn.request("GET", path)
        resp = conn.getresponse()
        resp.read()
        conn.close()
        return resp

    def test_root_redirects(self):
        r = self.get("/")
        self.assertEqual(r.status, 302)
        self.assertEqual(r.getheader("Location"), "/MillDynamics3/")

    def test_no_trailing_slash_redirects(self):
        r = self.get("/MillDynamics3")
        self.assertEqual(r.status, 301)
        self.assertEqual(r.getheader("Location"), "/MillDynamics3/")

    def test_index(self):
        r = self.get("/MillDynamics3/")
        self.assertEqual(r.status, 200)
        self.assertTrue(r.getheader("Content-Type").startswith("text/html"))

    def test_mime_types(self):
        self.assertEqual(self.get("/MillDynamics3/a.wasm").getheader("Content-Type"), "application/wasm")
        self.assertTrue(self.get("/MillDynamics3/a.js").getheader("Content-Type").startswith("text/javascript"))

    def test_outside_prefix_404(self):
        self.assertEqual(self.get("/index.html").status, 404)


class EnsureStaticTest(unittest.TestCase):
    def test_extracts_zip(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            z = tmp / "b.zip"
            with zipfile.ZipFile(z, "w") as zf:
                zf.writestr("index.html", "hi")
                zf.writestr("assets/x.js", "1")
            dest = tmp / "static"
            ensure_static(dest, z.as_uri())
            self.assertEqual((dest / "index.html").read_text(encoding="utf-8"), "hi")
            self.assertTrue((dest / "assets" / "x.js").exists())

    def test_rejects_zip_slip(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            z = tmp / "evil.zip"
            with zipfile.ZipFile(z, "w") as zf:
                zf.writestr("index.html", "hi")
                zf.writestr("../evil.txt", "x")
            dest = tmp / "static"
            with self.assertRaises(RuntimeError):
                ensure_static(dest, z.as_uri())
            self.assertFalse(dest.exists())
            self.assertFalse((tmp / "evil.txt").exists())


if __name__ == "__main__":
    unittest.main()
