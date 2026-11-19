import csv
import hashlib
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


class DownloadIntegrity(unittest.TestCase):
    def test_recheck_and_replace_same_size_file(self):
        script = Path(__file__).resolve().parents[1] / "scripts" / "00_download_data.py"
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "source"
            source.write_bytes(b"good")
            dest = root / "cache"
            dest.mkdir()
            target = dest / "matrix"
            target.write_bytes(b"oops")
            row = dict(tier="core", directory="example", name="matrix", relative_path="matrix",
                       url=source.as_uri(), size_bytes=4, md5=hashlib.md5(b"good").hexdigest())
            manifest = root / "manifest.csv"
            with manifest.open("w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=row)
                writer.writeheader()
                writer.writerow(row)
            cmd = [sys.executable, str(script), "--dest", str(dest), "--manifest", str(manifest)]
            subprocess.run(cmd, check=True, capture_output=True)
            self.assertEqual(target.read_bytes(), b"good")
            target.unlink()
            source.write_bytes(b"oops")
            self.assertNotEqual(subprocess.run(cmd, capture_output=True).returncode, 0)
            self.assertFalse(target.exists())


if __name__ == "__main__":
    unittest.main()
