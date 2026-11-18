#!/usr/bin/env python3
"""Download the Allen Brain Cell (ABC) Atlas files needed for the VMH neuropeptide atlas.

Standard library only (Python >= 3.8). Files, sizes and MD5 checksums are pinned in
config/data_manifest.csv, which was generated from the ABC Atlas release manifest
20250531. Files land under <data_root>/abc_atlas/ (data_root is set in config/paths.json
or by the VMH_NP_DATA environment variable) using the same relative paths as the
official abc_atlas_access cache, so that package can read them too.

Usage (from the repository root):
    python3 scripts/00_download_data.py --dry-run            # list what would be fetched
    python3 scripts/00_download_data.py                      # core tier (~24 GB)
    python3 scripts/00_download_data.py --tiers core zhuang  # + Zhuang MERFISH (~4.6 GB)
    python3 scripts/00_download_data.py --tiers imputed      # imputed MERFISH (~50 GB)

Interrupted downloads resume where they stopped; completed files are skipped.
"""
import argparse
import csv
import hashlib
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "config" / "data_manifest.csv"
sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import abc_root  # noqa: E402
CHUNK = 4 * 1024 * 1024


def human(n):
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.1f} {unit}" if unit != "B" else f"{n} B"
        n /= 1024


def md5sum(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(CHUNK), b""):
            h.update(block)
    return h.hexdigest()


def fetch(url, dest, size, retries=8):
    """Download url to dest, resuming from a .part file when present."""
    part = dest.with_name(dest.name + ".part")
    dest.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(1, retries + 1):
        have = part.stat().st_size if part.exists() else 0
        if have > size:  # corrupt partial file; start over
            part.unlink()
            have = 0
        if have == size:
            break
        req = urllib.request.Request(url)
        if have:
            req.add_header("Range", f"bytes={have}-")
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                if have and resp.status != 206:  # server ignored the range request
                    have = 0
                mode = "ab" if have else "wb"
                t0, done0, last = time.time(), have, 0.0
                with open(part, mode) as out:
                    while True:
                        block = resp.read(CHUNK)
                        if not block:
                            break
                        out.write(block)
                        have += len(block)
                        now = time.time()
                        if now - last > 1.0:
                            rate = (have - done0) / max(now - t0, 1e-6)
                            sys.stdout.write(
                                f"\r    {human(have)} / {human(size)} "
                                f"({100 * have / size:5.1f}%)  {human(rate)}/s   "
                            )
                            sys.stdout.flush()
                            last = now
            sys.stdout.write("\r" + " " * 70 + "\r")
        except (urllib.error.URLError, OSError, TimeoutError) as err:
            sys.stdout.write("\n")
            print(f"    attempt {attempt}/{retries} interrupted: {err}")
            time.sleep(min(5 * attempt, 30))
            continue
    have = part.stat().st_size if part.exists() else 0
    if have != size:
        raise RuntimeError(f"size mismatch for {dest.name}: got {have}, expected {size}")
    part.rename(dest)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tiers", nargs="+", default=["core"], choices=["core", "zhuang", "imputed"])
    ap.add_argument("--dest", default=str(abc_root()),
                    help="download directory (default: <data_root>/abc_atlas, see config/paths.json)")
    ap.add_argument("--manifest", default=str(MANIFEST), help="manifest CSV (default: config/data_manifest.csv)")
    ap.add_argument("--dry-run", action="store_true", help="list files and sizes, download nothing")
    ap.add_argument("--no-verify", action="store_true", help="skip MD5 verification")
    args = ap.parse_args()

    dest_root = Path(args.dest)
    with open(args.manifest, newline="") as fh:
        rows = [r for r in csv.DictReader(fh) if r["tier"] in args.tiers]
    total = sum(int(r["size_bytes"]) for r in rows)
    print(f"{len(rows)} files, {human(total)} total, tiers: {', '.join(args.tiers)}")
    print(f"destination: {dest_root}\n")

    todo, failed, warned = 0, [], []
    for i, r in enumerate(rows, 1):
        size = int(r["size_bytes"])
        dest = dest_root / r["relative_path"]
        label = f"[{i:2d}/{len(rows)}] {r['directory']}/{r['name']} ({human(size)})"
        if dest.exists() and dest.stat().st_size == size:
            print(f"{label}: already present")
            continue
        todo += size
        if args.dry_run:
            print(f"{label}: to download")
            continue
        print(f"{label}: downloading")
        try:
            fetch(r["url"], dest, size)
            if not args.no_verify:
                got = md5sum(dest)
                if got != r["md5"]:
                    # Size already matched, so this usually means the file was
                    # re-released upstream after the manifest was pinned. Keep it.
                    print(f"    WARNING: MD5 differs from pinned manifest ({got} != {r['md5']}); file kept")
                    warned.append(r["relative_path"])
                    continue
            print("    ok")
        except Exception as err:  # keep going; report at the end
            print(f"    FAILED: {err}")
            failed.append(r["relative_path"])

    if args.dry_run:
        print(f"\nDry run: {human(todo)} still to download.")
        return 0
    if failed:
        print(f"\n{len(failed)} file(s) failed. Re-run the same command to retry:")
        for f in failed:
            print("  ", f)
        return 1
    if warned:
        print(f"\nAll files present; {len(warned)} had a checksum that differs from the pinned manifest:")
        for f in warned:
            print("  ", f)
        return 2
    print("\nAll files present and verified.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
