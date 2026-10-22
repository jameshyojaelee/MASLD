#!/usr/bin/env python3
"""H26: fetch the gene-burden members of 'Raw summary statistics.zip' (Zenodo 10.5281/zenodo.14804239,
version of concept record 14804238; Chen R et al. 2025) with HTTP range requests, without downloading the
1.87 GB archive. Standard library only. Run with `python -I`. Downloaded files are untrusted data.

Usage: python -I fetch_zenodo_gene_burden.py <new_or_empty_out_dir> [--list-only]
"""
import hashlib
import io
import sys
import urllib.request
import zipfile
from pathlib import Path

URL = "https://zenodo.org/api/records/14804239/files/Raw%20summary%20statistics.zip/content"
WANT_PREFIX = "gb_"   # gene-burden members only


class HttpRangeFile(io.RawIOBase):
    """Read-only seekable file over HTTP range requests, with a small block cache."""

    def __init__(self, url, block=4 << 20):
        self.url, self.block, self.pos, self.cache = url, block, 0, {}
        req = urllib.request.Request(url, headers={"Range": "bytes=0-0"})
        with urllib.request.urlopen(req, timeout=120) as r:
            self.size = int(r.headers["Content-Range"].split("/")[1])

    def seekable(self):
        return True

    def readable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, off, whence=0):
        self.pos = off if whence == 0 else self.pos + off if whence == 1 else self.size + off
        return self.pos

    def _block(self, i):
        if i not in self.cache:
            a, b = i * self.block, min((i + 1) * self.block, self.size) - 1
            req = urllib.request.Request(self.url, headers={"Range": f"bytes={a}-{b}"})
            with urllib.request.urlopen(req, timeout=300) as r:
                self.cache[i] = r.read()
            if len(self.cache) > 8:
                self.cache.pop(next(iter(self.cache)))
        return self.cache[i]

    def read(self, n=-1):
        if n < 0:
            n = self.size - self.pos
        out = bytearray()
        while n > 0 and self.pos < self.size:
            i, o = divmod(self.pos, self.block)
            chunk = self._block(i)[o:o + n]
            out += chunk
            self.pos += len(chunk)
            n -= len(chunk)
        return bytes(out)

    def readinto(self, b):
        d = self.read(len(b))
        b[:len(d)] = d
        return len(d)


def main():
    out = Path(sys.argv[1])
    list_only = "--list-only" in sys.argv
    zf = zipfile.ZipFile(io.BufferedReader(HttpRangeFile(URL), buffer_size=1 << 20))
    rows = []
    for info in zf.infolist():
        name = Path(info.filename).name
        print(f"{info.file_size}\t{info.compress_size}\t{info.filename}")
        if list_only or info.is_dir() or "__MACOSX" in info.filename or not name.startswith(WANT_PREFIX):
            continue
        dest = out / name          # flat, basename only: no path from the archive is trusted
        if dest.exists():
            sys.exit(f"refusing to overwrite {dest}")
        h = hashlib.sha256()
        with zf.open(info) as src, open(dest, "wb") as dst:
            while True:
                buf = src.read(1 << 22)
                if not buf:
                    break
                h.update(buf)
                dst.write(buf)
        rows.append(f"{h.hexdigest()}  {name}  crc32={info.CRC:08x}  bytes={info.file_size}  member={info.filename}")
    if rows:
        (out / "raw_gene_burden_members.sha256").write_text("\n".join(rows) + "\n")


if __name__ == "__main__":
    main()
