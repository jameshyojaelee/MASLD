"""H49: write the list of Open TG-GATEs human in vitro CEL archives to fetch.

Usage: python3 -I h49_url_list.py <open_tggates_main.csv> <out_urls.tsv>

Reads the TG-GATEs compound table (untrusted download) and converts the
`Human - in vitro` FTP links to the HTTPS mirror on dbarchive.biosciencedbc.jp.
Writes two columns: file name and URL.
"""
import sys
import urllib.parse

import pandas as pd

FTP_PREFIX = "ftp://ftp.biosciencedbc.jp/archive/open-tggates/LATEST/Human/in_vitro/"
HTTPS_PREFIX = "https://dbarchive.biosciencedbc.jp/data/open-tggates/LATEST/Human/in_vitro/"


def main(main_csv: str, out_tsv: str) -> None:
    tab = pd.read_csv(main_csv, dtype=str, encoding="latin-1")
    links = tab["Human - in vitro"].dropna()
    rows = []
    for link in links:
        if not link.startswith(FTP_PREFIX):
            raise SystemExit(f"unexpected link: {link}")
        encoded = link[len(FTP_PREFIX):]
        name = urllib.parse.unquote(encoded)
        if "/" in name or name.startswith(".") or not name.endswith(".Human.in_vitro.Liver.zip"):
            raise SystemExit(f"unexpected file name: {name}")
        rows.append((name, HTTPS_PREFIX + encoded))
    if len(rows) != 158:
        raise SystemExit(f"expected 158 archives, found {len(rows)}")
    with open(out_tsv, "w") as fh:
        for name, url in rows:
            fh.write(f"{name}\t{url}\n")
    print(f"wrote {len(rows)} URLs to {out_tsv}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
