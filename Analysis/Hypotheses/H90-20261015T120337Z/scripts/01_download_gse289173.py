"""H90 job 6300 (a): download GSE289173_RAW.tar, record sha256, untar, gate G0."""
import gzip, hashlib, subprocess, sys, tarfile
from pathlib import Path
import pandas as pd
sys.path.insert(0, str(Path(__file__).parent))
import h90_common as C

C.RAW.mkdir(parents=True, exist_ok=True)
if not C.TAR.exists() or C.TAR.stat().st_size != C.TAR_BYTES:
    r = subprocess.run(["wget", "-c", "-q", "-O", str(C.TAR), C.TAR_URL])
    if r.returncode != 0:
        sys.exit("STOP: wget failed on this node; rerun the same wget on the login node (file transfer only)")
size = C.TAR.stat().st_size
tar_sha = C.sha256(C.TAR)
shaf = C.RAW / "GSE289173_RAW.tar.sha256"
if shaf.exists():
    if shaf.read_text().split()[0] != tar_sha:
        sys.exit("STOP: tar sha256 differs from the one recorded at download")
else:
    shaf.write_text(f"{tar_sha}  {C.TAR.name}\n")
with tarfile.open(C.TAR) as t:
    names = t.getnames()
    todo = [m for m in t.getmembers() if not (C.RAW / m.name).exists()]
    if todo:
        t.extractall(C.RAW, members=todo, filter="data")

rows = []
def add(check, value, expected, ok):
    rows.append(dict(check=check, value=value, expected=expected, pass_=bool(ok)))

add("tar_bytes", size, C.TAR_BYTES, size == C.TAR_BYTES)
add("tar_sha256", tar_sha, "recorded at download", True)
mtx = sorted(n for n in names if n.endswith("matrix.mtx.gz"))
add("n_matrix_mtx_gz", len(mtx), 48, len(mtx) == 48)
nps = pd.read_csv(C.GATE / "nuclei_per_sample.tsv", sep="\t", header=None, names=["file", "n"])
nps["gsm"] = nps.file.str.split("_").str[0]
design = pd.read_csv(C.WRITE / "donor_design.tsv", sep="\t")
fl514 = C.read_gz_lines(C.FEATURES_FL514)
per = []
for gsm in design.gsm:
    bc = [n for n in names if n.startswith(gsm + "_") and n.endswith("barcodes.tsv.gz")]
    ft = [n for n in names if n.startswith(gsm + "_") and n.endswith("features.tsv.gz")]
    mx = [n for n in names if n.startswith(gsm + "_") and n.endswith("matrix.mtx.gz")]
    ok_files = len(bc) == len(ft) == len(mx) == 1
    nbc = len(C.read_gz_lines(C.RAW / bc[0])) if ok_files else -1
    fsha = C.sha256(C.RAW / ft[0]) if ok_files else ""
    fsame_bytes = fsha == C.FEATURES_SHA
    fsame_content = fsame_bytes or (ok_files and C.read_gz_lines(C.RAW / ft[0]) == fl514)
    with gzip.open(C.RAW / mx[0], "rt") as fh:
        line = fh.readline()
        while line.startswith("%"):
            line = fh.readline()
        nr, nc, nnz = map(int, line.split())
    exp = int(nps.loc[nps.gsm == gsm, "n"].iloc[0])
    per.append(dict(gsm=gsm, barcodes=nbc, expected=exp, mtx_rows=nr, mtx_cols=nc, nnz=nnz,
                    features_sha_equal=fsame_bytes, features_content_equal=fsame_content,
                    ok=ok_files and nbc == exp and nc == nbc and nr == len(fl514) and fsame_content))
per = pd.DataFrame(per)
per.to_csv(C.OUT / "gate_G0_per_sample.tsv", sep="\t", index=False)
add("per_sample_barcodes_match", int(per.ok.sum()), 48, per.ok.all())
add("total_barcodes", int(per.barcodes.sum()), 283890, per.barcodes.sum() == 283890)
add("features_identical_to_FL514", int(per.features_content_equal.sum()), 48, per.features_content_equal.all())
g = pd.DataFrame(rows).rename(columns={"pass_": "pass"})
g.to_csv(C.OUT / "gate_G0_download.tsv", sep="\t", index=False)
print(g.to_string(index=False))
if not g["pass"].all():
    (C.OUT / "STOP_G0.txt").write_text("untestable: deposited matrices do not match the metadata\n")
    sys.exit("STOP G0: untestable: deposited matrices do not match the metadata")
print("G0 pass")
