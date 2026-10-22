"""Run 06_analysis.py on simulated counts to check that it executes and recovers a planted effect. No study data.
Usage: python -I test_analysis_synthetic.py WORK_DIR      (WORK_DIR must not exist)
Planted truth: 14 women in two cohorts; skew in biopsy 1 is 0.8; in biopsy 2 it moves by 0.04 per stage step
(so |Delta| rises with |stage change|), plus one woman with a flipped dominant X."""
import sys, os, gzip, json, subprocess
import numpy as np

work = sys.argv[1]
if os.path.exists(work):
    sys.exit("refusing to overwrite " + work)
S = os.path.dirname(os.path.abspath(__file__))
os.makedirs(os.path.join(work, "inputs")); os.makedirs(os.path.join(work, "counts")); os.makedirs(os.path.join(work, "gate"))
rng = np.random.default_rng(1)
NX, NA = 120, 600
sites = [("chrX", 1000 + 10 * i, "A", "G", "X_inactive" if i % 3 == 0 else "X_other") for i in range(NX)] + \
        [("chr19", 1000 + 10 * i, "A", "G", "AUTO") for i in range(NA)]
with open(os.path.join(work, "inputs", "sites.tsv"), "w") as fh:
    fh.write("chrom\tpos\tref\talt\trsid\tcls\tgenes\tstatus\n")
    for c, p, r, a, k in sites:
        fh.write("%s\t%d\t%s\t%s\trs0\t%s\tG\tx\n" % (c, p, r, a, k))


def counts(n, p):
    n0 = rng.binomial(n, 0.5); r0 = rng.binomial(n0, p); r1 = rng.binomial(n - n0, p)
    return r0, n0 - r0, r1, (n - n0) - r1


rows = []
hdr = ["person_id", "cohort", "set", "include", "lib1", "lib2", "stage1", "stage2", "abs_stage_change"]
gate_rows = []
for w in range(14):
    cohort = "C1" if w < 4 else "C2"
    stage = int(rng.integers(0, 3)); d_stage = int(rng.integers(0, 3))
    th1 = 0.8; th2 = th1 - 0.04 * d_stage if w != 13 else 0.2
    phase = rng.random(NX + NA) < 0.5
    libs = ["L%d_1" % w, "L%d_2" % w]
    for b, lib in enumerate(libs):
        th = th1 if b == 0 else th2
        with gzip.open(os.path.join(work, "counts", lib + ".counts.tsv.gz"), "wt") as fh:
            fh.write("chrom\tpos\tref_p0\talt_p0\tref_p1\talt_p1\tother_p0\tother_p1\n")
            for i, (c, p, r, a, k) in enumerate(sites):
                n = 20 + rng.poisson(45)
                pr = (th if phase[i] else 1 - th) if c == "chrX" else (0.5 + (0.1 if phase[i] else -0.1))
                fh.write("%s\t%d\t%d\t%d\t%d\t%d\t0\t0\n" % ((c, p) + counts(n, pr)))
    rows.append(["P%d" % w, cohort, "primary", "True", libs[0], libs[1], stage, stage + d_stage, d_stage])
    gate_rows.append(("P%d" % w, NX))
with open(os.path.join(work, "inputs", "n2a_pairs.tsv"), "w") as fh:
    fh.write("\t".join(hdr) + "\n")
    for r in rows:
        fh.write("\t".join(map(str, r)) + "\n")
with open(os.path.join(work, "gate", "gate_per_woman.tsv"), "w") as fh:
    fh.write("person_id\tk_x_informative\n")
    for p, k in gate_rows:
        fh.write("%s\t%d\n" % (p, k))
json.dump({"primary": {"gate_pass": True}}, open(os.path.join(work, "gate", "gate_summary.json"), "w"))
r = subprocess.run([sys.executable, "-I", os.path.join(S, "06_analysis.py"), work, os.path.join(work, "gate"), os.path.join(work, "out"), "primary"],
                   capture_output=True, text=True)
print(r.stdout[-4000:]); print(r.stderr[-3000:])
sys.exit(r.returncode)
