#!/usr/bin/env bash
# End-to-end check from a COLD START, with NO external data.
#
# Everything is built here from a seeded synthetic fixture whose ground truth we planted,
# including the cases the method must REFUSE. The checks can fail: three of them flip a
# single threshold and require the verdict to change, so a criterion that was never
# consulted would be caught.
#
# Needs numpy and a python3. h5py, scipy and torch are all optional.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK="${PRETRAINING_EXPOSURE_WORK:-$(mktemp -d "${TMPDIR:-/tmp}/pexp.XXXXXX")}"
export PYTHONNOUSERSITE=1

# ---- find a python with numpy --------------------------------------------------------
PY=""
for c in "${PRETRAINING_EXPOSURE_PY:-}" python3 \
         "$HERE/../../executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python"; do
    [ -n "$c" ] || continue
    if command -v "$c" >/dev/null 2>&1 || [ -x "$c" ]; then
        if "$c" -c "import numpy" >/dev/null 2>&1; then PY="$c"; break; fi
    fi
done
[ -n "$PY" ] || { echo "FAIL: no python3 with numpy found. Set \$PRETRAINING_EXPOSURE_PY."; exit 1; }

VOCAB="${EXPOSURE_VOCABULARY:-$HERE/../encoder-benchmark-v1/exposure.py}"
[ -f "$VOCAB" ] || {
    echo "FAIL: the exposure state vocabulary is not at $VOCAB."
    echo "      This tool IMPORTS the vocabulary rather than restating it. Point"
    echo "      \$EXPOSURE_VOCABULARY at an exposure.py, or place one beside this file."
    exit 1; }
export EXPOSURE_VOCABULARY="$VOCAB"

echo "python   $PY"
echo "vocab    $VOCAB"
echo "workdir  $WORK"
echo

# ---- 1. build the fixtures -----------------------------------------------------------
echo "== 1. build three synthetic fixtures with planted ground truth =="
for s in verbatim crosspipeline nothing_planted; do
    "$PY" -s -E "$HERE/make_fixture.py" --scenario "$s" --out "$WORK/$s"
done
echo

# ---- 2. scan --------------------------------------------------------------------------
run () {   # run <scenario> <outfile> [extra args...]
    local s="$1" o="$2"; shift 2
    "$PY" -s -E "$HERE/exposure_scan.py" \
        --query "$WORK/$s/query.tsv" --cohorts "$WORK/$s/cohorts.tsv" \
        --corpus "$WORK/$s/corpus.npy" --corpus-genes "$WORK/$s/corpus_genes.txt" \
        --out "$o" --device numpy --n-background 400 --model-id "fixture_$s" \
        --quiet "$@"
}

echo "== 2. scan each fixture, streaming the corpus in chunks =="
run verbatim        "$WORK/verbatim.json"        --chunk-rows 512 \
                    --emit-matches "$WORK/verbatim_matches_512.npz"
run crosspipeline   "$WORK/crosspipeline.json"   --chunk-rows 512
run nothing_planted "$WORK/nothing.json"         --chunk-rows 512
echo "   scanned 3 fixtures"
echo

echo "== 3. streaming invariance and backend agreement =="
run verbatim "$WORK/verbatim_c256.json" --chunk-rows 256 \
             --emit-matches "$WORK/verbatim_matches_256.npz"
run verbatim "$WORK/verbatim_c700.json" --chunk-rows 700 \
             --emit-matches "$WORK/verbatim_matches_700.npz"
if "$PY" -c "import torch" >/dev/null 2>&1; then
    "$PY" -s -E "$HERE/exposure_scan.py" \
        --query "$WORK/verbatim/query.tsv" --cohorts "$WORK/verbatim/cohorts.tsv" \
        --corpus "$WORK/verbatim/corpus.npy" --corpus-genes "$WORK/verbatim/corpus_genes.txt" \
        --out "$WORK/verbatim_torch.json" --device torch-cpu --n-background 400 \
        --chunk-rows 512 --model-id fixture_verbatim --quiet
    echo "   torch-cpu backend also run"
else
    echo "   torch absent: backend-agreement check will report SKIP"
fi
echo

# ---- 3b. the same corpus in other containers must give the same verdicts --------------
echo "== 3b. re-container the corpus and require identical verdicts =="
"$PY" -s -E - "$WORK" <<'PYEOF'
import gzip, sys
from pathlib import Path
import numpy as np

d = Path(sys.argv[1]) / "crosspipeline"
X = np.load(d / "corpus.npy")
genes = [l.strip() for l in open(d / "corpus_genes.txt")]
with gzip.open(d / "corpus.tsv.gz", "wt") as fh:
    fh.write("row\t" + "\t".join(genes) + "\n")
    for i, r in enumerate(X):
        fh.write(f"r{i}\t" + "\t".join(f"{v:.7g}" for v in r) + "\n")
made = ["tsv.gz"]
try:
    import h5py
except ImportError:
    h5py = None
if h5py is not None:
    # AnnData dense: gene ids in var/_index as a CATEGORICAL group, not a string array.
    with h5py.File(d / "corpus_dense.h5ad", "w") as f:
        f.create_dataset("X", data=X)
        g = f.create_group("var").create_group("_index")
        g.create_dataset("categories", data=np.array(genes, dtype="S16"))
        g.create_dataset("codes", data=np.arange(len(genes), dtype="i4"))
    # AnnData CSR: an exact round trip, so the reader is isolated from any content change.
    n, m = X.shape
    with h5py.File(d / "corpus_csr.h5ad", "w") as f:
        g = f.create_group("X")
        g.attrs["shape"] = np.array([n, m])
        g.create_dataset("data", data=X.ravel().astype("float32"))
        g.create_dataset("indices", data=np.tile(np.arange(m, dtype="i4"), n))
        g.create_dataset("indptr", data=np.arange(0, n * m + 1, m, dtype="i8"))
        f.create_group("var").create_dataset("_index", data=np.array(genes, dtype="S16"))
    made += ["dense h5ad (categorical var)", "CSR h5ad"]
print("   built: " + ", ".join(made))
PYEOF
CQ=(--query "$WORK/crosspipeline/query.tsv" --cohorts "$WORK/crosspipeline/cohorts.tsv"
    --device numpy --n-background 400 --chunk-rows 512 --quiet)
"$PY" -s -E "$HERE/exposure_scan.py" "${CQ[@]}" --corpus "$WORK/crosspipeline/corpus.tsv.gz" \
      --corpus-genes "$WORK/crosspipeline/corpus_genes.txt" --out "$WORK/cont_tsv.json"
if "$PY" -c "import h5py" >/dev/null 2>&1; then
    "$PY" -s -E "$HERE/exposure_scan.py" "${CQ[@]}" \
          --corpus "$WORK/crosspipeline/corpus_dense.h5ad" --out "$WORK/cont_h5.json"
    "$PY" -s -E "$HERE/exposure_scan.py" "${CQ[@]}" \
          --corpus "$WORK/crosspipeline/corpus_csr.h5ad" --out "$WORK/cont_csr.json"
else
    echo "   h5py absent: the two h5ad checks will report SKIP"
fi
echo

# ---- 4. the three threshold flips: prove each leg is load-bearing ---------------------
echo "== 4. flip one threshold at a time; the verdict must change =="
run verbatim "$WORK/flip_window.json" --chunk-rows 512 --window-over-n-max 1e9
run verbatim "$WORK/flip_band.json"   --chunk-rows 512 --band-lower 0.0
run verbatim "$WORK/flip_inj.json"    --chunk-rows 512 --injectivity-min 1.01
echo "   3 single-threshold reruns done"
echo

# ---- 5. assertions --------------------------------------------------------------------
echo "== 5. assertions =="
"$PY" -s -E - "$HERE" "$WORK" <<'PYEOF'
import json, sys, itertools
from pathlib import Path
import numpy as np

here, work = Path(sys.argv[1]), Path(sys.argv[2])
sys.path.insert(0, str(here))
import exposure_scan as X

exp = json.load(open(here / "expected.json"))
load = lambda p: json.load(open(work / p))
V, C, N = load("verbatim.json"), load("crosspipeline.json"), load("nothing.json")
truthV = json.load(open(work / "verbatim/truth.json"))
truthC = json.load(open(work / "crosspipeline/truth.json"))

fails, n = [], itertools.count(1)


def check(name, ok, detail=""):
    i = next(n)
    print(f"  [{'PASS' if ok else 'FAIL'}] {i:>2}. {name}" + (f"  --  {detail}" if detail else ""))
    if not ok:
        fails.append(name)


def skip(name, detail=""):
    print(f"  [SKIP] {next(n):>2}. {name}" + (f"  --  {detail}" if detail else ""))


# --- the planted positive: verbatim contiguous copy ----------------------------------
for tag, R, T in (("verbatim", V, truthV), ("crosspipeline", C, truthC)):
    c = R["cohorts"]["planted_contiguous"]
    st = R["exposure_state_measured"]["planted_contiguous"]
    want = T["cohorts"]["planted_contiguous"]["constant_offset"]
    check(f"{tag}: planted contiguous copy is encoder_seen",
          st == "encoder_seen", f"state={st}")
    check(f"{tag}: exact corpus offset recovered",
          c["constant_offset"] == want, f"got {c['constant_offset']}, planted {want}")
    check(f"{tag}: every sample on that offset",
          c["largest_aligned_run"] == c["n"],
          f"{c['largest_aligned_run']}/{c['n']}")
    check(f"{tag}: 0 exceedances in {c['n_permutations']} permutation draws",
          c["n_exceedances"] == 0, f"{c['n_exceedances']} exceedances")
    b = R["block_confirmation"]["planted_contiguous"]
    check(f"{tag}: mutual argmax over both rows and columns",
          b["mutual_argmax"] == b["n_query"], f"{b['mutual_argmax']}/{b['n_query']}")
    check(f"{tag}: diagonal r exceeds off-diagonal r",
          b["diag_mean"] > b["offdiag_mean"],
          f"diag {b['diag_mean']:.4f} vs offdiag {b['offdiag_mean']:.4f}")

# --- the false-positive guard: the more important half -------------------------------
for tag, R in (("verbatim", V), ("crosspipeline", C)):
    for coh in ("similar_only", "absent"):
        st = R["exposure_state_measured"][coh]
        check(f"{tag}: '{coh}' is NOT called encoder_seen",
              st != "encoder_seen", f"state={st}")

# --- the two legs are not redundant ---------------------------------------------------
sc = V["cohorts"]["planted_scattered_rows"]
check("scattered copy: injectivity PASSES but the window FAILS",
      sc["distinct_fraction"] >= V["criterion"]["injectivity_min"]
      and sc["window_width_over_n"] > V["criterion"]["window_over_n_max"]
      and not sc["leg1_structure_passes"] and sc["leg2_magnitude_passes"],
      f"injectivity {sc['distinct_fraction']}, window {sc['window_width_over_n']}x")
check("scattered copy is therefore NOT encoder_seen",
      V["exposure_state_measured"]["planted_scattered_rows"] != "encoder_seen")

ps = V["cohorts"]["paired_sibling_scrambled"]
check("paired sibling: structure PASSES but magnitude FAILS",
      ps["leg1_structure_passes"] and not ps["leg2_magnitude_passes"],
      f"injectivity {ps['distinct_fraction']}, window {ps['window_width_over_n']}x, "
      f"median r {ps['max_correlation_median']:.4f} vs band "
      f"{V['same_sample_band']['lower_bound_applied']:.4f}")
check("paired sibling is therefore NOT encoder_seen",
      V["exposure_state_measured"]["paired_sibling_scrambled"] != "encoder_seen")

sh = V["cohorts"]["planted_shuffled_order"]
check("order-scrambled copy is still found (leg 1 is order-free)",
      V["exposure_state_measured"]["planted_shuffled_order"] == "encoder_seen"
      and sh["aligned_fraction"] < 0.5 and not sh["calibrates_band"],
      f"aligned run {sh['largest_aligned_run']}/{sh['n']}, judged on legs 1+2")

# --- flipping one threshold must change the verdict ----------------------------------
fw, fb, fi = load("flip_window.json"), load("flip_band.json"), load("flip_inj.json")
check("removing the window cap FLIPS the scattered copy to encoder_seen",
      fw["exposure_state_measured"]["planted_scattered_rows"] == "encoder_seen",
      "the window test is what rejected it")
check("removing the band FLIPS the paired sibling to encoder_seen",
      fb["exposure_state_measured"]["paired_sibling_scrambled"] == "encoder_seen",
      "the magnitude test is what rejected it")
check("an unreachable injectivity floor calls NOTHING",
      not any(s == "encoder_seen" for s in fi["exposure_state_measured"].values()),
      "the criterion is consulted, not decorative")

# --- the permutation null behaves -----------------------------------------------------
c = V["cohorts"]["planted_contiguous"]
a, m = c["analytic_expected_max_run"], c["perm_null_mean_run"]
check("permutation null mean sits near the analytic expectation",
      abs(m - a) <= exp["null_mean_vs_analytic_abs_tolerance"],
      f"empirical {m:.3f} vs analytic {a:.3f}, |diff| {abs(m - a):.3f}")
check("the observed aligned run exceeds the null mean by a wide margin",
      c["largest_aligned_run"] > 5 * m, f"{c['largest_aligned_run']} vs {m:.3f}")
ab = V["cohorts"]["absent"]
check("the null is not degenerate: an absent cohort does NOT clear it",
      ab["n_exceedances"] > 0.05 * ab["n_permutations"],
      f"{ab['n_exceedances']}/{ab['n_permutations']} exceedances")
perm_keys = [k for k in V["cohorts"]["planted_contiguous"] if k.startswith("perm_null")]
check("no null MAXIMUM is reported as a bar",
      all("max" not in k for k in perm_keys),
      f"a null max grows with the draw count; reported instead: {perm_keys}")

# --- cross-pipeline noise at the measured magnitude ----------------------------------
band = C["same_sample_band"]
lo, hi = exp["reference_same_sample_band"]["min"], exp["reference_same_sample_band"]["max"]
check("cross-pipeline fixture band matches the measured magnitude",
      lo <= band["median"] <= hi,
      f"fixture median {band['median']:.4f} inside reference [{lo}, {hi}]")
check("detection survives that attenuation",
      C["exposure_state_measured"]["planted_contiguous"] == "encoder_seen"
      and C["exposure_state_measured"]["planted_shuffled_order"] == "encoder_seen")

# --- refusal when the band cannot be calibrated ---------------------------------------
check("nothing planted -> band is NOT CALIBRABLE, not a default",
      N["band_calibration"] == "not_calibrable" and N["same_sample_band"] is None,
      f"band_calibration={N['band_calibration']}")
check("nothing planted -> no cohort is called",
      not any(s == "encoder_seen" for s in N["exposure_state_measured"].values()))

# --- streaming ------------------------------------------------------------------------
c512, c256, c700 = V, load("verbatim_c256.json"), load("verbatim_c700.json")
check("chunk count matches ceil(n_rows / chunk_rows): the corpus really streamed",
      c512["corpus"]["chunks_read"] == -(-c512["corpus"]["n_rows"] // 512)
      and c256["corpus"]["chunks_read"] == -(-c256["corpus"]["n_rows"] // 256)
      and c700["corpus"]["chunks_read"] == -(-c700["corpus"]["n_rows"] // 700),
      f"{c256['corpus']['chunks_read']} / {c512['corpus']['chunks_read']} / "
      f"{c700['corpus']['chunks_read']} chunks at 256 / 512 / 700 rows")
z = [np.load(work / f"verbatim_matches_{k}.npz") for k in (256, 512, 700)]
rows = [k for k in z[0].files if k.endswith("__corpus_row")]
rs = [k for k in z[0].files if k.endswith("__max_r")]
# The matched ROW must be bit-identical: it is the answer. The correlation itself is a
# BLAS reduction whose blocking depends on the chunk width, so it is checked to 1e-6.
same_rows = all(np.array_equal(z[0][k], z[i][k]) for k in rows for i in (1, 2))
same_r = max(float(np.abs(z[0][k] - z[i][k]).max()) for k in rs for i in (1, 2))
check("matched corpus rows are IDENTICAL across chunk sizes", same_rows)
check("max correlations agree across chunk sizes to 1e-6", same_r < 1e-6,
      f"largest disagreement {same_r:.2e}")

# --- container formats ------------------------------------------------------------------
for tag, label in (("tsv", "gzipped TSV"), ("h5", "dense h5ad, categorical var index"),
                   ("csr", "AnnData CSR h5ad")):
    p = work / f"cont_{tag}.json"
    name = f"{label} corpus gives identical verdicts"
    if not p.exists():
        skip(name, "h5py not importable")
        continue
    K = json.load(open(p))
    check(name, K["exposure_state_measured"] == C["exposure_state_measured"]
          and K["corpus"]["n_rows"] == C["corpus"]["n_rows"], K["corpus"]["source"].split("/")[-1])

# --- backend --------------------------------------------------------------------------
tp = work / "verbatim_torch.json"
if tp.exists():
    T = json.load(open(tp))
    check("torch-cpu backend agrees with numpy on every state",
          T["exposure_state_measured"] == V["exposure_state_measured"]
          and T["backend"] != V["backend"],
          f"{V['backend']} vs {T['backend']}")
else:
    skip("torch-cpu backend agrees with numpy", "torch not importable")

# --- the vocabulary is imported, never restated --------------------------------------
vocab, vpath = X.load_exposure_vocabulary()
check("the state vocabulary is imported from a real file",
      Path(vpath).is_file() and "encoder_seen" in vocab.CONFOUNDED_STATES, vpath)
check("exposure_scan.py defines no state sets of its own",
      not hasattr(X, "CLEAN_STATES") and not hasattr(X, "CONFOUNDED_STATES"),
      "a convention copied instead of imported is a defect no reading can find")
check("`unknown` is NEVER folded into clean",
      vocab.disposition("unknown") == "unresolved"
      and all(v["clean"] == [] for v in (V["ledger"], C["ledger"], N["ledger"]))
      and all(d != "clean" for R in (V, C, N) for d in R["disposition"].values()),
      "a corpus with no accessions cannot yield clean_declared")
check("every unknown cohort is reported unresolved, none dropped",
      all(sorted(R["ledger"]["clean"] + R["ledger"]["confounded"]
                 + R["ledger"]["unresolved"]) == sorted(R["exposure_state_measured"])
          for R in (V, C, N)))
check("a negative reading travels with every output",
      all("no near-duplicate" in R["negative_reading"] for R in (V, C, N)))

# --- every planted expectation, in one sweep -----------------------------------------
for tag, R, T in (("verbatim", V, truthV), ("crosspipeline", C, truthC),
                  ("nothing_planted", N, json.load(open(work / "nothing_planted/truth.json")))):
    got = {k: R["exposure_state_measured"][k] for k in T["cohorts"]}
    want = {k: T["cohorts"][k]["expected_state"] for k in T["cohorts"]}
    check(f"{tag}: all {len(want)} planted expectations hold", got == want,
          "" if got == want else f"differs at {[k for k in want if got[k] != want[k]]}")

print()
if fails:
    print(f"{len(fails)} CHECK(S) FAILED: {fails}")
    raise SystemExit(1)
print("ALL CHECKS PASSED")
PYEOF

echo
if [ -z "${PRETRAINING_EXPOSURE_WORK:-}" ]; then rm -rf "$WORK"; fi
