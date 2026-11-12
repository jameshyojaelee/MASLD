#!/usr/bin/env python3
"""Stage D summaries and candidate panels (PRESPEC_AG2X_stageD.md sections 4-6). No hosted call.

Per variant:
  tracks: observed-track mean of predicted REF and ALT signal, +/- 1 kb; local change per condition =
    mean over observed tracks of log2(sum ALT) - log2(sum REF) over the 501 bp centred on the variant.
  ISM (REF, 64 bp; ALT background, 32 bp): score per substitution = mean over observed tracks of the scorer raw
    score; position importance = mean of its 3 substitution scores; rank of the variant position by |importance|;
    top-decile positions (7 of 64, 4 of 32) and the motif hits (A3 PWMs, A3 thresholds) that contain them;
    for each motif gained (lost) at the variant, ratio of mean |importance| inside the best covering hit to the
    mean outside it, in ALT-background (REF) hepatocyte ISM.
  check: REF-ISM entry for the variant's own substitution against the Stage B raw score of the same readout.
Panels: one PDF per (variant, readout); fonttype 42; 6-pt text; tracks: ALT filled #C9265E, REF line #9E9E9E, +/-1 kb and
+/-150 bp files; ISM colour scale blue #1565C0 to magenta #C9265E. CAPTION_<variant>.md beside them. Every table row carries candidate_not_adopted = True.
"""
from __future__ import annotations

import json
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import ag2x_lib as L  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

D = L.PKG / "stageD"
RAW = D / "raw"
TAB = D / "tables"
PAN = D / "panels"
B = L.PKG / "stageB"
REF_C, ALT_C, CTRL_C = "#1565C0", "#C9265E", "#9E9E9E"
FS = 6
plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "Helvetica", "font.size": FS,
                     "axes.labelsize": FS, "xtick.labelsize": FS, "ytick.labelsize": FS, "axes.linewidth": 0.4,
                     "xtick.major.width": 0.4, "ytick.major.width": 0.4, "axes.spines.top": False,
                     "axes.spines.right": False, "savefig.dpi": 300})
CMAP = LinearSegmentedColormap.from_list("bwm", [REF_C, "#FFFFFF", ALT_C])
TRACK_PANELS = [("bulk_liver", "ATAC"), ("hepatocyte", "SC_ATAC"), ("monocyte_ctrl", "SC_ATAC")]
LOCAL_READOUTS = [("bulk_liver", "ATAC"), ("bulk_liver", "DNASE"), ("hepatocyte", "SC_ATAC"), ("monocyte_ctrl", "SC_ATAC")]
ISM_READOUTS = [("ism_ref", "bulk_liver", "ATAC"), ("ism_ref", "bulk_liver", "DNASE"), ("ism_ref", "hepatocyte", "SC_ATAC"),
                ("ism_alt", "hepatocyte", "SC_ATAC")]
ISM_PANELS = [("ism_ref", "bulk_liver", "ATAC"), ("ism_ref", "hepatocyte", "SC_ATAC"), ("ism_alt", "hepatocyte", "SC_ATAC")]
HALF = {"ism_ref": 32, "ism_alt": 16}
LOCAL_W = 501
STAGEB_COL = {("bulk_liver", "ATAC"): "bulk_liver|ATAC", ("bulk_liver", "DNASE"): "bulk_liver|DNASE",
              ("hepatocyte", "SC_ATAC"): "hepatocyte|SC_ATAC", ("monocyte_ctrl", "SC_ATAC"): "monocyte_ctrl|SC_ATAC"}


def safe(u):
    return u.replace(":", "_")


def flip(u):
    c, p, r, a = u.split(":")
    return f"{c}:{p}:{a}:{r}"


def load_tracks(u, cname):
    f = RAW / "tracks" / f"{safe(u)}__{cname}.npz"
    if not f.exists():
        return None, None
    return np.load(f), json.load(open(f.with_suffix(".json")))


def local_change(u, cname, otype):
    z, meta = load_tracks(u, cname)
    if z is None or f"{otype}_ref" not in z:
        return np.nan, 0, None
    m = meta["outputs"][otype]
    res, vb = m["resolution"], m["variant_bin"]
    h = (LOCAL_W // 2) // res
    ref, alt = z[f"{otype}_ref"], z[f"{otype}_alt"]
    sr = ref[vb - h:vb + h + 1].sum(axis=0)
    sa = alt[vb - h:vb + h + 1].sum(axis=0)
    vals = np.log2(sa) - np.log2(sr)
    return float(np.mean(vals)), int(ref.shape[1]), meta


def load_ism(u, step, cname, otype):
    fs = sorted((RAW / step).glob(f"{safe(u)}__{cname}__*.parquet")) if (RAW / step).exists() else []
    if not fs:
        return None
    x = pd.concat([pd.read_parquet(f) for f in fs], ignore_index=True)
    x = x[(x.output_type == otype) & (x.imputation_mode == "observed")]
    g = x.groupby("variant_id").raw_score.mean().reset_index()
    vv = g.variant_id.str.extract(r"^(chr[^:]+):(\d+):([ACGT]+)>([ACGT]+)$")
    g["pos"] = vv[1].astype(int)
    g["ref_base"] = vv[2]
    g["alt_base"] = vv[3]
    g["n_tracks"] = x.groupby("variant_id").track_name.nunique().reindex(g.variant_id).to_numpy()
    return g


def ism_matrix(g, u, step):
    c, p, r, a = L.parse_uid(u)
    h = HALF[step]
    start1 = p - h + 1  # 1-based first position
    n = 2 * h
    M = np.full((4, n), np.nan)
    for _, x in g.iterrows():
        j = x.pos - start1
        if 0 <= j < n:
            M["ACGT".index(x.alt_base), j] = x.raw_score
    off = h - 1
    return M, off


def runs(idx):
    idx = sorted(idx)
    out, cur = [], [idx[0]] if idx else []
    for i in idx[1:]:
        if i == cur[-1] + 1:
            cur.append(i)
        else:
            out.append(cur)
            cur = [i]
    if cur:
        out.append(cur)
    return out


def ism_summary(M, off, seq, hits, allele, top_k):
    imp = np.nanmean(M, axis=0)
    order = np.argsort(-np.abs(imp))
    rank_var = int(np.where(order == off)[0][0]) + 1
    top = sorted(order[:top_k].tolist())
    h = hits[hits.allele == allele]
    desc = []
    for rr in runs(top):
        s_rel, e_rel = rr[0] - off, rr[-1] - off
        ov = h[(h.start_rel <= e_rel) & (h.end_rel >= s_rel)]
        fams = ";".join(sorted(set(ov.name.astype(str)))[:8]) if len(ov) else "no motif hit"
        desc.append(f"{s_rel:+d}..{e_rel:+d} ({seq[rr[0]:rr[-1] + 1]}; mean importance {imp[rr].mean():+.3f}; motifs: {fams})")
    return imp, rank_var, desc


def hit_ratio(imp, off, hits, allele, motif_id, n):
    h = hits[(hits.allele == allele) & (hits.motif_id == motif_id) & hits.covers_variant]
    if not len(h):
        return np.nan, ""
    b = h.sort_values("score", ascending=False).iloc[0]
    s, e = int(b.start_rel) + off, int(b.end_rel) + off
    s, e = max(s, 0), min(e, n - 1)
    inside = np.zeros(n, bool)
    inside[s:e + 1] = True
    a = np.abs(imp)
    if inside.all() or not inside.any():
        return np.nan, f"{b.start_rel:+d}..{b.end_rel:+d}"
    return float(a[inside].mean() / a[~inside].mean()), f"{int(b.start_rel):+d}..{int(b.end_rel):+d} ({b.strand})"


def plot_tracks(u, cname, otype, out, half_bp=1000):
    """ALT drawn as a filled magenta profile, REF as a gray (#9E9E9E) line on top, so both stay visible."""
    z, meta = load_tracks(u, cname)
    if z is None:
        return False
    m = meta["outputs"][otype]
    res, vb = m["resolution"], m["variant_bin"]
    ref, alt = z[f"{otype}_ref"].mean(axis=1), z[f"{otype}_alt"].mean(axis=1)
    k = half_bp // res
    xs = (np.arange(vb - k, vb + k + 1) - vb) * res
    r_, a_ = ref[vb - k:vb + k + 1], alt[vb - k:vb + k + 1]
    fig, ax = plt.subplots(figsize=(3.0, 1.1))
    ax.axvline(0, color="#E0E0E0", lw=0.4, zorder=0)
    ax.fill_between(xs, 0, a_, color=ALT_C, alpha=0.3, lw=0, zorder=1)
    ax.plot(xs, a_, color=ALT_C, lw=0.4, zorder=2)
    ax.plot(xs, r_, color=CTRL_C, lw=0.6, zorder=3)
    ax.set_xlim(-half_bp, half_bp)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("position relative to variant (bp)")
    ax.set_ylabel(f"predicted {otype.replace('_', ' ')}")
    fig.tight_layout(pad=0.3)
    fig.savefig(out, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    return True


def plot_ism(M, off, seq, alt_base, out):
    n = M.shape[1]
    v = np.nanmax(np.abs(M))
    fig, ax = plt.subplots(figsize=(0.07 * n + 0.8, 0.9))
    im = ax.imshow(np.ma.masked_invalid(M), cmap=CMAP, vmin=-v, vmax=v, aspect="auto", interpolation="nearest")
    ax.set_yticks(range(4))
    ax.set_yticklabels(list("ACGT"))
    ax.set_xticks(range(n))
    ax.set_xticklabels(list(seq))
    ax.get_xticklabels()[off].set_color(ALT_C)
    ax.tick_params(length=1, pad=1)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.01)
    cb.ax.tick_params(labelsize=FS, length=1)
    cb.set_label("predicted change (log2)", fontsize=FS)
    cb.outline.set_linewidth(0.3)
    fig.tight_layout(pad=0.3)
    fig.savefig(out, bbox_inches="tight", pad_inches=0.02)  # keeps the full colour-bar label inside the page
    plt.close(fig)


def main():
    TAB.mkdir(parents=True, exist_ok=True)
    PAN.mkdir(parents=True, exist_ok=True)
    v = pd.read_csv(D / "inputs/stageD_variants.tsv", sep="\t")
    me = pd.read_csv(D / "measured/measured_at_variants.tsv", sep="\t").set_index("uid")
    mc = pd.read_csv(D / "motifs/motif_changes_by_cluster.tsv", sep="\t")
    hits = pd.read_csv(D / "motifs/motif_hits_64bp.tsv.gz", sep="\t")
    seqs = pd.read_csv(D / "motifs/sequences_64bp.tsv", sep="\t").set_index("uid")
    vf = pd.read_csv(B / "tables/stageB_variant_fields.tsv.gz", sep="\t").set_index("uid")
    loc_rows, ism_rows, motif_rows, check_rows = [], [], [], []
    for _, r in v.iterrows():
        u = r.uid
        c, p, ref_b, alt_b = L.parse_uid(u)
        gene = str(r.genes_stageB).split(",")[0].strip() if isinstance(r.genes_stageB, str) else ""
        folder = PAN / f"{r.locus_id}_{gene}".replace(" ", "")
        folder.mkdir(parents=True, exist_ok=True)
        hu = hits[hits.uid == u]
        refseq = seqs.at[u, "ref_seq"]
        altseq = seqs.at[u, "alt_seq"]
        # local predicted change
        for cname, ot in LOCAL_READOUTS:
            val, ntr, meta = local_change(u, cname, ot)
            col = STAGEB_COL[(cname, ot)]
            loc_rows.append({"uid": u, "condition": cname, "output": ot, "local_log2_alt_minus_ref_501bp": val,
                             "n_observed_tracks": ntr, "context_bp": meta["context_bp"] if meta else np.nan,
                             "stageB_raw": vf.at[u, f"raw:{col}"] if u in vf.index else np.nan,
                             "stageB_percentile": vf.at[u, f"pct:{col}"] if u in vf.index else np.nan,
                             "predicted": True, "candidate_not_adopted": True})
        # tracks panels
        for cname, ot in TRACK_PANELS:
            plot_tracks(u, cname, ot, folder / f"{safe(u)}_tracks_{cname}_{ot}.pdf")
            plot_tracks(u, cname, ot, folder / f"{safe(u)}_tracks_{cname}_{ot}_zoom150.pdf", half_bp=150)
        # ISM
        for step, cname, ot in ISM_READOUTS:
            g = load_ism(u, step, cname, ot)
            if g is None or not len(g):
                continue
            M, off = ism_matrix(g, u, step)
            h = HALF[step]
            seq = (refseq if step == "ism_ref" else altseq)[32 - h:32 + h]
            hh = hu.copy()
            allele = "ref" if step == "ism_ref" else "alt"
            imp, rank_var, desc = ism_summary(M, off, seq, hh, allele, 7 if step == "ism_ref" else 4)
            own = np.nan
            if step == "ism_ref":
                own = M["ACGT".index(alt_b), off]
                col = STAGEB_COL.get((cname, ot))
                sb = vf.at[u, f"raw:{col}"] if (col and u in vf.index) else np.nan
                check_rows.append({"uid": u, "readout": f"{cname}|{ot}", "ism_own_substitution": own, "stageB_raw": sb,
                                   "abs_diff": abs(own - sb) if np.isfinite(sb) else np.nan,
                                   "candidate_not_adopted": True})
            ism_rows.append({"uid": u, "step": step, "condition": cname, "output": ot, "n_substitutions": int(np.isfinite(M).sum()),
                             "n_observed_tracks": int(g.n_tracks.max()), "variant_position_rank_abs_importance": rank_var,
                             "of_positions": M.shape[1], "variant_position_importance": float(imp[off]),
                             "own_substitution_score": own, "max_abs_substitution": float(np.nanmax(np.abs(M))),
                             "top_positions": " | ".join(desc), "predicted": True, "candidate_not_adopted": True})
            if (step, cname, ot) in ISM_PANELS:
                plot_ism(M, off, seq, alt_b, folder / f"{safe(u)}_{step}_{cname}_{ot}.pdf")
            # motif hit ratios for gained / lost motifs (one row per cluster)
            if cname == "hepatocyte" and ot == "SC_ATAC":
                for _, m in mc[mc.uid == u].iterrows():
                    want = "alt" if m.direction == "gained" else "ref"
                    if allele != want:
                        continue
                    ratio, span = hit_ratio(imp, off, hh, allele, m.motif_id, M.shape[1])
                    motif_rows.append({"uid": u, "motif_id": m.motif_id, "name": m["name"], "cluster_family": m.cluster_family,
                                       "direction": m.direction, "delta": m.delta_alt_minus_ref,
                                       "tf_expressed_hepatocyte": m.tf_expressed_hepatocyte, "ism": step,
                                       "hit_span_rel": span, "importance_ratio_inside_vs_outside": ratio,
                                       "candidate_not_adopted": True})
    loc = pd.DataFrame(loc_rows)
    ism = pd.DataFrame(ism_rows)
    mot = pd.DataFrame(motif_rows)
    chk = pd.DataFrame(check_rows)
    loc.to_csv(TAB / "stageD_local_change.tsv", sep="\t", index=False)
    ism.to_csv(TAB / "stageD_ism_summary.tsv", sep="\t", index=False)
    mot.to_csv(TAB / "stageD_motif_ism_overlap.tsv", sep="\t", index=False)
    chk.to_csv(TAB / "stageD_ism_vs_stageB_check.tsv", sep="\t", index=False)
    # sign agreement with measured data
    sg = []
    for _, r in v.iterrows():
        u = r.uid
        bl = loc[(loc.uid == u) & (loc.condition == "bulk_liver") & (loc.output == "ATAC")].local_log2_alt_minus_ref_501bp
        hp = loc[(loc.uid == u) & (loc.condition == "hepatocyte") & (loc.output == "SC_ATAC")].local_log2_alt_minus_ref_501bp
        bl = float(bl.iloc[0]) if len(bl) else np.nan
        hp = float(hp.iloc[0]) if len(hp) else np.nan
        ca = str(me.at[u, "currin_caqtl_beta_alt_values"]) if u in me.index else "nan"
        ca_vals = [float(x) for x in ca.split(";") if x not in ("", "nan")]
        ai = me.at[u, "hep_ai_pooled_log2_alt_ref"] if (u in me.index and "hep_ai_pooled_log2_alt_ref" in me) else np.nan
        sg.append({"uid": u, "pred_bulk_liver_atac": bl, "caqtl_beta_alt": ";".join(f"{x:.3f}" for x in ca_vals),
                   "caqtl_sign_agrees": ("not measured" if not ca_vals else
                                         ("yes" if all(np.sign(x) == np.sign(bl) for x in ca_vals) else
                                          "no" if all(np.sign(x) != np.sign(bl) for x in ca_vals) else "mixed")),
                   "pred_hepatocyte_scatac": hp, "hep_ai_pooled_log2_alt_ref": ai,
                   "hep_ai_sign_agrees": ("not measured" if not np.isfinite(ai) else
                                          ("yes" if np.sign(ai) == np.sign(hp) else "no")),
                   "candidate_not_adopted": True})
    sg = pd.DataFrame(sg)
    sg.to_csv(TAB / "stageD_measured_sign_agreement.tsv", sep="\t", index=False)
    L.log(f"local {len(loc)} rows; ism {len(ism)}; motif overlap {len(mot)}; check max abs diff "
          f"{chk.abs_diff.max() if len(chk) else 'NA'}")


if __name__ == "__main__":
    main()
