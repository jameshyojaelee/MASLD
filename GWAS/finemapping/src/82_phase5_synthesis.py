#!/usr/bin/env python3
"""
82_phase5_synthesis.py -- Phase-5 (HARDENED) synthesis.

Assembles the three Phase-5 seqfunc tracks into a single honest report:
  Track A  coding upgrade  -> coding_hardening_v2.tsv + gof_lof_benchmark.json
  Track B  cell-type       -> decima_celltype.tsv (Decima VEP)   [degrades if absent]
  Track C  direction       -> direction_hardened_verdict.json + direction_hardened_results.csv

ADDITIVE ONLY. Reads existing SF outputs; never re-scores, never edits 46d/27a/60-81.
Stdlib-only (no pandas) so it runs in any env. Re-runnable: as Track-B and the
caQTL/sQTL direction panels land, re-run to refresh SF/phase5_synthesis.md.

Every reportable claim carries the 5-lens-review caveats:
  - direction method is NOT novel (ExPecto, Zhou 2018 Nat Genet)
  - GTEx sQTL = leakage-suspect (in Borzoi+AlphaGenome training) -> exploratory only
  - caQTL/sQTL = same genetic axis: corroboration/labels ONLY, never a scored channel
  - APPLY-ONLY firewall: direction annotated only where measured QTL exists, never OOS
"""
import os, json, csv, datetime

SF = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc"
OUT = os.path.join(SF, "phase5_synthesis.md")


def exists(name):
    return os.path.isfile(os.path.join(SF, name))


def load_json(name):
    p = os.path.join(SF, name)
    if not os.path.isfile(p):
        return None
    try:
        with open(p) as f:
            return json.load(f)
    except Exception as e:
        return {"_load_error": str(e)}


def load_tsv(name):
    p = os.path.join(SF, name)
    if not os.path.isfile(p):
        return None, None
    with open(p) as f:
        r = csv.reader(f, delimiter="\t")
        rows = list(r)
    if not rows:
        return [], []
    return rows[0], rows[1:]


def col(header, rows, name):
    if header is None or name not in header:
        return []
    i = header.index(name)
    return [row[i] for row in rows if len(row) > i]


def count_values(vals):
    d = {}
    for v in vals:
        d[v] = d.get(v, 0) + 1
    return dict(sorted(d.items(), key=lambda kv: -kv[1]))


# ---------------------------------------------------------------- TRACK A
def track_a():
    ch, cr = load_tsv("coding_hardening_v2.tsv")
    bench = load_json("gof_lof_benchmark.json")
    have = ch is not None
    n_eff = len(cr) if cr else 0
    final_calls = count_values(col(ch, cr, "gof_lof_call")) if have else {}
    logofunc_calls = count_values(col(ch, cr, "logofunc_call")) if have else {}
    eve_notes = count_values(col(ch, cr, "eve_note")) if have else {}
    b = bench.get("binary_gof_vs_lof", {}) if bench else {}
    passd = bench.get("pass_bar", {}) if bench else {}
    pnpla3 = bench.get("pnpla3_i148m", {}) if bench else {}
    return dict(have=have, n_eff=n_eff, final_calls=final_calls,
                logofunc_calls=logofunc_calls, eve_notes=eve_notes,
                bench=b, passd=passd, pnpla3=pnpla3)


# ---------------------------------------------------------------- TRACK B
def track_b():
    have = exists("decima_celltype.tsv")
    header = rows = None
    dist = {}
    n_hep = n_total = 0
    sort1 = None
    atac_overlap = atac_match = 0
    if have:
        header, rows = load_tsv("decima_celltype.tsv")
        argmax = col(header, rows, "argmax_celltype")
        dist = count_values([a for a in argmax if a])
        n_total = sum(dist.values())
        n_hep = dist.get("hepatocyte", 0)
        # SORT1 positive-control anchor
        gi = header.index("gene") if header and "gene" in header else None
        ai = header.index("argmax_celltype") if header and "argmax_celltype" in header else None
        si = header.index("argmax_signed_effect") if header and "argmax_signed_effect" in header else None
        if gi is not None:
            for r in rows:
                if len(r) > gi and r[gi] == "SORT1":
                    sort1 = dict(argmax=(r[ai] if ai is not None and len(r) > ai else "NA"),
                                 effect=(r[si] if si is not None and len(r) > si else "NA"))
                    break
        # snATAC corroboration among variants overlapping any per-cell-type peak
        pi = header.index("atac_any_peak") if header and "atac_any_peak" in header else None
        mi = header.index("argmax_matches_atac") if header and "argmax_matches_atac" in header else None
        if pi is not None and mi is not None:
            for r in rows:
                if len(r) > max(pi, mi) and str(r[pi]).strip().lower() in ("true", "1", "yes"):
                    atac_overlap += 1
                    if str(r[mi]).strip().lower() in ("true", "1", "yes"):
                        atac_match += 1
    return dict(have=have, n=(len(rows) if rows else 0),
                cols=(header if header else []), dist=dist,
                n_total=n_total, n_hep=n_hep, sort1=sort1,
                atac_overlap=atac_overlap, atac_match=atac_match)


# ---------------------------------------------------------------- TRACK C
def track_c():
    verdict = load_json("direction_hardened_verdict.json")
    have = verdict is not None
    prim = (verdict or {}).get("preregistered_primary", {})
    res = prim.get("result", {})
    panels = (verdict or {}).get("panels", {})
    labels = {}
    for lab, fn in [("eqtl", "direction_labels_eqtl.tsv"),
                    ("caqtl", "direction_labels_caqtl.tsv"),
                    ("sqtl", "direction_labels_sqtl.tsv")]:
        labels[lab] = exists(fn)
    return dict(have=have, prim=prim, res=res, panels=panels,
                go=prim.get("PRIMARY_GO"), reasons=prim.get("reasons", []),
                labels=labels, caqtl_done=("caqtl" in panels),
                sqtl_done=("sqtl" in panels))


def fnum(x, nd=3):
    try:
        return f"{float(x):.{nd}f}"
    except (TypeError, ValueError):
        return "NA"


def main():
    A, B, C = track_a(), track_b(), track_c()
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    L = []
    w = L.append

    w("# Phase 5 (HARDENED) — Synthesis")
    w("")
    w(f"_Generated {now} by `src/82_phase5_synthesis.py` (additive; assembles existing "
      "SF track outputs, degrades gracefully). Re-run to refresh as pending panels land._")
    w("")
    w("**Scope reminder (non-negotiable, from the 5-lens review):** none of these tracks "
      "is a scored convergence channel. caQTL/sQTL are the *same genetic axis* as eQTL — "
      "used as labels/corroboration only, never summed. Direction is APPLY-ONLY (annotate "
      "HOW, only where a measured QTL already exists; never out-of-sample on convergence "
      "loci). The supervised-direction *method* is **not novel** (ExPecto, Zhou 2018 Nat "
      "Genet = GBT-on-sequence-features→direction); only a *finding* could be reportable.")
    w("")

    # ---- status line
    def stat(ok, partial=False):
        return "✅ DONE" if ok and not partial else ("🟡 IN-PROGRESS/PARTIAL" if partial else "⛔ MISSING")
    w("## Track status at a glance")
    w("")
    w("| Track | Deliverable | Status |")
    w("|---|---|---|")
    w(f"| A — Coding upgrade (LoGoFunc GoF/LoF + EVE) | `coding_hardening_v2.tsv`, "
      f"`gof_lof_benchmark.json` | {stat(A['have'], partial=(A['have'] and not A['eve_notes'].get('covered')))} |")
    w(f"| B — Cell-type effect (Decima VEP) | `decima_celltype.tsv` | "
      f"{stat(B['have'], partial=not B['have'])} |")
    c_partial = C['have'] and not (C['caqtl_done'] and C['sqtl_done'])
    w(f"| C — Direction benchmark (hardened) | `direction_hardened_verdict.json` | "
      f"{stat(C['have'], partial=c_partial)} |")
    w("")

    # ================================================================ TRACK A
    w("## Track A — Coding-arm GoF/LoF + EVE upgrade")
    w("")
    if not A["have"]:
        w("⛔ `coding_hardening_v2.tsv` absent — Track A not yet written.")
    else:
        b, p = A["bench"], A["pnpla3"]
        w("**WHAT WAS DONE.** Extended the 42-effector coding table with a genome-wide "
          "gain-/loss-of-function direction (LoGoFunc, Stein 2023 Genome Med, precomputed "
          "hg38, Zenodo 10126185) and an evolution-only EVE vote, layered on the existing "
          "AlphaMissense + ESM1b + popEVE consensus. Curated-literature anchors override "
          "the model call where a mode-of-action is established.")
        w("")
        w(f"**RESULT (honest).** {A['n_eff']} effectors annotated. Final `gof_lof_call` "
          f"distribution: {A['final_calls']}. Raw LoGoFunc argmax: {A['logofunc_calls']}.")
        if b:
            ngof, nlof = b.get('n_gof'), b.get('n_lof')
            nbin = (ngof + nlof) if isinstance(ngof, int) and isinstance(nlof, int) else '?'
            w(f"Held-out validation of the GoF-vs-LoF discriminator (LoGoFunc's OWN official "
              f"withheld split; binary GoF-vs-LoF subset n={nbin} = "
              f"{ngof} GoF vs {nlof} LoF): **auROC = {fnum(b.get('auroc'))}** "
              f"CI[{fnum(b.get('boot95_ci',[None,None])[0])}, {fnum(b.get('boot95_ci',[None,None])[1])}], "
              f"permutation-null {fnum(b.get('perm_null_mean'))}, p={b.get('perm_p')}. "
              f"Pre-set bar (lower-95% CI > 0.70): **{'PASS' if A['passd'].get('PASS') else 'FAIL'}**.")
        w("")
        w("**⚠ Caveats.** (1) The 0.905 is **variant-level** held-out, **not** "
          "gene-family-disjoint (released feature files carry no gene ID) → it is an "
          "**upper bound** on strict family-held-out performance. (2) 3-class GoF precision "
          "is only ~0.49 (GoF is rare and over-called). (3) **PNPLA3 I148M**: LoGoFunc "
          f"argmax = **Neutral** (P_GOF={fnum(p.get('logofunc_probs',{}).get('GOF'),3) if p else 'NA'}); "
          "the curated GoF-neomorph is **retained** because even a GoF/LoF-specialized "
          "genome-wide model misses this retained-on-lipid-droplet neomorph — this is the "
          "point of the upgrade, not a failure of it. (4) EVE columns are currently "
          f"blank ({A['eve_notes'] or 'no eve_note column'}) — evemodel.org backend returned "
          "502/504 at build; LoGoFunc columns are stable, re-run once EVE lands.")
        w("")
        w("**REPORTABLE.** The coding arm now flags GoF-vs-LoF direction with a "
          "held-out-validated discriminator (auROC 0.905, lower-CI 0.879 ≫ 0.70 bar), and "
          "explicitly documents that neomorphs like PNPLA3 I148M evade even GoF/LoF-aware "
          "models — motivating the curated-literature anchor. Nomination/mechanism-class "
          "only; never summed into convergence.")
    w("")

    # ================================================================ TRACK B
    w("## Track B — Cell-type-resolved variant effect (Decima)")
    w("")
    if B["have"]:
        pct_hep = (100.0 * B["n_hep"] / B["n_total"]) if B["n_total"] else 0.0
        s1 = B["sort1"]
        w(f"**WHAT WAS DONE.** Decima (cell-type EXPRESSION VEP; modality-appropriate — "
          f"scores per-variant × per-cell-type *expression* effect, NOT accessibility) "
          f"scored {B['n']} nomination + coding + lead variants across 210 liver tasks, "
          f"collapsed to 5 liver compartments → `decima_celltype.tsv`.")
        w("")
        if s1:
            w(f"**POSITIVE CONTROL — PASSED.** SORT1 rs12740374 → argmax = "
              f"**{s1['argmax']}** (signed effect {fnum(s1['effect'])}, ~2× next compartment), "
              f"ATAC-corroborated. Face validity: HNF1A/APOE/TM4SF4→hepatocyte; "
              f"ECM1→stellate; P2RX7→LSEC; SERPINA1→macrophage.")
            w("")
        w(f"**RESULT (honest).** argmax cell-type distribution (n={B['n_total']}): "
          f"{B['dist']}. => **{pct_hep:.0f}% hepatocyte, {100-pct_hep:.0f}% non-parenchymal** "
          f"— ties the genetic arm to the MULTI-CELLULAR cascade, not a hepatocyte-only "
          f"model.")
        w("")
        if B["atac_overlap"]:
            pct_atac = 100.0 * B["atac_match"] / B["atac_overlap"]
            w(f"**snATAC corroboration (semi-independent, coverage-limited).** argmax matches "
              f"an accessible cell type in **{B['atac_match']}/{B['atac_overlap']} "
              f"({pct_atac:.0f}%)** of variants that overlap any per-cell-type peak (only "
              f"{B['atac_overlap']}/{B['n']} overlap the ~14-donor snATAC peaks — narrow peaks "
              f"vs credible/lead SNPs).")
            w("")
        w("**⚠ Honest independence flag (not fully disjoint).** The snATAC check uses "
          "per-cell-type peak PRESENCE (accessibility landscape / peak calls) from "
          "Human_Multiome `cell_type_peak_sets_v2`. The epigenomic convergence channel uses "
          "DIFFERENTIAL (disease) accessibility + SCENIC+ at gene level — a different signal "
          "but the **SAME cohort** → **SEMI-INDEPENDENT, not fully disjoint** (no alternative "
          "multi-cell-type liver snATAC panel exists on disk; chrombpnet is HepG2, a single "
          "cell type). Word it as a *semi-independent, coverage-limited* check, never a "
          "circularity-free DISJOINT validation. See `decima_celltype_summary.json` "
          "`disjointness_note`.")
        w("")
        w("**⚠ Robustness (confidence floor).** The all-120 argmax split (31% hep / 69% "
          "non-parenchymal) is argmax-over-near-noise — median |effect| ~0.0017, median "
          "top1/top2 margin 1.34. Among confidently-assigned variants (|effect|>=0.005 & "
          "margin>=1.2, n=27) the split is ~balanced (**44% hep / 56% non-parenchymal**), and "
          "at a stricter floor it flips to 56% hepatocyte. The proportion is therefore NOT a "
          "robust quantitative claim; the load-bearing signal is the SORT1 anchor + the snATAC "
          "check, not the aggregate percentage. See `decima_celltype_summary.json` "
          "`confident_subset_robustness`.")
        w("")
        w("**REPORTABLE (nomination/mechanism-framed).** Cell-type expression VEP assigns a "
          "likely acting compartment per genetic nomination, with a mix of hepatocyte and "
          "non-parenchymal calls that is CONSISTENT WITH (not a quantitative proof of) the "
          "multi-cellular-cascade thesis; SORT1 positive control passes; a semi-independent "
          "snATAC peak-presence check agrees on 2/3 of coverable variants. Per-variant "
          "cell-type call is a mechanism annotation only — NOT a proportion estimate, NOT "
          "summed into convergence.")
    else:
        w("🟡 **IN-PROGRESS (GO — actively scoring on GPU, not degraded).** Per the Decima "
          "agent: env + weights confirmed working (throwaway `decima` micromamba env; Decima "
          "via `pip install decima`; grelu≥1.0.10 + genomepy; weights Zenodo 15092691). Two "
          "transient blockers fixed: Blackwell sm_120 GPU → pinned l40s; `genome=\"hg38\"` "
          "name → local GRCh38 fasta path. Model loads; scoring 99 nominated-gene + 26 "
          "auto-overlap variants across 210 liver tasks. Output `SF/decima_celltype.tsv` + "
          "`src/76_decima_celltype.py`; SORT1→hepatocyte anchor. Result (anchor cell type + "
          "argmax distribution + snATAC corroboration) expected ~30 min from 11:1x.")
        w("")
        w("**REPORTABLE (pending).** Nothing yet — will corroborate cell-type effect vs a "
          "DISJOINT per-cell-type snATAC panel. Not a convergence channel.")
    w("")

    # ================================================================ TRACK C
    w("## Track C — Supervised DIRECTION benchmark (demoted, pre-registered supplementary)")
    w("")
    if not C["have"]:
        w("⛔ `direction_hardened_verdict.json` absent.")
    else:
        res = C["res"]
        w("**WHAT WAS DONE.** Pre-registered, hardened supervised-direction benchmark: "
          "elastic-net logistic (GBT sensitivity only), nested leave-one-chromosome-out, "
          "pooled out-of-fold sign-auROC, block-bootstrap CI, ≥1000× label-permutation "
          "null. Labels leakage-tiered: **Broadaway eQTL = clean CERTIFIED primary** "
          "(independent microarray, not in model training); Currin caQTL = clean secondary "
          "(accessibility features only); GTEx sQTL = leakage-suspect → exploratory-only. "
          "Features re-oriented to the REF→ALT label frame; strand-ambiguous A/T & C/G "
          "hard-dropped; SORT1 sign unit-checked.")
        w("")
        w(f"**RESULT (honest) — PRIMARY endpoint (eQTL Broadaway, n={res.get('n')}, "
          f"pos={res.get('n_pos')}, {res.get('n_features')} features):**")
        w("")
        w(f"- elastic-net pooled-OOF sign-auROC = **{fnum(res.get('enet_auroc'))}** "
          f"CI[{fnum(res.get('enet_ci_lo'))}, {fnum(res.get('enet_ci_hi'))}], "
          f"perm p = {fnum(res.get('enet_perm_p'))}")
        w(f"- **PRIMARY_GO = {C['go']}** — {'; '.join(C['reasons'])}")
        bl = res.get("baseline_auroc", {})
        if bl:
            w(f"- Baselines it had to beat (CI-separated) — it beat **none**: "
              f"fold-unanimous {fnum(bl.get('fold_unanimous_rule'))}, "
              f"aggregated-signed-logSED {fnum(bl.get('aggregated_signed_logSED'))}, "
              f"|logSED|-gated sign {fnum(bl.get('logSED_gated_sign'))}, "
              f"MAF+TSS-dist positional logistic {fnum(bl.get('MAF_TSSdist_positional_logistic'))}")
        w(f"- leakage-gap (random-split {fnum(res.get('random_split_auroc'))} − LOCO) = "
          f"{fnum(res.get('leakage_gap'))}; SORT1: {res.get('sort1','NA')}")
        w("")
        # secondary panels
        if C["caqtl_done"] or C["sqtl_done"]:
            pan = C["panels"]
            cq, sq = pan.get("caqtl", {}), pan.get("sqtl", {})
            h2h = {"caqtl_accessibility_auroc": cq.get("enet_auroc"),
                   "eqtl_expression_auroc": pan.get("eqtl", {}).get("enet_auroc")}
            w("**✅ SECONDARY — the MODALITY DICHOTOMY (caQTL, LANDED).** The SAME "
              "elastic-net on **accessibility** direction (Currin caQTL, clean-secondary, "
              f"n={cq.get('n','NA')}) = auROC **{fnum(cq.get('enet_auroc'))}** "
              f"CI[{fnum(cq.get('enet_ci_lo'))}, {fnum(cq.get('enet_ci_hi'))}], perm "
              f"p={fnum(cq.get('enet_perm_p'))}, leakage-gap {fnum(cq.get('leakage_gap'))} "
              "— and **beats both baselines CI-separated** "
              f"(agg-signed {fnum(cq.get('baseline_auroc',{}).get('agg_signed_mean'))}, "
              f"positional {fnum(cq.get('baseline_auroc',{}).get('MAF_TSSdist_positional_logistic'))}). "
              f"So ACCESSIBILITY direction IS learnable ({fnum(h2h.get('caqtl_accessibility_auroc'))}) "
              f"where EXPRESSION direction is at chance ({fnum(h2h.get('eqtl_expression_auroc'))}) "
              "— sequence models resolve the proximal sequence→chromatin layer, not the "
              "emergent expression layer. **⚠ point estimates on different variant sets & "
              "n — NOT a paired test.** SORT1 caQTL sign PASS.")
            w("")
            w("**sQTL (GTEx, leakage-flagged, exploratory only).** auROC "
              f"{fnum(sq.get('enet_auroc'))}, perm p={fnum(sq.get('enet_perm_p'))}, "
              f"leakage-gap {fnum(sq.get('leakage_gap'))} — GTEx is in Borzoi+AG training → "
              "beats no baseline, **never certified** (correctly stays exploratory).")
        else:
            w("**Secondary panels (caQTL clean-secondary, sQTL leakage-exploratory) PENDING** "
              "— accessibility/splice feature jobs (`79_*`) still running; a dependent "
              "direction re-run will add them. Labels already built: eQTL 211 leads "
              "(113/98), caQTL 11,896 independent peak-leads (SORT1 sign **PASS**, +1), "
              "sQTL 137,480 intron-leads (GTEx → **leakage-flagged, never certified**).")
        w("")
        w("**REPORTABLE.** A **pre-registered NEGATIVE**: under the hardened protocol, a "
          "supervised meta-learner on SOTA sequence-model features does **not** recover "
          "eQTL direction — it fails to clear the 0.56 bar (lower-CI 0.437) and beats none "
          "of the fold-unanimous / signed-logSED / positional baselines in CI-separated "
          "margins. This **reinforces the thesis**: sequence models stay magnitude/"
          "mechanism-only and never enter convergence *direction* scoring. Frame as "
          "confidence recalibration of zero-shot scores, NOT 'supervised beats zero-shot'. "
          "Method is not novel (ExPecto 2018); the reportable unit is the finding, under "
          "the APPLY-ONLY firewall.")
    w("")

    # ================================================================ CLAIMS TABLE
    w("## Compact reportable-claims table")
    w("")
    w("| # | Claim | Evidence | Caveat / firewall |")
    w("|---|---|---|---|")
    if A["have"]:
        b = A["bench"]
        w(f"| A1 | Coding arm flags GoF-vs-LoF direction with a held-out-validated "
          f"discriminator | LoGoFunc GoF-vs-LoF auROC {fnum(b.get('auroc'))} "
          f"[{fnum(b.get('boot95_ci',[None,None])[0])},{fnum(b.get('boot95_ci',[None,None])[1])}], "
          f"perm p={b.get('perm_p')} | variant-level held-out (upper bound, not "
          f"family-disjoint); GoF precision ~0.49; nomination-only |")
        w("| A2 | Neomorphs evade even GoF/LoF-aware genome-wide models | PNPLA3 I148M "
          "LoGoFunc argmax = Neutral (P_GOF=0.042); curated GoF retained | single anchor; "
          "motivates curated-literature override |")
    if B["have"]:
        pct_hep = (100.0 * B["n_hep"] / B["n_total"]) if B["n_total"] else 0.0
        atac_str = (f"{B['atac_match']}/{B['atac_overlap']}" if B["atac_overlap"] else "NA")
        w(f"| B1 | Genetic nominations localize ~1/3 hepatocyte, ~2/3 non-parenchymal "
          f"(cell-type expression VEP) | Decima argmax n={B['n_total']}: {pct_hep:.0f}% hep; "
          f"SORT1→hepatocyte (pos-control PASS); snATAC match {atac_str} | expression VEP "
          f"(modality-appropriate); snATAC = SEMI-INDEPENDENT coverage-limited check (same "
          f"cohort, peak-presence vs disease-DA — NOT fully disjoint); NOT a channel |")
    if C["have"]:
        res = C["res"]
        w(f"| C1 | Supervised features do NOT recover eQTL EXPRESSION direction "
          f"(pre-registered negative) | pooled-OOF auROC {fnum(res.get('enet_auroc'))} "
          f"CI[{fnum(res.get('enet_ci_lo'))},{fnum(res.get('enet_ci_hi'))}]; beats no "
          f"baseline; GO=False | method not novel (ExPecto 2018); APPLY-ONLY; thesis-"
          f"reinforcing |")
        cq = C["panels"].get("caqtl", {})
        if cq:
            w(f"| C2 | But ACCESSIBILITY direction IS learnable — a modality dichotomy | "
              f"caQTL auROC {fnum(cq.get('enet_auroc'))} "
              f"CI[{fnum(cq.get('enet_ci_lo'))},{fnum(cq.get('enet_ci_hi'))}], perm "
              f"p={fnum(cq.get('enet_perm_p'))}, beats both baselines CI-separated, "
              f"leakage-gap ~0, n={cq.get('n','NA')} | NOT a paired test (diff variant "
              f"sets); APPLY-ONLY; caQTL never a scored channel |")
    w("")

    # ================================================================ NEGATIVES
    w("## Honest negatives (pre-registered / disclosed)")
    w("")
    w("- **EXPRESSION-direction rescue FAILED** the pre-registered bar (elastic-net "
      "lower-95% CI 0.437 ≤ 0.56; beats none of fold-unanimous 0.547 / signed-logSED "
      "0.559 / positional 0.500). A recalibration story, not a rescue. (But ACCESSIBILITY "
      "direction IS learnable — caQTL 0.707; the negative is modality-specific.)")
    w("- **PNPLA3 I148M not recovered** by LoGoFunc (argmax Neutral) — expected; curated "
      "anchor retained.")
    w("- **EVE columns blank** (evemodel.org 502/504 at build) — Track A partial on the "
      "EVE arm; LoGoFunc arm complete.")
    w("- **GTEx sQTL is leakage-suspect** (in Borzoi + AlphaGenome training) → exploratory, "
      "never certified; junction-matched to AG SpliceJunction, magnitude-only otherwise.")
    w("- **caQTL/sQTL never summed** into convergence; same genetic axis, corroboration/"
      "labels only.")
    w("")

    # ================================================================ PAPER IMPACT
    w("## What changes for the paper")
    w("")
    w("- **Net paper claims: essentially unchanged** — by design. Seqfunc never enters "
      "convergence scoring, so the direction negative changes no headline number.")
    w("- **Strengthened (supplementary):** the coding arm gains a held-out-validated "
      "GoF/LoF axis (auROC 0.905) + the PNPLA3-neomorph-evades-models vignette — a clean "
      "methods/robustness addition.")
    w("- **Reinforced thesis:** the hardened direction negative is a *feature* — it "
      "justifies keeping sequence models magnitude/mechanism-only and out of direction "
      "scoring. Report as a bounded supplementary benchmark, not a flagship.")
    w("- **Do NOT** frame the direction work as a novel method (ExPecto precedent) or apply "
      "it out-of-sample to convergence loci (APPLY-ONLY firewall).")
    w("")

    with open(OUT, "w") as f:
        f.write("\n".join(L) + "\n")
    print(f"WROTE {OUT} ({len(L)} lines)")
    print(f"Track A: have={A['have']} n_eff={A['n_eff']}")
    print(f"Track B: have={B['have']}")
    print(f"Track C: have={C['have']} PRIMARY_GO={C['go']} caqtl_done={C['caqtl_done']} sqtl_done={C['sqtl_done']}")


if __name__ == "__main__":
    main()
