#!/usr/bin/env python3
"""
90_ag_member_finemechanism.py — Phase-6 WS-3 Tier-2: AlphaGenome ZERO-SHOT fine
regulatory sub-mechanism for the NON-LEAD credible-set members (critique G1 completeness).

Tier-1 (src/89) PIP-weighted the COARSE mechanism class (coding/splice/regulatory) over all
credible-set members from FREE VEP annotation and showed it is robust to credible-set
uncertainty. This Tier-2 step extends that to the FINE regulatory sub-mechanism
(promoter/TSS, enhancer-disrupting, TF-footprint, splice-altering) so we can ask whether the
sub-mechanism is ALSO consistent across a diffuse credible set, or only defined at the lead.

SCOPE: the ~505 diffuse/moderate-CS members that carry an hg38 coordinate and are not already
AG-scored as leads (64 leads live in ag_mechanism_profile.tsv). Concentrated CS (eff #vars ~1)
are already answered by their single lead.

MODEL/ToS: AlphaGenome, ONE score_variant call/variant with the 11 scalar liver-masked heads
(reusing the VALIDATED src/72 fingerprint + assign_mechanism_class logic via importlib).
ZERO-SHOT ANNOTATION ONLY — never used to train a model; non-commercial (AlphaGenome FAQ).
The failed contact which-gene arm is NOT run (contact_high=False); TF-footprint requires
corroboration we do not assert here, so it stays conservative.

FRAMING: the fine taxonomy is a HAND-BUILT, UNVALIDATED classifier (critique M6) → the PIP-
weighted fine mechanism is HYPOTHESIS-LEVEL. The load-bearing G1 answer is Tier-1's coarse
posterior (src/89). APPLY-ONLY: annotation only; never a convergence channel; never a direction claim.

NOTE (2026-07-12): an earlier run failed with grpc `Expected str, not tuple` — that was a BUG
HERE (load_api_key returns a (key, source) TUPLE that must be unpacked), NOT env drift. The
alphagenome env is fine. Fixed. Tier-1 (src/89) answers critique G1 at the mechanism-CLASS level
WITHOUT AlphaGenome; this fine-mechanism layer is optional hypothesis-level completeness.
"""
import argparse
import importlib.util
import os
import time

import numpy as np
import pandas as pd

ROOT = os.environ.get("MASLD_PROJECT_ROOT",
                      "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FIN = os.path.join(ROOT, "GWAS/finemapping")
SF = os.path.join(FIN, "results/seqfunc")
SRC = os.path.join(FIN, "src")


def load_ag72():
    spec = importlib.util.spec_from_file_location(
        "ag72", os.path.join(SRC, "72_ag_multimodal_mechanism.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)             # module-level only; alphagenome imported inside main()
    return mod


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def build_targets():
    """Members of diffuse/moderate CS with hg38 coords + the CS's nominated target gene."""
    mem = pd.read_csv(os.path.join(SF, "nominated_cs_members.tsv"), sep="\t")
    percs = pd.read_csv(os.path.join(SF, "pip_weighted_mechanism_per_cs.tsv"), sep="\t")
    keep_cs = set(percs[percs["concentration"].isin(["diffuse", "moderate"])]["cs_key"])
    mem = mem[mem["cs_key"].isin(keep_cs)].copy()
    sub = pd.read_csv(os.path.join(SF, "variant_substrate_hg38.tsv"), sep="\t")
    sub = sub[["variant_id_hg19", "chr", "pos_hg38", "ref_hg38", "alt_hg38", "gene", "ensembl"]]
    m = mem.merge(sub, left_on="variant_id", right_on="variant_id_hg19", how="left")
    m = m[m["pos_hg38"].notna()].copy()
    # CS target gene = the nominated best gene at the locus (from nominations), else substrate gene
    nom = pd.read_csv(os.path.join(SF, "borzoi_magnitude_nominations.tsv"), sep="\t")
    leads = pd.read_csv(os.path.join(SF, "borzoi_magnitude_leads.tsv"), sep="\t")
    loc_gene = leads.merge(nom[["locus_id", "best_gene", "best_gene_ensembl"]],
                           on="locus_id", how="left")
    # map each CS to a locus target gene via the nominated lead in that CS
    lead2loc = dict(zip(leads["variant_id_hg19"].astype(str), leads["locus_id"].astype(str)))
    loc2gene = dict(zip(loc_gene["locus_id"].astype(str),
                        zip(loc_gene["best_gene"].astype(str), loc_gene["best_gene_ensembl"].astype(str))))
    # unique members (score once): pick the target gene from the member's substrate gene as fallback
    u = m.drop_duplicates("variant_id").copy()
    tg, te = [], []
    for _, r in u.iterrows():
        g, e = str(r.get("gene", "")), str(r.get("ensembl", ""))
        tg.append(g if g and g != "nan" else "")
        te.append(e if e and e != "nan" else "")
    u["target_sym"], u["target_ens"] = tg, te
    return u


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(SF, "ag_member_finemechanism.tsv"))
    ap.add_argument("--key-file", default=os.path.expanduser("~/.alphagenome_key"))
    ap.add_argument("--seq-length", default="SEQUENCE_LENGTH_1MB")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--sleep", type=float, default=0.2)
    ap.add_argument("--max-retries", type=int, default=6)
    args = ap.parse_args()

    ag = load_ag72()
    api_key, key_src = ag.load_api_key(args.key_file)   # returns (key, source) tuple
    if api_key is None:
        raise SystemExit("no AlphaGenome API key (env ALPHA_GENOME_API_KEY or ~/.alphagenome_key)")
    log(f"API key loaded from {key_src} (len={len(api_key)}; not shown)")

    from alphagenome.data import genome
    from alphagenome.models import dna_client, variant_scorers
    import grpc
    seq_len = dna_client.SUPPORTED_SEQUENCE_LENGTHS[args.seq_length]
    R = variant_scorers.RECOMMENDED_VARIANT_SCORERS
    scalar_scorers = [R[k] for k in ag.SCALAR_SCORER_KEYS if k in R]
    model = dna_client.create(api_key)
    log(f"dna_client ready; {len(scalar_scorers)} scalar scorers; seq_len={seq_len:,}")

    def _retry(fn, tag):
        last = None
        for attempt in range(args.max_retries):
            try:
                return fn()
            except grpc.RpcError as e:
                last = e
                if e.code() not in (grpc.StatusCode.UNAVAILABLE, grpc.StatusCode.RESOURCE_EXHAUSTED,
                                    grpc.StatusCode.DEADLINE_EXCEEDED, grpc.StatusCode.INTERNAL):
                    raise
                wait = min(45.0, 2.0 ** attempt)
                log(f"  transient {e.code()} on {tag} (attempt {attempt+1}); backoff {wait:.0f}s")
                time.sleep(wait)
        raise last

    def score_one(chrom, pos1, ref, alt, vid, tens, tsym):
        c = chrom if str(chrom).startswith("chr") else f"chr{chrom}"
        v = genome.Variant(chromosome=c, position=int(pos1), reference_bases=str(ref),
                           alternate_bases=str(alt), name=vid)
        iv = v.reference_interval.resize(seq_len)
        sc = _retry(lambda: model.score_variant(
            interval=iv, variant=v, variant_scorers=scalar_scorers,
            organism=dna_client.Organism.HOMO_SAPIENS), f"fp {vid}")
        df = variant_scorers.tidy_scores(sc)
        fp = ag.extract_fingerprint(df, tens or None, tsym or None)
        # fingerprint-only fine mechanism: no contact (failed arm), no TF corroboration, no chromHMM
        prim, sec, conf, scores = ag.assign_mechanism_class(
            fp, tss_dist=None, contact_high=False, tf_corroborated=False, chromhmm_state="")
        return fp, prim, sec, conf, scores

    # anchors first (format + sanity validation: SORT1 should light accessibility/enhancer)
    rows = []
    for a in ag.ANCHORS:
        try:
            fp, prim, sec, conf, sc = score_one(a["chr"], a["pos_hg38"], a["hg38_ref"],
                                                a["hg38_alt"], a["variant_id"], a["ensembl"], a["gene"])
            log(f"ANCHOR {a['anchor_name']}: fine_mech={prim} (sec={sec}, margin={conf}) scores={ {k:round(v,3) for k,v in sc.items()} }")
            rows.append(dict(variant_id=a["variant_id_hg19"], is_anchor=True, anchor=a["anchor_name"],
                             target_sym=a["gene"], fine_mech=prim, fine_mech_secondary=sec,
                             class_margin=conf, **{f"q_{k}": fp.get(k) for k in
                             ("q_atac", "q_dnase", "q_chip_tf", "q_cage", "q_procap",
                              "q_splice_sites", "q_splice_site_usage", "q_rna_gene")}))
        except Exception as e:
            log(f"ANCHOR {a['anchor_name']} FAILED: {e}")

    tgt = build_targets()
    if args.limit:
        tgt = tgt.head(args.limit)
    log(f"scoring {len(tgt)} unique diffuse/moderate-CS members ...")
    for i, (_, r) in enumerate(tgt.iterrows(), 1):
        vid = str(r["variant_id"])
        try:
            fp, prim, sec, conf, sc = score_one(r["chr"], r["pos_hg38"], r["ref_hg38"],
                                                r["alt_hg38"], vid, r["target_ens"], r["target_sym"])
            rows.append(dict(variant_id=vid, is_anchor=False, anchor="",
                             target_sym=r["target_sym"], fine_mech=prim, fine_mech_secondary=sec,
                             class_margin=conf, **{f"q_{k}": fp.get(k) for k in
                             ("q_atac", "q_dnase", "q_chip_tf", "q_cage", "q_procap",
                              "q_splice_sites", "q_splice_site_usage", "q_rna_gene")}))
        except Exception as e:
            log(f"  member {vid} FAILED: {e}")
            rows.append(dict(variant_id=vid, is_anchor=False, anchor="", target_sym=r["target_sym"],
                             fine_mech="ERROR", fine_mech_secondary="", class_margin=np.nan))
        if i % 25 == 0:
            log(f"  {i}/{len(tgt)} members scored")
        time.sleep(args.sleep)

    out = pd.DataFrame(rows)
    out.to_csv(args.out, sep="\t", index=False)
    log(f"wrote {args.out} ({len(out)} rows; {(out['fine_mech']!='ERROR').sum()} scored OK)")
    log("fine_mech distribution: " + str(out[out['fine_mech'] != 'ERROR']['fine_mech'].value_counts().to_dict()))


if __name__ == "__main__":
    main()
