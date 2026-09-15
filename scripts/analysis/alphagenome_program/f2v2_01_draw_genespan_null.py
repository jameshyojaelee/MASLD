#!/usr/bin/env python3
"""f2-haplotype-v2 step 1: draw the 100 gene-span-matched chr6 null pairs. No network use.

Defect 1 of v1: the six GNMT-region pairs read RNA out over the 3,194 bp GENCODE v49 GNMT span, while all
100 of their matched chr6 null draws read RNA out over a central 20 kb window, because a null draw sits at
a random chr6 position and so never contains the GNMT span. This script draws a null whose RNA readout IS a
gene span, by placing each draw inside a real chr6 gene of matched width and handing score_pair that gene's
Ensembl id. Nothing in f2_recipe.py or in score_pair changes; only which (gene, pair) goes in.

Rules are those stamped in prespec/f2v2_00_prespec.json before this ran.
"""

from __future__ import annotations

import csv
import gzip
import json
import os
import pathlib
import re
import sys

import numpy as np
import pysam

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import f2_recipe as R  # noqa: E402

PROJECT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT = pathlib.Path(os.environ["F2V2_OUT_ROOT"])
V1 = PROJECT / "GWAS/finemapping/results/alphagenome_program/f2-haplotype-20260914T230433Z"
PANEL = V1 / "tables/null_panel/chr6.tsv.gz"
P5_PAIRS = PROJECT / "GWAS/finemapping/results/alphagenome_atlas/p5-haplotypes-20260910T004716Z/tables/haplotype_pairs.tsv"

SEED = 20260915
N_DRAWS = 100
GNMT_ENSEMBL = "ENSG00000124713"
GNMT_SPAN = (42960690, 42963883)          # GENCODE v49, 1-based inclusive; width 3,194 bp
GNMT_PRIMARY = (42961020, 42963523)       # rs2296805 T>G, rs2296804 C>G
WIDTH_LO, WIDTH_HI = 2472, 6388           # 2x GNMT's span at the top; the separation bin's floor at the bottom
HOST_GENE_TYPES = {"protein_coding"}      # GNMT's own gene_type; see prespec/f2v2_00_amendment_01.json
N_OBS = 200


def decile_bins(seps: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    edges = np.quantile(seps, [i / 10 for i in range(11)])
    idx = np.clip(np.searchsorted(edges[1:-1], seps, side="right"), 0, 9)
    return edges, idx


def main() -> None:
    (OUT / "tables").mkdir(parents=True, exist_ok=True)

    # the separation bin must be the SAME bin v1 used, so recompute the edges from the same 200-pair list
    with open(P5_PAIRS) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    rows.sort(key=lambda r: (-float(r["weight_product"]), r["signal_uid"], r["v1"], r["v2"]))
    seps = np.array([int(r["separation_bp"]) for r in rows[:N_OBS]], float)
    edges, _ = decile_bins(seps)
    gnmt_sep = abs(GNMT_PRIMARY[1] - GNMT_PRIMARY[0])
    gbin = int(np.clip(np.searchsorted(edges[1:-1], gnmt_sep, side="right"), 0, 9))
    lo, hi = float(edges[gbin]), float(edges[gbin + 1])
    print(f"decile edges {np.round(edges, 1).tolist()}")
    print(f"GNMT primary separation {gnmt_sep} bp -> bin {gbin} [{lo:.1f}, {hi:.1f}]")

    # ---- host genes: chr6, width matched, not GNMT and not overlapping it
    genes = []
    with gzip.open(R.GTF, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            p = line.split("\t", 9)
            if p[2] != "gene" or p[0] != "chr6":
                continue
            gid = re.search(r'gene_id "([^"]+)"', p[8])
            gtype = re.search(r'gene_type "([^"]+)"', p[8])
            gname = re.search(r'gene_name "([^"]+)"', p[8])
            if not gid:
                continue
            a, b = int(p[3]), int(p[4])
            genes.append({"ensembl": gid.group(1).split(".")[0], "start": a, "end": b, "width": b - a + 1,
                          "gene_type": gtype.group(1) if gtype else "", "gene_name": gname.group(1) if gname else ""})
    print(f"chr6 gene rows in GENCODE v49: {len(genes)}")
    cand = [g for g in genes
            if WIDTH_LO <= g["width"] <= WIDTH_HI
            and g["gene_type"] in HOST_GENE_TYPES
            and g["ensembl"] != GNMT_ENSEMBL
            and not (g["end"] >= GNMT_SPAN[0] and g["start"] <= GNMT_SPAN[1])]
    print(f"{sorted(HOST_GENE_TYPES)}, width in [{WIDTH_LO}, {WIDTH_HI}], not overlapping GNMT: {len(cand)} genes")

    # ---- panel and reference
    pos, al1, al2 = [], [], []
    with gzip.open(PANEL, "rt") as handle:
        handle.readline()
        for line in handle:
            q = line.rstrip("\n").split("\t")
            pos.append(int(q[1])); al1.append(q[2]); al2.append(q[3])
    order = np.argsort(np.asarray(pos))
    pos = np.asarray(pos)[order]; al1 = np.asarray(al1)[order]; al2 = np.asarray(al2)[order]
    print(f"chr6 common-EUR panel: {pos.size} variants")

    fasta = pysam.FastaFile(R.FASTA)
    chrlen = fasta.get_reference_length("chr6")

    # ---- every qualifying (gene, pair): both variants inside the span, separation in the bin,
    #      REF equal to a panel allele, 1 Mb window inside chr6
    hosts = []
    n_pairs_total = 0
    for g in cand:
        i0 = int(np.searchsorted(pos, g["start"], "left"))
        i1 = int(np.searchsorted(pos, g["end"], "right"))
        if i1 - i0 < 2:
            continue
        idx = np.arange(i0, i1)
        pairs = []
        for a in idx:
            j0 = int(np.searchsorted(pos, pos[a] + lo, "left"))
            j1 = int(np.searchsorted(pos, pos[a] + hi, "right"))
            for b in range(max(j0, a + 1), min(j1, i1)):
                x1, x2 = int(pos[a]), int(pos[b])
                if not (lo <= x2 - x1 <= hi):
                    continue
                mid = (x1 + x2) // 2
                if mid - R.WINDOW // 2 < 0 or mid - R.WINDOW // 2 + R.WINDOW > chrlen:
                    continue
                placed, ok = [], True
                for q, x in ((a, x1), (b, x2)):
                    base = fasta.fetch("chr6", x - 1, x).upper()
                    u, v = str(al1[q]), str(al2[q])
                    if base == u:
                        placed.append((x, u, v))
                    elif base == v:
                        placed.append((x, v, u))
                    else:
                        ok = False
                        break
                if ok:
                    pairs.append(placed)
        if pairs:
            hosts.append({"gene": g, "pairs": pairs})
            n_pairs_total += len(pairs)
    print(f"host genes with >=1 qualifying pair: {len(hosts)}; qualifying pairs: {n_pairs_total}")
    if not hosts:
        raise RuntimeError("no host gene carries a qualifying pair")

    # amendment 01: one pair from every host gene, then a second, position-disjoint pair from as many hosts
    # as the draw count still needs, those hosts drawn without replacement.
    rng = np.random.default_rng(SEED)
    first = rng.permutation(len(hosts))
    plan: list[tuple[int, int]] = []          # (host index, pair index)
    used: dict[int, list[int]] = {}
    for h in first:
        if len(plan) >= N_DRAWS:
            break
        k = int(rng.integers(0, len(hosts[int(h)]["pairs"])))
        plan.append((int(h), k))
        used[int(h)] = [k]
    need = N_DRAWS - len(plan)
    if need > 0:
        def disjoint_options(h: int) -> list[int]:
            taken = used[h][0]
            p_taken = hosts[h]["pairs"][taken]
            pos_taken = {p_taken[0][0], p_taken[1][0]}
            return [k for k, pr in enumerate(hosts[h]["pairs"])
                    if k != taken and not ({pr[0][0], pr[1][0]} & pos_taken)]

        eligible = [int(h) for h in first if disjoint_options(int(h))]
        print(f"hosts carrying a position-disjoint second pair: {len(eligible)}; {need} needed")
        if len(eligible) < need:
            raise RuntimeError(f"{need} second pairs needed but only {len(eligible)} hosts can give one")
        for h in rng.permutation(np.asarray(eligible))[:need]:
            h = int(h)
            options = disjoint_options(h)
            k = int(options[int(rng.integers(0, len(options)))])
            plan.append((h, k))
            used[h].append(k)
    draws = []
    for n, (h, k) in enumerate(plan):
        host = hosts[h]
        (p1, r1, a1), (p2, r2, a2) = host["pairs"][k]
        g = host["gene"]
        draws.append({"tag": f"gsnull{n}", "chrom": "chr6", "ensembl": g["ensembl"], "gene_name": g["gene_name"],
                      "gene_type": g["gene_type"], "span_start": g["start"], "span_end": g["end"],
                      "span_width_bp": g["width"], "decile_bin": gbin, "bin_lo": lo, "bin_hi": hi,
                      "pos1": p1, "ref1": r1, "alt1": a1, "pos2": p2, "ref2": r2, "alt2": a2,
                      "separation_bp": p2 - p1,
                      "frac_pos1_in_span": round((p1 - g["start"]) / (g["width"] - 1), 4),
                      "frac_pos2_in_span": round((p2 - g["start"]) / (g["width"] - 1), 4),
                      "n_qualifying_pairs_in_gene": len(host["pairs"]),
                      "pair_index_in_gene": k,
                      "is_first_pair_from_this_gene": bool(used[h][0] == k),
                      "base_tag": "gnmt_primary"})
    cols = list(draws[0])
    with open(OUT / "tables" / "gnmt_genespan_null_pairs.tsv", "w") as handle:
        w = csv.DictWriter(handle, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for d in draws:
            w.writerow(d)

    W = np.array([d["span_width_bp"] for d in draws], float)
    S = np.array([d["separation_bp"] for d in draws], float)
    meta = {"seed": SEED, "n_draws": len(draws), "chromosome": "chr6",
            "separation_bin": {"bin": gbin, "lo": lo, "hi": hi, "gnmt_primary_separation_bp": gnmt_sep},
            "width_band_bp": [WIDTH_LO, WIDTH_HI], "host_gene_types": sorted(HOST_GENE_TYPES),
            "n_distinct_host_genes": len({d["ensembl"] for d in draws}),
            "n_draws_that_are_a_second_pair_from_a_gene_already_drawn": int(sum(1 for d in draws if not d["is_first_pair_from_this_gene"])),
            "gnmt_span": {"start": GNMT_SPAN[0], "end": GNMT_SPAN[1], "width_bp": GNMT_SPAN[1] - GNMT_SPAN[0] + 1},
            "candidate_genes": len(cand), "host_genes_with_a_qualifying_pair": len(hosts),
            "qualifying_pairs_total": n_pairs_total,
            "achieved_span_width_bp": {"min": float(W.min()), "median": float(np.median(W)), "max": float(W.max())},
            "achieved_separation_bp": {"min": float(S.min()), "median": float(np.median(S)), "max": float(S.max())},
            "observed_gnmt_primary_frac_in_span": [round((GNMT_PRIMARY[0] - GNMT_SPAN[0]) / (GNMT_SPAN[1] - GNMT_SPAN[0]), 4),
                                                   round((GNMT_PRIMARY[1] - GNMT_SPAN[0]) / (GNMT_SPAN[1] - GNMT_SPAN[0]), 4)],
            "gene_type_counts": {t: int(sum(1 for d in draws if d["gene_type"] == t))
                                 for t in sorted({d["gene_type"] for d in draws})}}
    json.dump(meta, (OUT / "tables" / "gnmt_genespan_null_design.json").open("w"), indent=1)
    print(json.dumps(meta, indent=1))


if __name__ == "__main__":
    main()
