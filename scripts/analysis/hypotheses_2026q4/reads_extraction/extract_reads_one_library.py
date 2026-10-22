#!/usr/bin/env python3
"""One pass over one STAR BAM (GRCh38, GENCODE v49) writing small tables for
H30 (FOXO1 hotspots), H35 (A-to-I editing), H36 (junctions), H38 (IGH constant
region) and H3/H31 (allele counts at named variants).

Read filters for editing and allele counts: mapped, primary, not duplicate or QC-fail,
STAR unique (MAPQ 255), base quality >= --min-baseq, base not within --end-trim of
either read end (editing only). The two SERPINA1 variants also accept primary
alignments with MAPQ 3 (NH = 2), because the locus is duplicated on alt contig
KI270847.1 and STAR splits those reads between the copies; both copies are counted.
Overlapping mates are counted once at the named variants (including the HSD17B13
insertion) and the FOXO1 window; at editing sites both mates are counted (as in the Alu
Editing Index). IGH counts use every alignment of a fragment (multi-mappers included).
Junction counts come from the STAR SJ.out.tab of the same run, not from the BAM.

Usage: extract_reads_one_library.py --library L --cohort C --layout paired|single \
  --bam B --sj SJ --rpg RPG --log-final LOG --resources DIR --out-dir DIR
"""
import argparse
import gzip
import hashlib
import json
import os
import shutil
import sys
import time
from array import array
from collections import Counter
from bisect import bisect_left

import numpy as np
import pysam

ap = argparse.ArgumentParser()
ap.add_argument("--library", required=True)
ap.add_argument("--cohort", required=True)
ap.add_argument("--layout", required=True, choices=["paired", "single"])
ap.add_argument("--bam", required=True)
ap.add_argument("--sj", required=True)
ap.add_argument("--rpg", required=True)
ap.add_argument("--log-final", required=True)
ap.add_argument("--resources", required=True)
ap.add_argument("--out-dir", required=True)
ap.add_argument("--threads", type=int, default=2)
ap.add_argument("--min-baseq", type=int, default=20)
ap.add_argument("--end-trim", type=int, default=5)
ap.add_argument("--min-site-cov", type=int, default=3)
ap.add_argument("--max-records", type=int, default=0, help="stop after N records (smoke tests only)")
ap.add_argument("--seed", type=int, default=20260923)
args = ap.parse_args()

t_start = time.time()
SCRIPT_SHA256 = hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest()  # hash of the code that runs, taken at start
if os.path.exists(args.out_dir):
    print(f"already extracted: {args.out_dir}")
    sys.exit(0)
tmp_dir = f"{args.out_dir}.partial-{os.environ.get('SLURM_JOB_ID', os.getpid())}-{os.environ.get('SLURM_ARRAY_TASK_ID', '0')}"
if os.path.exists(tmp_dir):
    shutil.rmtree(tmp_dir)
os.makedirs(tmp_dir)

BASE = {"A": 0, "C": 1, "G": 2, "T": 3}
RES = args.resources


def read_tsv(path):
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(head, line.rstrip("\n").split("\t"))) for line in fh]


# ------------------------------------------------------------------ resources
variants = read_tsv(os.path.join(RES, "variants_hg38.tsv"))
hotspots = read_tsv(os.path.join(RES, "hotspots_foxo1.tsv"))
window = read_tsv(os.path.join(RES, "foxo1_window.tsv"))
ig_rows = read_tsv(os.path.join(RES, "ig_constant_exons.tsv"))
panel = np.load(os.path.join(RES, "editing_panel_hg38.npz"))

af = pysam.AlignmentFile(args.bam, "rb", threads=args.threads, check_sq=False)
tid_of = {n: i for i, n in enumerate(af.references)}

# pileup targets: SNVs (plus their alt-contig copies) and the FOXO1 window (unique positions)
targets = {}  # (chrom, pos1) -> index
target_list = []
variant_targets = {}  # rsid -> target indices (main contig first, then alt-contig copy)
multi_ok = set()  # targets that also accept STAR primary alignments with MAPQ 3 (NH = 2): main/alt contig pairs


def add_target(chrom, pos):
    key = (chrom, pos)
    if key not in targets:
        targets[key] = len(target_list)
        target_list.append(key)
    return targets[key]


for v in variants:
    if v["kind"] == "snv":
        idxs = [add_target(v["chrom"], int(v["pos"]))]
        if v["alt_contig"] != "NA":
            idxs.append(add_target(v["alt_contig"], int(v["alt_pos"])))
            multi_ok.update(idxs)
        variant_targets[v["rsid"]] = idxs
for w in window:
    add_target(w["chrom"], int(w["pos"]))
tg_pos, tg_idx = {}, {}
for (c, p), i in sorted(targets.items(), key=lambda kv: (kv[0][0], kv[0][1])):
    t = tid_of[c]
    tg_pos.setdefault(t, array("i")).append(p)
    tg_idx.setdefault(t, []).append(i)
tg_cnt = array("i", bytes(4 * len(target_list) * 2 * 6))  # [target][frag strand][A,C,G,T,del,other]
tg_seen = [dict() for _ in target_list]  # qname -> 1, to count overlapping mates once
hot_keys = {(h["chrom"], int(h["pos"])) for h in hotspots}
hot_alt = {}  # (pos) -> {base: set of fragment signatures}
hot_all = {}  # pos -> set of fragment signatures (any base)
for c, p in hot_keys:
    hot_alt[p] = {b: set() for b in "ACGT"}
    hot_all[p] = set()

# insertion variant window (rs72613567)
ins_v = [v for v in variants if v["kind"] == "insertion"][0]
INS_CHROM, INS_ANCHOR = ins_v["chrom"], int(ins_v["pos"])
INS_L, INS_R = INS_ANCHOR - 3, INS_ANCHOR + 4  # flanks, 1-based
ref_win = pysam.FastaFile(os.environ.get("MASLD_FASTA", "/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/fasta/genome.fa"))
INS_REF = ref_win.fetch(INS_CHROM, INS_L - 1, INS_R).upper()
INS_ALT = INS_REF[: INS_ANCHOR - INS_L + 1] + "A" + INS_REF[INS_ANCHOR - INS_L + 1:]
ins_tid = tid_of[INS_CHROM]
ins_cnt = {"ref": [0, 0], "alt": [0, 0], "other": [0, 0]}  # by fragment strand
ins_seen = set()  # qnames already classified, so overlapping mates count once

# editing panel
ed_pos, ed_site = {}, {}
for c in [f"chr{i}" for i in range(1, 23)] + ["chrX"]:
    p = panel[f"{c}_pos"]
    if len(p):
        t = tid_of[c]
        ed_pos[t] = array("i", p.astype(np.int32).tobytes())
        ed_site[t] = array("i", panel[f"{c}_site"].astype(np.int32).tobytes())
ed_strand = panel["strand"]
n_sites = len(ed_strand)
ed_nonalu_np = panel["is_nonalu"]
ed_ref, ed_alt, ed_nonalu = (bytes(panel[k].astype(np.uint8).tolist()) for k in ("ref_code", "alt_code", "is_nonalu"))  # fast int indexing
ed_cnt = array("i", bytes(4 * 4 * n_sites))  # [site][ref+, alt+, ref-, alt-] by fragment strand
ed_other = [0, 0]  # bases at panel sites that are neither ref nor alt: [ALU, non-ALU]

# IGH constant region on chr14 and its alt contig KI270846.1 (every constant gene is duplicated there)
ig_by_tid, ig_range = {}, {}
for r_ in ig_rows:
    t_ = tid_of[r_["chrom"]]
    ig_by_tid.setdefault(t_, []).append((int(r_["start"]), int(r_["end"]), r_["gene"]))
for t_, lst in ig_by_tid.items():
    ig_range[t_] = (min(e[0] for e in lst), max(e[1] for e in lst))
ig_frag = {}  # qname -> [genes, his, nh]

# ------------------------------------------------------------------ the pass
n_rec = n_unique = n_nh2 = n_edit_reads = n_target_reads = 0
end_trim, min_bq = args.end_trim, args.min_baseq
it = iter(af)
for r in it:
    n_rec += 1
    if args.max_records and n_rec > args.max_records:
        break
    flag = r.flag
    if flag & 4:
        continue
    tid = r.reference_id
    start = r.reference_start
    end = r.reference_end
    rng = ig_range.get(tid)
    if rng is not None and start < rng[1] and end > rng[0] and not flag & 0x200:
        genes = {g for es, ee, g in ig_by_tid[tid] if es <= end and ee >= start + 1}
        if genes:
            rec = ig_frag.get(r.query_name)
            if rec is None:
                rec = ig_frag[r.query_name] = [set(), set(), r.get_tag("NH") if r.has_tag("NH") else 1]
            rec[0].update(genes)
            rec[1].add(r.get_tag("HI") if r.has_tag("HI") else 1)
    mq = r.mapping_quality
    if flag & 0x700 or (mq != 255 and mq != 3):
        continue
    uniq = mq == 255
    if uniq:
        n_unique += 1
    else:
        n_nh2 += 1
    fs = 1 if bool(flag & 16) != bool(flag & 128) else 0  # strand of read 1 (single-end: the read): 1 = minus
    ep = ed_pos.get(tid)
    tp = tg_pos.get(tid)
    need_edit = need_tg = False
    if ep is not None and uniq:
        i0 = bisect_left(ep, start + 1)
        if i0 < len(ep) and ep[i0] <= end:
            need_edit = True
    if tp is not None:
        t0 = bisect_left(tp, start + 1)
        if t0 < len(tp) and tp[t0] <= end:
            need_tg = True
    need_ins = uniq and tid == ins_tid and start < INS_R and end >= INS_L
    if not (need_edit or need_tg or need_ins):
        continue
    seq = r.query_sequence
    qual = r.query_qualities
    cig = r.cigartuples
    rlen = len(seq)
    # aligned blocks: (ref0, q0, len); deletions: (ref0, len)
    blocks, dels = [], []
    rp, qp = start, 0
    for op, ln in cig:
        if op == 0 or op == 7 or op == 8:
            blocks.append((rp, qp, ln))
            rp += ln
            qp += ln
        elif op == 1 or op == 4:
            qp += ln
        elif op == 2:
            dels.append((rp, ln))
            rp += ln
        elif op == 3:
            rp += ln
    if need_edit:
        n_edit_reads += 1
        sp = ed_site[tid]
        j = i0
        n = len(ep)
        while j < n and ep[j] <= end:
            p0 = ep[j] - 1
            for bs, bq, bl in blocks:
                if bs <= p0 < bs + bl:
                    q = bq + (p0 - bs)
                    if q >= end_trim and q < rlen - end_trim and qual[q] >= min_bq:
                        s = sp[j]
                        b = BASE.get(seq[q], 4)
                        if b == ed_ref[s]:
                            ed_cnt[4 * s + 2 * fs] += 1
                        elif b == ed_alt[s]:
                            ed_cnt[4 * s + 2 * fs + 1] += 1
                        elif b < 4:
                            ed_other[ed_nonalu[s]] += 1
                    break
            j += 1
    if need_tg:
        n_target_reads += 1
        qn = r.query_name
        j = t0
        n = len(tp)
        ti = tg_idx[tid]
        frag_sig = None
        while j < n and tp[j] <= end:
            p1 = tp[j]
            p0 = p1 - 1
            t = ti[j]
            if not uniq and t not in multi_ok:
                j += 1
                continue
            code = None
            for bs, bq, bl in blocks:
                if bs <= p0 < bs + bl:
                    q = bq + (p0 - bs)
                    if qual[q] >= min_bq:
                        code = BASE.get(seq[q], 5)
                    break
            if code is None:
                for ds, dl in dels:
                    if ds <= p0 < ds + dl:
                        code = 4
                        break
            if code is not None and qn not in tg_seen[t]:
                tg_seen[t][qn] = 1
                tg_cnt[(t * 2 + fs) * 6 + code] += 1
                if p1 in hot_all:
                    if frag_sig is None:
                        frag_sig = (start, r.next_reference_start if flag & 1 else end, fs)
                    hot_all[p1].add(frag_sig)
                    if code < 4:
                        hot_alt[p1]["ACGT"[code]].add(frag_sig)
            j += 1
    if need_ins:
        pairs = r.get_aligned_pairs()
        qmap = {rp_: qp_ for qp_, rp_ in pairs if rp_ is not None and qp_ is not None}
        if (INS_L - 1) in qmap and (INS_R - 1) in qmap:
            ql, qr = qmap[INS_L - 1], qmap[INS_R - 1]
            if r.query_name not in ins_seen:
                ins_seen.add(r.query_name)
                local = seq[ql:qr + 1].upper()
                kind = "ref" if local == INS_REF else ("alt" if local == INS_ALT else "other")
                ins_cnt[kind][fs] += 1
af.close()
t_pass = time.time() - t_start

# ------------------------------------------------------------------ library QC inputs
rpg = np.loadtxt(args.rpg, dtype=str, delimiter="\t", usecols=(0, 1, 2, 3)) if os.path.getsize(args.rpg) else None
special = {r[0]: [int(x) for x in r[1:]] for r in rpg if r[0].startswith("N_")}
gene_rows = np.array([[int(x) for x in r[1:]] for r in rpg if not r[0].startswith("N_")])
col2, col3, col4 = (int(gene_rows[:, k].sum()) for k in (0, 1, 2))
frac_sense = col3 / max(col3 + col4, 1)
strand_call = "forward" if frac_sense >= 0.8 else ("reverse" if frac_sense <= 0.2 else "unstranded")
log = {}
with open(args.log_final) as fh:
    for line in fh:
        if "|" in line:
            k, v = line.split("|", 1)
            log[k.strip()] = v.strip()

# ------------------------------------------------------------------ junctions (STAR SJ.out.tab of the same run)
strand_ch = {"0": ".", "1": "+", "2": "-"}
n_j = tot_u = tot_m = 0
with open(args.sj) as fh, gzip.open(os.path.join(tmp_dir, "junctions.tsv.gz"), "wt") as oj:
    oj.write("chrom\tintron_start\tintron_end\tstrand\tmotif\tannotated\tn_unique\tn_multi\tmax_overhang\n")
    for line in fh:
        f = line.rstrip("\n").split("\t")
        oj.write("\t".join([f[0], f[1], f[2], strand_ch[f[3]], f[4], f[5], f[6], f[7], f[8]]) + "\n")
        n_j += 1
        tot_u += int(f[6])
        tot_m += int(f[7])

# ------------------------------------------------------------------ variants.tsv
tg = np.frombuffer(tg_cnt, dtype=np.int32).reshape(len(target_list), 2, 6)


def tcounts(chrom, pos):
    return tg[targets[(chrom, pos)]]


with open(os.path.join(tmp_dir, "variants.tsv"), "w") as fv:
    fv.write("library\trsid\tgene\tchrom\tpos\tref\talt\tn_ref_plus\tn_ref_minus\tn_alt_plus\tn_alt_minus\tn_other\tn_del\tdepth\talt_fraction\tdepth_on_alt_contig\n")
    for v in variants:
        d_alt_contig = 0
        if v["kind"] == "snv":
            idxs = variant_targets[v["rsid"]]
            c = tg[idxs].sum(axis=0)  # [strand][A,C,G,T,del,other], main contig plus alt-contig copy
            rc, ac = BASE[v["ref"]], BASE[v["alt"]]
            n_other = int(sum(c[:, k].sum() for k in range(4) if k not in (rc, ac)))
            nr, na = int(c[:, rc].sum()), int(c[:, ac].sum())
            nd = int(c[:, 4].sum())
            row = [int(c[0, rc]), int(c[1, rc]), int(c[0, ac]), int(c[1, ac]), n_other, nd]
            if len(idxs) > 1:
                ca = tg[idxs[1:]].sum(axis=0)
                d_alt_contig = int(sum(ca[:, k].sum() for k in range(4)))
        else:
            nr, na = sum(ins_cnt["ref"]), sum(ins_cnt["alt"])
            n_other = sum(ins_cnt["other"])
            row = [ins_cnt["ref"][0], ins_cnt["ref"][1], ins_cnt["alt"][0], ins_cnt["alt"][1], n_other, 0]
        depth = nr + na + n_other
        fv.write("\t".join(map(str, [args.library, v["rsid"], v["gene"], v["chrom"], v["pos"], v["ref"], v["alt"]] + row
                          + [depth, f"{na / depth:.6g}" if depth else "NA", d_alt_contig])) + "\n")

# ------------------------------------------------------------------ FOXO1 window and hotspot fragments
with open(os.path.join(tmp_dir, "foxo1_pileup.tsv"), "w") as fp:
    fp.write("library\tchrom\tpos\tref_plus\tA_p\tC_p\tG_p\tT_p\tdel_p\tother_p\tA_m\tC_m\tG_m\tT_m\tdel_m\tother_m\n")
    for w in window:
        c = tcounts(w["chrom"], int(w["pos"]))
        fp.write("\t".join([args.library, w["chrom"], w["pos"], w["ref_plus"]] + [str(int(x)) for x in c[0]] + [str(int(x)) for x in c[1]]) + "\n")
with open(os.path.join(tmp_dir, "foxo1_hotspots.tsv"), "w") as fh_:
    fh_.write("library\thotspot\tpos\tref_plus\talt_plus\tn_ref_reads\tn_alt_reads\tn_other_reads\tdepth\tn_distinct_fragments_all\tn_distinct_fragments_alt\n")
    for h in hotspots:
        p = int(h["pos"])
        c = tcounts(h["chrom"], p)
        rc, ac = BASE[h["ref_plus"]], BASE[h["alt_plus"]]
        nr, na = int(c[:, rc].sum()), int(c[:, ac].sum())
        no = int(sum(c[:, k].sum() for k in range(4) if k not in (rc, ac)))
        fh_.write("\t".join(map(str, [args.library, h["name"], p, h["ref_plus"], h["alt_plus"], nr, na, no, nr + na + no,
                                       len(hot_all[p]), len(hot_alt[p][h["alt_plus"]])])) + "\n")

# ------------------------------------------------------------------ IGH constant region
ig_count = Counter()
for qn, (genes, his, nh) in ig_frag.items():
    ig_count[(",".join(sorted(genes)), int(len(his) >= nh), "unique" if nh == 1 else "multi")] += 1
with open(os.path.join(tmp_dir, "ig_fragments.tsv"), "w") as fi:
    fi.write("library\tgene_set\tall_alignments_in_constant_exons\tmapping\tn_fragments\n")
    for (gs, comp, mp), n in sorted(ig_count.items()):
        fi.write(f"{args.library}\t{gs}\t{comp}\t{mp}\t{n}\n")

# ------------------------------------------------------------------ editing
ec = np.frombuffer(ed_cnt, dtype=np.int32).reshape(n_sites, 4)
tot = ec.sum(axis=1)
keep = np.nonzero(tot >= args.min_site_cov)[0]
with gzip.open(os.path.join(tmp_dir, "editing_sites.tsv.gz"), "wt") as fe:
    fe.write("site\tref_p\talt_p\tref_m\talt_m\n")
    for s in keep:
        fe.write(f"{s}\t{ec[s, 0]}\t{ec[s, 1]}\t{ec[s, 2]}\t{ec[s, 3]}\n")


def sense_cols(ec_, strand_arr):
    """Counts on the transcript strand of each site (needs a stranded library)."""
    if strand_call == "unstranded":
        return None
    # forward library: fragment strand = transcript strand; reverse library: opposite
    site_plus = strand_arr == 1
    want_frag_minus = site_plus if strand_call == "reverse" else ~site_plus
    ref = np.where(want_frag_minus, ec_[:, 2], ec_[:, 0])
    alt = np.where(want_frag_minus, ec_[:, 3], ec_[:, 1])
    return ref, alt


with open(os.path.join(tmp_dir, "editing_summary.tsv"), "w") as fs_:
    fs_.write("library\tsite_class\tmode\tn_sites_any\tn_sites_cov10\tref_reads\talt_reads\tpooled_index\tmean_site_fraction_cov10\n")
    for cname, cmask in (("ALU", ed_nonalu_np == 0), ("NONALU", ed_nonalu_np == 1), ("ALL", np.ones(n_sites, bool))):
        modes = {"all_fragments": (ec[:, 0] + ec[:, 2], ec[:, 1] + ec[:, 3])}
        sc = sense_cols(ec, ed_strand)
        if sc is not None:
            modes["transcript_strand"] = sc
        for mname, (rf, al) in modes.items():
            m = cmask & ((rf + al) > 0)
            m10 = cmask & ((rf + al) >= 10)
            rsum, asum = int(rf[m].sum()), int(al[m].sum())
            frac = (al[m10] / (rf[m10] + al[m10])).mean() if m10.any() else float("nan")
            idx = asum / (rsum + asum) if (rsum + asum) else float("nan")
            fs_.write(f"{args.library}\t{cname}\t{mname}\t{int(m.sum())}\t{int(m10.sum())}\t{rsum}\t{asum}\t{idx:.6g}\t{frac:.6g}\n")

# ------------------------------------------------------------------ library QC and run info
t_total = time.time() - t_start
qc = {
    "library": args.library, "cohort": args.cohort, "layout": args.layout,
    "bam_bytes": os.path.getsize(args.bam), "n_records_scanned": n_rec, "n_unique_primary": n_unique, "n_primary_mapq3": n_nh2,
    "n_reads_with_editing_site": n_edit_reads, "n_reads_with_variant_target": n_target_reads,
    "star_input_reads": log.get("Number of input reads", "NA"), "star_unique_pct": log.get("Uniquely mapped reads %", "NA"),
    "unique_records_vs_star_log": round(n_unique / (int(log["Uniquely mapped reads number"]) * (2 if args.layout == "paired" else 1)), 5)
    if int(log.get("Uniquely mapped reads number", 0)) else "NA",
    "star_avg_input_length": log.get("Average input read length", "NA"),
    "rpg_unstranded_total": col2, "rpg_sense_total": col3, "rpg_antisense_total": col4,
    "frac_sense": round(frac_sense, 4), "strand_call": strand_call,
    "n_junctions": n_j, "junction_unique_reads": tot_u, "junction_multi_reads": tot_m,
    "editing_other_base_alu": ed_other[0], "editing_other_base_nonalu": ed_other[1],
    "seconds_pass": round(t_pass, 1), "seconds_total": round(t_total, 1),
}
with open(os.path.join(tmp_dir, "library_qc.tsv"), "w") as fq:
    fq.write("\t".join(qc) + "\n" + "\t".join(str(v) for v in qc.values()) + "\n")
with open(os.path.join(tmp_dir, "run_info.json"), "w") as fr:
    json.dump({
        "args": vars(args), "python": sys.version.split()[0], "pysam": pysam.__version__, "numpy": np.__version__,
        "script_sha256": SCRIPT_SHA256, "seed": args.seed,
        "slurm_job": os.environ.get("SLURM_JOB_ID"), "slurm_array_task": os.environ.get("SLURM_ARRAY_TASK_ID"),
        "ins_ref_window": INS_REF, "ins_alt_window": INS_ALT, "resources_realpath": os.path.realpath(RES),
    }, fr, indent=1)
os.rename(tmp_dir, args.out_dir)
print(f"done {args.library}: {n_rec} records, {n_unique} unique, {t_total:.0f} s")
