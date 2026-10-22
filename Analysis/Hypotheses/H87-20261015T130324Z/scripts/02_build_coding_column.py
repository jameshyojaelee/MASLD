#!/usr/bin/env python3
"""H87 step 2: per-gene Gene Catalog column for fine-mapped protein-altering or
canonical-splice credible-set members (PRESPEC_H87.md, sha256 58085a41...).

Gates G1-G5, stop rules S1-S7, controls PC-1..PC-3 and NC-1..NC-3, diagnostics D1-D9,
then the gene column. No test, no p value, no family. The verdict is computed by the
prespec's Reading rule.
"""

from __future__ import annotations

import argparse
import csv
import json
import random
import re
import sys
import time
import traceback
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import h87_common as C  # noqa: E402

THRESH_COUNT = 0.1
THRESH_FINE = 0.5

PC1 = {  # symbol: (expected max PIP, expected top site, expected top study)
    "MTARC1": (0.757099466781, "1:220970028", "2022_36402844_PDFF_EUR"),
    "TM6SF2": (1.0, "19:19379549", None),
    "PNPLA3": (0.706314979086, "22:44324727", "MVP_NAFLD_AMR"),
}
S3_GENES = ["PNPLA3", "TM6SF2", "MTARC1", "GCKR"]
PC2 = {"PNPLA3": 0.706, "TM6SF2": 1.0, "MTARC1": 0.757, "GCKR": 0.968, "APOE": 1.0, "SERPINA1": 0.598}
PC3_IDS = {"ENSG00000258555": "SPECC1L-ADORA2A", "ENSG00000286231": "", "ENSG00000227184": "EPPK1"}
NC1 = ["MBOAT7", "TRIB1"]
NC2 = dict(symbol="HSD17B13", ensembl_id="ENSG00000170509", max_pip=0.0415769171263,
           state="coding_cs_member_pip_lt_0.1", n_sites_any=1, n_studies_any=2)
FIG6C = ["GLP1R", "HKDC1", "MTARC1", "THRB"]
G2_EXPECTED = dict(rows=597, sites=145, variant_ids=147, studies=29, sites_pip01=53,
                   studies_pip01=24, sites_pip05=35, studies_pip05=18)
CONTROL_SYMBOLS = set(PC1) | set(PC2) | set(NC1) | {"HSD17B13", "EPPK1"} | set(FIG6C)


class Stop(Exception):
    def __init__(self, rule, message):
        super().__init__(f"{rule}: {message}")
        self.rule = rule
        self.message = message


# ----------------------------------------------------------------------------------------
# Logging and diagnostics
# ----------------------------------------------------------------------------------------
class Diag:
    def __init__(self, out: Path):
        self.rows = []
        self.out = out
        self.logf = open(out / "build.log", "x")

    def log(self, msg=""):
        print(msg, flush=True)
        self.logf.write(str(msg) + "\n")
        self.logf.flush()

    def add(self, section, key, value):
        self.rows.append((section, key, value))
        self.log(f"[{section}] {key}: {value}")

    def write(self, name="diagnostics.tsv"):
        with open(self.out / name, "x", newline="") as handle:
            w = csv.writer(handle, delimiter="\t", lineterminator="\n")
            w.writerow(["section", "key", "value"])
            for row in self.rows:
                w.writerow([row[0], row[1], row[2] if not isinstance(row[2], (dict, list)) else json.dumps(row[2], default=str)])


# ----------------------------------------------------------------------------------------
# GTF
# ----------------------------------------------------------------------------------------
ATTR = re.compile(r'(\S+) "([^"]*)"')


def parse_attrs(text):
    tags, out = [], {}
    for k, v in ATTR.findall(text):
        if k == "tag":
            tags.append(v)
        elif k not in out:
            out[k] = v
    out["tags"] = tags
    return out


def gtf_gene_records(path: Path):
    genes = {}
    with open(path) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 8)
            if f[2] != "gene":
                continue
            a = parse_attrs(f[8])
            gid_full = a["gene_id"]
            if gid_full.endswith("_PAR_Y"):
                continue
            gid = gid_full.split(".")[0]
            genes[gid] = dict(gene_id=gid, gene_id_versioned=gid_full, gene_name=a.get("gene_name", ""),
                              gene_type=a.get("gene_type", ""), tags=a["tags"], chrom=f[0],
                              start=int(f[3]), end=int(f[4]), strand=f[6])
    return genes


def gtf_structures(path: Path, gene_ids: set[str]):
    """CDS intervals and exons per transcript for the given unversioned gene IDs."""
    cds = defaultdict(list)
    exons = defaultdict(lambda: defaultdict(list))
    with open(path) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 8)
            if f[2] not in ("CDS", "exon"):
                continue
            i = f[8].find('gene_id "')
            gid_full = f[8][i + 9:f[8].find('"', i + 9)]
            if gid_full.endswith("_PAR_Y"):
                continue
            gid = gid_full.split(".")[0]
            if gid not in gene_ids:
                continue
            s, e = int(f[3]), int(f[4])
            if f[2] == "CDS":
                cds[gid].append((f[0], s, e))
            else:
                j = f[8].find('transcript_id "')
                tid = f[8][j + 15:f[8].find('"', j + 15)]
                exons[gid][tid].append((f[0], s, e))
    return cds, exons


def splice_windows(exon_by_tx):
    """Intronic bases 1-2 next to each internal exon boundary, per transcript."""
    win = []
    for tid, ex in exon_by_tx.items():
        ex = sorted(ex, key=lambda x: x[1])
        for k, (ch, s, e) in enumerate(ex):
            if k > 0:
                win.append((ch, s - 2, s - 1))
            if k < len(ex) - 1:
                win.append((ch, e + 1, e + 2))
    return win


# ----------------------------------------------------------------------------------------
# Gene assignment rule A2-A6
# ----------------------------------------------------------------------------------------
class Rules:
    def __init__(self, v49):
        self.v49 = v49
        self.sym_pc = defaultdict(list)
        for gid, g in v49.items():
            if g["gene_type"] == "protein_coding":
                self.sym_pc[g["gene_name"]].append(gid)

    def a3(self, gid37, symbol):
        gid = gid37.split(".")[0]
        if gid in self.v49:
            return gid, "direct_v49_id"
        cands = self.sym_pc.get(symbol, []) if symbol else []
        if len(cands) == 1:
            return cands[0], "id_remapped_by_symbol"
        return None, "unmapped_v49"

    def a4(self, v49id):
        g = self.v49[v49id]
        reasons = []
        if "readthrough_gene" in g["tags"]:
            reasons.append("readthrough_gene")
        if g["gene_name"] in (v49id, g["gene_id_versioned"]) or g["gene_name"] == "":
            reasons.append("no_symbol")
        if g["gene_type"] != "protein_coding":
            reasons.append(f"gene_type={g['gene_type']}")
        return reasons

    @staticmethod
    def a2(item):
        """GRCh37 genes with a T term on a protein_coding transcript."""
        genes = {}
        for t in item.get("transcript_consequences", []) or []:
            terms = set(t.get("consequence_terms", []))
            if t.get("biotype") != "protein_coding" or not (terms & C.T):
                continue
            g = genes.setdefault(t["gene_id"], dict(gene_id37=t["gene_id"], symbol37=t.get("gene_symbol", "") or "",
                                                     terms=set(), canonical_pa=False, transcripts=[]))
            g["terms"] |= terms & C.T
            g["transcripts"].append(t.get("transcript_id"))
        # canonical PA check across all transcripts of each gene (A5)
        for t in item.get("transcript_consequences", []) or []:
            if t.get("gene_id") in genes and t.get("canonical") == 1 and set(t.get("consequence_terms", [])) & C.PA:
                genes[t["gene_id"]]["canonical_pa"] = True
        return genes

    def assign(self, item):
        res = dict(a2_genes=[], dropped=[], unmapped_v49=[], kept=[], multi_gene_site=False,
                   a5_applied=False, unmapped_reason="")
        g37 = self.a2(item)
        res["a2_genes"] = sorted(g37)
        if not g37:
            res["unmapped_reason"] = "no_protein_coding_T_term"
            return res
        v = {}
        for gid37, g in sorted(g37.items()):
            v49id, path = self.a3(gid37, g["symbol37"])
            if v49id is None:
                res["unmapped_v49"].append(f"{gid37}({g['symbol37']})")
                continue
            reasons = self.a4(v49id)
            if reasons:
                res["dropped"].append(f"{gid37}->{v49id}({self.v49[v49id]['gene_name']}):{'+'.join(reasons)}")
                continue
            rec = v.setdefault(v49id, dict(v49_id=v49id, v49_symbol=self.v49[v49id]["gene_name"], gene_ids37=[],
                                           symbols37=[], paths=set(), terms=set(), canonical_pa=False))
            rec["gene_ids37"].append(gid37)
            rec["symbols37"].append(g["symbol37"])
            rec["paths"].add(path)
            rec["terms"] |= g["terms"]
            rec["canonical_pa"] = rec["canonical_pa"] or g["canonical_pa"]
        if not v:
            if res["unmapped_v49"] and not res["dropped"]:
                res["unmapped_reason"] = "no_v49_match"
            elif res["dropped"] and not res["unmapped_v49"]:
                res["unmapped_reason"] = "only_dropped_genes"
            else:
                res["unmapped_reason"] = "only_dropped_genes_or_no_v49_match"
            return res
        kept = list(v.values())
        if len(kept) > 1:
            canon = [r for r in kept if r["canonical_pa"]]
            if len(canon) == 1:
                kept = canon
                res["a5_applied"] = True
            else:
                res["multi_gene_site"] = True
        res["kept"] = kept
        return res


# ----------------------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------------------
def fnum(x):
    return float(x) if x not in ("", None) else np.nan


def load_json(path):
    with open(path) as handle:
        return json.load(handle)


def map_pos(ens: Path, chrom, pos):
    p = ens / "map" / f"{chrom}_{pos}.json"
    if not p.exists():
        return None
    m = load_json(p).get("mappings", [])
    if len(m) != 1:
        return None
    mp = m[0]["mapped"]
    if mp["start"] != mp["end"]:
        return None
    return ("chr" + str(mp["seq_region_name"]), int(mp["start"]))


def affected_grch37(chrom, pos, ref, alt):
    pos = int(pos)
    if len(ref) == len(alt):
        return [(chrom, pos)]
    if len(ref) > len(alt):
        return [(chrom, p) for p in range(pos + 1, pos + len(ref))]
    return [(chrom, pos), (chrom, pos + 1)]


def tier_rank(scope):
    return 1 if scope.startswith("tier1") else 2


def fmt(x, nd=12):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return ""
    if isinstance(x, float):
        return f"{x:.{nd}g}"
    return str(x)


# ----------------------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------------------
def run(args, dg: Diag, state: dict):
    random.seed(args.seed)
    np.random.seed(args.seed)
    out, ens = args.out, args.ensembl
    debug = args.debug_controls_only
    dg.add("RUN", "mode", "DEBUG_CONTROLS_ONLY" if debug else ("NULL_DEBUG" if args.null_debug else "PRODUCTION"))
    dg.add("RUN", "seed", args.seed)

    # ---- D1 / G1 / S1 ------------------------------------------------------------------
    try:
        C.verify_prespec()
    except SystemExit as err:
        raise Stop("S1", str(err))
    dg.add("D1", "prespec_sha256", C.PRESPEC_SHA)
    t0 = time.time()
    chk = C.check_inputs()
    for r in chk:
        dg.add("D1", f"input:{r['role']}", f"ok={r['ok']} bytes={r['observed_bytes']} sha256={r['observed_sha256']}")
    dg.add("D1", "hash_seconds", round(time.time() - t0, 1))
    bad = [r["role"] for r in chk if not r["ok"]]
    if bad:
        raise Stop("S1", f"input sha256 or size differs: {bad}")

    # ---- G2 / G3 / S2 ------------------------------------------------------------------
    mem = pd.read_csv(C.INPUTS["members"][0], sep="\t", dtype=str, keep_default_na=False)
    mem["pip"] = mem["normalized_susie_pip"].map(fnum)
    mem["site"] = mem["chromosome"] + ":" + mem["position"]
    gated = mem[(mem.pip_sum_gate_passed == "true") & mem.consequence_category.isin(C.CODING_CATS)].copy()
    obs = dict(rows=len(gated), sites=gated.site.nunique(), variant_ids=gated.variant_id.nunique(),
               studies=gated.study.nunique(),
               sites_pip01=gated[gated.pip >= 0.1].site.nunique(), studies_pip01=gated[gated.pip >= 0.1].study.nunique(),
               sites_pip05=gated[gated.pip >= 0.5].site.nunique(), studies_pip05=gated[gated.pip >= 0.5].study.nunique())
    count_out = Path(C.INPUTS["count_output"][0]).read_text()
    line = [ln for ln in count_out.splitlines() if ln.startswith("gate=true rows_")][0]
    kv = dict(tok.rsplit("=", 1) for tok in line.split()[1:])
    from_file = dict(rows=int(kv["rows_protein_altering"]) + int(kv["rows_canonical_splice"]), sites=int(kv["sites"]),
                     variant_ids=int(kv["variant_ids"]), studies=int(kv["studies"]), sites_pip01=int(kv["sites_pip>=0.1"]),
                     studies_pip01=int(kv["studies_pip>=0.1"]), sites_pip05=int(kv["sites_pip>=0.5"]),
                     studies_pip05=int(kv["studies_pip>=0.5"]))
    dg.add("D1", "G2_observed", obs)
    dg.add("D1", "G2_count_n.out", from_file)
    repo = pd.read_csv(C.INPUTS["repo_vep"][0], sep="\t", dtype=str, keep_default_na=False).set_index("variant_id")
    g3 = int(gated.variant_id.drop_duplicates().isin(repo.index).sum())
    dg.add("D1", "G3_variant_ids_with_oriented_alleles", f"{g3} of {gated.variant_id.nunique()}")
    if obs != G2_EXPECTED or from_file != G2_EXPECTED:
        raise Stop("S2", f"G2 counts differ: observed {obs}, count_n.out {from_file}, prespec {G2_EXPECTED}")
    if g3 != 147:
        raise Stop("S2", f"G3: {g3} of 147 variant_ids have oriented alleles")

    # ---- S6 part 1: repo canonical_splice definition ------------------------------------
    src = Path(C.INPUTS["repo_splice_rule"][0]).read_text()
    m = re.search(r'def consequence_category\(.*?\n(.*?)\n\s*return "unresolved"', src, re.S)
    block = m.group(0) if m else ""
    ms = re.search(r'if tokens & \{([^}]*)\}:\s*\n\s*return "canonical_splice"', block)
    splice_terms = set(re.findall(r'"([^"]+)"', ms.group(1))) if ms else set()
    dg.add("D2", "repo_canonical_splice_rule",
           f"noncoding_genetics.py consequence_category: not protein_altering and tokens & {sorted(splice_terms)}")
    if splice_terms != C.SPL:
        raise Stop("S6", f"repo canonical_splice terms {sorted(splice_terms)} != {sorted(C.SPL)}")

    # ---- GTF pass 1: gene records; PC-3a --------------------------------------------------
    t0 = time.time()
    v49 = gtf_gene_records(C.INPUTS["gtf"][0])
    dg.add("D1", "gtf_gene_records", f"{len(v49)} genes (PAR_Y copies skipped) in {time.time() - t0:.0f}s")
    rules = Rules(v49)
    pc3a = {}
    for gid, sym in PC3_IDS.items():
        v, path = rules.a3(gid, sym)
        drop = rules.a4(v) if v else ["unmapped_v49"]
        pc3a[gid] = dict(v49=v, path=path, a4=drop)
    pc3a_ok = (pc3a["ENSG00000258555"]["v49"] == "ENSG00000258555" and "readthrough_gene" in pc3a["ENSG00000258555"]["a4"]
               and pc3a["ENSG00000286231"]["v49"] == "ENSG00000286231" and "readthrough_gene" in pc3a["ENSG00000286231"]["a4"]
               and "no_symbol" in pc3a["ENSG00000286231"]["a4"]
               and pc3a["ENSG00000227184"]["v49"] == "ENSG00000261150"
               and pc3a["ENSG00000227184"]["path"] == "id_remapped_by_symbol" and pc3a["ENSG00000227184"]["a4"] == [])
    dg.add("PC-3a", "fixed_input_result", pc3a)
    dg.add("PC-3a", "pass", pc3a_ok)
    state["pc3a"] = (pc3a_ok, pc3a)
    if not pc3a_ok:
        raise Stop("PC-3a", f"A3/A4 rule code gives {pc3a}")

    # ---- Fetch status, G4, G5, S4, S7 ----------------------------------------------------
    fs = load_json(ens / "fetch_status.json")
    dg.add("D2", "fetch_host", fs.get("host"))
    dg.add("D2", "G4_ping", fs.get("g4_ping"))
    if fs.get("stop"):
        rule = fs["stop"].split(":")[0]
        raise Stop(rule, f"fetch step: {fs['stop']}")
    man = pd.read_csv(ens / "fetch_manifest.tsv", sep="\t", dtype=str)
    bad = [r.file for r in man.itertuples() if C.sha256_file(ens / r.file) != r.sha256]
    dg.add("D2", "fetched_files", f"{len(man)} files, sha256 re-check failures {len(bad)}")
    if bad:
        raise Stop("S4", f"stored Ensembl responses changed after fetch: {bad[:5]}")
    info_sw = load_json(ens / "info_software.json") if (ens / "info_software.json").exists() else {}
    info_dt = load_json(ens / "info_data.json") if (ens / "info_data.json").exists() else {}
    dg.add("D2", "ensembl_grch37_rest_software_release", info_sw.get("release"))
    dg.add("D2", "ensembl_grch37_rest_data_releases", info_dt.get("releases"))
    dg.add("D2", "map_endpoint_servers", sorted(set(man[man.endpoint.str.startswith("/map")].endpoint)))

    def load_vep(label):
        items = {}
        for f in sorted(ens.glob(f"vep_{label}_batch*.json")):
            if f.name.endswith("_request.json"):
                continue
            for it in load_json(f):
                items[str(it["id"])] = it
        return items

    prim_ids = sorted(gated.variant_id.unique())
    vep = load_vep("primary")
    g5 = sorted(vep) == prim_ids
    dg.add("D2", "G5_identity", f"submitted {len(prim_ids)} returned {len(vep)} identical={g5}")
    if not g5:
        raise Stop("S7", f"VEP identity mismatch: missing {sorted(set(prim_ids) - set(vep))[:5]} extra {sorted(set(vep) - set(prim_ids))[:5]}")

    # ---- D2: per-variant agreement, S6 part 2 --------------------------------------------
    repo_cat = gated.drop_duplicates("variant_id").set_index("variant_id")["consequence_category"]
    agree, disagree = 0, []
    for vid in prim_ids:
        terms_any = {t for tc in vep[vid].get("transcript_consequences", []) or [] for t in tc.get("consequence_terms", [])}
        has_t = bool(terms_any & C.T)
        new_cat = "protein_altering" if terms_any & C.PA else ("canonical_splice" if terms_any & C.SPL else "none_from_T")
        if has_t:
            agree += 1
        if new_cat != repo_cat[vid]:
            disagree.append(dict(variant_id=vid, repo=repo_cat[vid], new=new_cat,
                                 new_terms=C.sort_terms(terms_any)))
    dg.add("D2", "S6_variant_ids_with_T_term_any_transcript", f"{agree} of {len(prim_ids)}")
    dg.add("D2", "category_disagreements", disagree)
    state["s6_agree"] = agree
    if agree < 140:
        raise Stop("S6", f"only {agree} of 147 variant_ids carry a T term in the new response")

    d2 = load_vep("d2indel")
    indel_sets = {}
    for form in ("4:88231392:T:TA_del", "4:88231392:T:TA_ins"):
        it = d2.get(form)
        if it is None:
            indel_sets[form] = None
            continue
        tc = [dict(gene=t.get("gene_symbol"), gene_id=t.get("gene_id"), biotype=t.get("biotype"),
                   transcript=t.get("transcript_id"), canonical=t.get("canonical"), terms=t.get("consequence_terms"))
              for t in it.get("transcript_consequences", []) or []]
        hsd_t = any((t["gene_id"] == NC2["ensembl_id"] or t["gene"] == "HSD17B13") and t["biotype"] == "protein_coding"
                    and set(t["terms"]) & C.T for t in tc)
        indel_sets[form] = dict(allele_string=it.get("allele_string"), most_severe=it.get("most_severe_consequence"),
                                transcripts=tc, hsd17b13_T_on_protein_coding=hsd_t)
    for form, s in indel_sets.items():
        dg.add("D2", f"indel_{form}", s)
    if indel_sets.get("4:88231392:T:TA_del") and indel_sets.get("4:88231392:T:TA_ins"):
        indel_ambig = (indel_sets["4:88231392:T:TA_del"]["hsd17b13_T_on_protein_coding"]
                       != indel_sets["4:88231392:T:TA_ins"]["hsd17b13_T_on_protein_coding"])
    else:
        indel_ambig = None
    dg.add("D2", "indel_orientation_ambiguous", indel_ambig)
    dg.add("D2", "repo_orientation_4:88231392:T:TA",
           f"REF {repo.loc['4:88231392:T:TA', 'reference_allele']} ALT {repo.loc['4:88231392:T:TA', 'alternate_allele']} "
           f"({repo.loc['4:88231392:T:TA', 'reference_orientation']})")

    # ---- A1-A6 on the primary set ----------------------------------------------------------
    assign = {vid: rules.assign(vep[vid]) for vid in prim_ids}
    # PC-3b: the three IDs in the real response
    pc3b = {}
    for gid in PC3_IDS:
        hits = []
        for vid in prim_ids:
            for t in vep[vid].get("transcript_consequences", []) or []:
                if t.get("gene_id") == gid:
                    hits.append(vid)
                    break
        res = []
        for vid in hits:
            a = assign[vid]
            kept = [k["v49_id"] for k in a["kept"] if gid in k["gene_ids37"]]
            dropped = [d for d in a["dropped"] if d.startswith(gid)]
            res.append(dict(variant_id=vid, in_A2=gid in a["a2_genes"], kept_as=kept, dropped=dropped,
                            unmapped_v49=[u for u in a["unmapped_v49"] if u.startswith(gid)]))
        pc3b[gid] = dict(appears=bool(hits), results=res)
    pc3b_viol = []
    for gid in ("ENSG00000258555", "ENSG00000286231"):
        for r in pc3b[gid]["results"]:
            if r["kept_as"]:
                pc3b_viol.append(gid)
    eppk = pc3b["ENSG00000227184"]
    eppk_pred = any(r["variant_id"].startswith("8:144942477") and r["kept_as"] == ["ENSG00000261150"] for r in eppk["results"])
    dg.add("PC-3b", "result", pc3b)
    dg.add("PC-3b", "readthrough_present_without_removal", pc3b_viol)
    dg.add("PC-3b", "EPPK1_8:144942477_remapped_to_ENSG00000261150", eppk_pred)
    state["pc3b"] = (pc3b, pc3b_viol, eppk_pred)
    if pc3b_viol:
        raise Stop("PC-3b", f"readthrough gene present without removal: {pc3b_viol}")

    # ---- Member -> gene map (primary, gated) ---------------------------------------------
    rows = []
    for r in gated.itertuples(index=False):
        a = assign[r.variant_id]
        base = dict(credible_set_uid=r.credible_set_uid, study=r.study, trait=r.trait, trait_scope=r.trait_scope,
                    tier=r.tier, ancestry=r.ancestry, variant_id=r.variant_id, site=r.site, chromosome=r.chromosome,
                    position=r.position, normalized_susie_pip=r.pip, repo_consequence_category=r.consequence_category,
                    a2_gene_ids37=";".join(a["a2_genes"]), a4_dropped=";".join(a["dropped"]),
                    a3_unmapped_v49=";".join(a["unmapped_v49"]), a5_canonical_rule_applied=a["a5_applied"],
                    multi_gene_site=a["multi_gene_site"])
        if not a["kept"]:
            rows.append(dict(base, ensembl_id="", symbol="", gene_ids37="", symbols37="", mapping_path="",
                             gene_T_terms="", canonical_PA=False, unmapped_reason=a["unmapped_reason"]))
        for k in a["kept"]:
            rows.append(dict(base, ensembl_id=k["v49_id"], symbol=k["v49_symbol"], gene_ids37=";".join(k["gene_ids37"]),
                             symbols37=";".join(k["symbols37"]), mapping_path=";".join(sorted(k["paths"])),
                             gene_T_terms=",".join(C.sort_terms(k["terms"])), canonical_PA=k["canonical_pa"],
                             unmapped_reason=""))
    mm = pd.DataFrame(rows)
    if args.null_debug:
        # Code exercise only: random many-to-one gene IDs and permuted PIPs, so no real outcome is computed.
        rng = np.random.default_rng(args.seed)
        cat0 = pd.read_parquet(C.INPUTS["catalog_reviewed"][0])[["ensembl_id", "symbol"]].astype(str)
        mapped = (mm.ensembl_id != "").to_numpy()
        pick = rng.integers(0, len(cat0), int(mapped.sum()))
        mm.loc[mapped, "ensembl_id"] = cat0.ensembl_id.to_numpy()[pick]
        mm.loc[mapped, "symbol"] = cat0.symbol.to_numpy()[pick]
        mm["normalized_susie_pip"] = rng.permutation(mm.normalized_susie_pip.to_numpy())
        dg.add("RUN", "null_debug_note", "NULL DEBUG: gene IDs replaced at random and PIPs permuted; every number below is meaningless")

    # ---- GTF pass 2 and D3 build check ----------------------------------------------------
    gene_ids = set(mm.ensembl_id[mm.ensembl_id != ""]) | {NC2["ensembl_id"]}
    t0 = time.time()
    cds, exons = gtf_structures(C.INPUTS["gtf"][0], gene_ids)
    dg.add("D3", "gtf_structures", f"{len(cds)} genes with CDS, {len(exons)} with exons, {time.time() - t0:.0f}s")
    site_pos38 = {}
    site_aff38 = {}
    for vid in prim_ids:
        rr = repo.loc[vid]
        ch, pos = rr["chromosome"], int(rr["position"])
        site_pos38[f"{ch}:{pos}"] = map_pos(ens, ch, pos)
        aff = [map_pos(ens, c2, p2) for c2, p2 in affected_grch37(ch, pos, rr["reference_allele"], rr["alternate_allele"])]
        site_aff38.setdefault(f"{ch}:{pos}", set()).update([a for a in aff if a is not None])

    def concordant(site, gid, terms):
        aff = site_aff38.get(site) or set()
        if not aff:
            return False, "no_grch38_mapping"
        ok_cds = ok_spl = False
        if terms & C.PA:
            ok_cds = any(c == ch and s <= p <= e for (ch, p) in aff for (c, s, e) in cds.get(gid, []))
        if terms & C.SPL:
            win = splice_windows(exons.get(gid, {}))
            ok_spl = any(c == ch and s <= p <= e for (ch, p) in aff for (c, s, e) in win)
        if ok_cds or ok_spl:
            return True, ("cds" if ok_cds else "splice_window")
        return False, "outside_v49_cds_or_splice_window"

    pairs = (mm[mm.ensembl_id != ""].groupby(["site", "ensembl_id", "symbol"])
             .agg(pair_pip=("normalized_susie_pip", "max"), terms=("gene_T_terms", lambda s: ",".join(sorted(set(",".join(s).split(","))))))
             .reset_index())
    conc = [concordant(r.site, r.ensembl_id, set(r.terms.split(","))) for r in pairs.itertuples()]
    pairs["concordant"] = [c[0] for c in conc]
    pairs["d3_reason"] = [c[1] for c in conc]
    pairs["pos_grch38"] = [f"{site_pos38[s][0]}:{site_pos38[s][1]}" if site_pos38.get(s) else "" for s in pairs.site]
    p01 = pairs[pairs.pair_pip >= THRESH_COUNT]
    p05 = pairs[pairs.pair_pip >= THRESH_FINE]
    dg.add("D3", "pairs_all", f"{int(pairs.concordant.sum())} of {len(pairs)} concordant")
    dg.add("D3", "pairs_pip_ge_0.1", f"{int(p01.concordant.sum())} of {len(p01)} concordant "
           f"({100 * p01.concordant.mean() if len(p01) else float('nan'):.1f}%; predicted >= 95%)")
    dg.add("D3", "pairs_pip_ge_0.5", f"{int(p05.concordant.sum())} of {len(p05)} concordant (predicted all)")
    dg.add("D3", "discordant_pairs", pairs[~pairs.concordant][["site", "symbol", "pair_pip", "terms", "pos_grch38", "d3_reason"]].to_dict("records"))
    lift = pd.read_csv(C.INPUTS["liftover"][0], dtype={"chromosome": str})
    lift = lift.drop_duplicates(["chromosome", "position"])
    lift["site"] = lift.chromosome + ":" + lift.position.astype(str)
    lift = lift[lift.site.isin(set(gated.site))]
    dg.add("D3", "gated_sites_in_variant_overlap_summary_without_hg38", int(lift.pos_hg38.isna().sum()))
    lift = lift.dropna(subset=["pos_hg38"])
    lcmp = []
    for r in lift.itertuples():
        e = site_pos38.get(r.site)
        lcmp.append(dict(site=r.site, repo=f"{r.chr_hg38}:{int(r.pos_hg38)}", ensembl=f"{e[0]}:{e[1]}" if e else "",
                         equal=bool(e) and e[0] == r.chr_hg38 and e[1] == int(r.pos_hg38)))
    lcmp = pd.DataFrame(lcmp)
    dg.add("D3", "repo_liftover_vs_ensembl_map",
           f"{int(lcmp.equal.sum()) if len(lcmp) else 0} of {len(lcmp)} gated sites in variant_overlap_summary.csv equal")
    if len(lcmp) and not lcmp.equal.all():
        dg.add("D3", "repo_liftover_mismatches", lcmp[~lcmp.equal].to_dict("records"))
    disc_pair = set(zip(pairs.site[~pairs.concordant], pairs.ensembl_id[~pairs.concordant]))
    mm["d3_concordant"] = [("" if not g else ((s, g) not in disc_pair)) for s, g in zip(mm.site, mm.ensembl_id)]
    mm["pos_grch38"] = [f"{site_pos38[s][0]}:{site_pos38[s][1]}" if site_pos38.get(s) else "" for s in mm.site]

    # ---- D4 mapping path counts --------------------------------------------------------------
    site_status = {}
    for vid in prim_ids:
        a = assign[vid]
        s = ":".join(vid.split(":")[:2])
        site_status.setdefault(s, []).append(a)
    n_mapped = sum(any(a["kept"] for a in v) for v in site_status.values())
    unm = defaultdict(list)
    for s, v in site_status.items():
        if not any(a["kept"] for a in v):
            unm[v[0]["unmapped_reason"]].append(s)
    remapped = sorted({f"{gid37}->{k['v49_id']}({k['v49_symbol']})" for vid in prim_ids for k in assign[vid]["kept"]
                       if "id_remapped_by_symbol" in k["paths"] for gid37 in k["gene_ids37"]})
    dropped = sorted({d for vid in prim_ids for d in assign[vid]["dropped"]})
    multi = sorted({":".join(vid.split(":")[:2]) for vid in prim_ids if assign[vid]["multi_gene_site"]})
    a5 = sorted({":".join(vid.split(":")[:2]) for vid in prim_ids if assign[vid]["a5_applied"]})
    unmapped_v49 = sorted({u for vid in prim_ids for u in assign[vid]["unmapped_v49"]})
    dg.add("D4", "sites_mapped", f"{n_mapped} of {len(site_status)}")
    dg.add("D4", "sites_unmapped_by_reason", {k: sorted(v) for k, v in unm.items()})
    dg.add("D4", "genes_remapped_by_symbol", remapped)
    dg.add("D4", "genes_dropped_by_A4", dropped)
    dg.add("D4", "genes_unmapped_v49", unmapped_v49)
    dg.add("D4", "sites_resolved_by_A5_canonical", a5)
    dg.add("D4", "multi_gene_sites", multi)

    # ---- Per-gene aggregation -------------------------------------------------------------------
    cat = pd.read_parquet(C.INPUTS["catalog_reviewed"][0])
    cat = cat[["ensembl_id", "symbol", "primary_evidence_class"]].astype(str)
    sym2id = dict(zip(cat.symbol, cat.ensembl_id))

    def top_terms(gid, vid):
        a = assign[vid]
        for k in a["kept"]:
            if k["v49_id"] == gid:
                return ",".join(C.sort_terms(k["terms"]))
        return ""

    def aggregate(mmap, pipcol="normalized_susie_pip"):
        res = {}
        for gid, g in mmap[mmap.ensembl_id != ""].groupby("ensembl_id"):
            mx = float(g[pipcol].max())
            g01 = g[g[pipcol] >= THRESH_COUNT]
            top = g[g[pipcol] == mx].copy()
            top["tr"] = top.trait_scope.map(tier_rank)
            top = top.sort_values(["tr", "study", "variant_id"]).iloc[0]
            if mx >= THRESH_FINE:
                st = "fine_mapped_coding_pip_ge_0.5"
            elif mx >= THRESH_COUNT:
                st = "coding_cs_member_pip_0.1_to_0.5"
            else:
                st = "coding_cs_member_pip_lt_0.1"
            flags = []
            if "id_remapped_by_symbol" in ";".join(g.mapping_path):
                flags.append("id_remapped_by_symbol")
            if g.multi_gene_site.any():
                flags.append("multi_gene_site")
            if (g.d3_concordant == False).any():  # noqa: E712
                flags.append("build_discordant")
            if gid == NC2["ensembl_id"] and indel_ambig and (g.variant_id == "4:88231392:T:TA").any():
                flags.append("indel_orientation_ambiguous")
            tiers = sorted(set(g01.trait_scope), key=tier_rank)
            res[gid] = dict(
                coding_cs_symbol_v49=g.symbol.iloc[0],
                coding_cs_max_pip=mx,
                coding_cs_n_studies=int(g01.study.nunique()),
                coding_cs_n_sites=int(g01.site.nunique()),
                coding_cs_n_studies_any=int(g.study.nunique()),
                coding_cs_n_sites_any=int(g.site.nunique()),
                coding_cs_top_variant=top.variant_id,
                coding_cs_top_study=top.study,
                coding_cs_top_consequence=top_terms(gid, top.variant_id) if top.variant_id in assign else "",
                coding_cs_top_trait_tier=top.trait_scope,
                coding_cs_top_pos_grch38=top.pos_grch38 if "pos_grch38" in top else "",
                coding_cs_trait_tiers=";".join(tiers),
                coding_cs_state=st,
                coding_cs_flags=";".join(flags),
            )
        return res

    if debug:
        keep_ids = {sym2id[s] for s in CONTROL_SYMBOLS if s in sym2id} | {"ENSG00000261150"}
        mm_used = mm[mm.ensembl_id.isin(keep_ids)].copy()
        dg.add("RUN", "debug_note", f"DEBUG: gene aggregation restricted to {len(keep_ids)} control genes; "
               "no Catalog-wide counts computed")
    else:
        mm_used = mm
    agg = aggregate(mm_used)

    # ---- S5 ---------------------------------------------------------------------------------
    s5_genes = sorted({r.ensembl_id for r in p05[~p05.concordant].itertuples() if agg.get(r.ensembl_id, {}).get("coding_cs_max_pip", 0) >= THRESH_FINE})
    dg.add("D3", "S5_genes_ge_0.5_with_discordant_pair_ge_0.5", [(g, agg[g]["coding_cs_symbol_v49"]) for g in s5_genes])
    state["s5"] = s5_genes

    def gene_row(sym):
        gid = sym2id.get(sym)
        return gid, agg.get(gid)

    # ---- PC-1 / S3 -------------------------------------------------------------------------------
    pc1 = {}
    s3_fail = []
    for sym, (val, site, study) in PC1.items():
        gid, r = gene_row(sym)
        if r is None:
            pc1[sym] = dict(observed=None, expected=val, ok=False)
            s3_fail.append(f"{sym} no entry")
            continue
        ok = bool(abs(r["coding_cs_max_pip"] - val) <= 1e-9)
        pc1[sym] = dict(observed=r["coding_cs_max_pip"], expected=val, ok=ok, top_variant=r["coding_cs_top_variant"],
                        top_study=r["coding_cs_top_study"], expected_site=site, expected_study=study,
                        top_site_match=r["coding_cs_top_variant"].startswith(site + ":"),
                        top_study_match=(study is None or r["coding_cs_top_study"] == study))
        if not ok:
            s3_fail.append(f"{sym} max PIP {r['coding_cs_max_pip']!r} != {val}")
    for sym in S3_GENES:
        gid, r = gene_row(sym)
        if r is None or r["coding_cs_max_pip"] < THRESH_FINE:
            s3_fail.append(f"{sym} has no entry at PIP >= 0.5")
    dg.add("PC-1", "result", pc1)
    dg.add("PC-1", "S3_failures", s3_fail)
    state["pc1"] = (pc1, s3_fail)

    # ---- PC-2 -------------------------------------------------------------------------------------
    pc2 = {}
    for sym, val in PC2.items():
        gid, r = gene_row(sym)
        if r is None:
            sites = gated[gated.variant_id.isin([v for v in prim_ids
                                                 if any(t.get("gene_symbol") == sym for t in vep[v].get("transcript_consequences", []) or [])])]
            pc2[sym] = dict(observed=None, predicted=val, state=None, pass_=False, value_matches_gate_to_3dp=False,
                            note=f"no entry; sites with this GRCh37 symbol in response: {sorted(sites.site.unique())}")
            continue
        pc2[sym] = dict(observed=round(r["coding_cs_max_pip"], 6), predicted=val, state=r["coding_cs_state"],
                        top=r["coding_cs_top_variant"], flags=r["coding_cs_flags"],
                        pass_=(r["coding_cs_state"] == "fine_mapped_coding_pip_ge_0.5"),
                        value_matches_gate_to_3dp=bool(abs(round(r["coding_cs_max_pip"], 3) - val) < 1e-9))
    dg.add("PC-2", "result", pc2)
    state["pc2"] = pc2

    # ---- NC-1 -------------------------------------------------------------------------------------
    nc1 = {}
    for sym in NC1:
        gid, r = gene_row(sym)
        nc1[sym] = "no entry (as predicted)" if r is None else dict(entry=r, members=mm[mm.ensembl_id == gid][["site", "study", "normalized_susie_pip", "gene_T_terms"]].to_dict("records"))
    dg.add("NC-1", "result", nc1)
    state["nc1"] = nc1

    # ---- NC-2 -------------------------------------------------------------------------------------
    r = agg.get(NC2["ensembl_id"])
    if r is None:
        it = vep.get("4:88231392:T:TA", {})
        nc2 = dict(entry=None, note="no entry", response_terms=[dict(gene=t.get("gene_symbol"), biotype=t.get("biotype"),
                                                                       terms=t.get("consequence_terms")) for t in it.get("transcript_consequences", []) or []],
                   assignment=assign.get("4:88231392:T:TA"))
        nc2_pass = False
    else:
        nc2 = dict(observed_max_pip=r["coding_cs_max_pip"], state=r["coding_cs_state"], n_sites_any=r["coding_cs_n_sites_any"],
                   n_studies_any=r["coding_cs_n_studies_any"], flags=r["coding_cs_flags"], top=r["coding_cs_top_variant"],
                   top_consequence=r["coding_cs_top_consequence"])
        nc2_pass = (abs(r["coding_cs_max_pip"] - NC2["max_pip"]) <= 1e-9 and r["coding_cs_state"] == NC2["state"]
                    and r["coding_cs_n_sites_any"] == NC2["n_sites_any"] and r["coding_cs_n_studies_any"] == NC2["n_studies_any"])
    nc2["pass"] = nc2_pass
    dg.add("NC-2", "result", nc2)
    state["nc2"] = nc2

    # ---- NC-3 -------------------------------------------------------------------------------------
    nc3_vep = load_vep("nc3")
    nc3_ids = sorted(mem[(mem.pip_sum_gate_passed == "true") & (mem.consequence_category == "synonymous_or_utr")
                         & (mem.pip >= 0.5)].variant_id.unique())
    nc3_assigned = {v: sorted(Rules.a2(nc3_vep[v])) for v in nc3_ids if v in nc3_vep and Rules.a2(nc3_vep[v])}
    nc3 = dict(n_variant_ids=len(nc3_ids), n_returned=sum(v in nc3_vep for v in nc3_ids), n_with_assigned_gene=len(nc3_assigned),
               assigned=nc3_assigned, predicted="<= 1", pass_=len(nc3_assigned) <= 1, s6_agreement=f"{agree} of 147")
    dg.add("NC-3", "result", nc3)
    state["nc3"] = nc3

    # ---- D6: ungated sensitivity ------------------------------------------------------------------
    ung = mem[(mem.pip_sum_gate_passed == "false") & mem.consequence_category.isin(C.CODING_CATS)].copy()
    ung["pip_d6"] = ung.susie_pip.map(fnum) / ung.raw_cs_pip_sum.map(fnum)
    d6_vep = load_vep("d6")
    d6_rows = []
    for r_ in ung.itertuples(index=False):
        it = vep.get(r_.variant_id) or d6_vep.get(r_.variant_id)
        if it is None:
            d6_rows.append(dict(variant_id=r_.variant_id, site=r_.site, study=r_.study, pip_d6=r_.pip_d6, ensembl_id="", symbol="", note="no VEP response"))
            continue
        a = rules.assign(it)
        for k in a["kept"] or [None]:
            d6_rows.append(dict(variant_id=r_.variant_id, site=r_.site, study=r_.study, pip_d6=r_.pip_d6,
                                ensembl_id=k["v49_id"] if k else "", symbol=k["v49_symbol"] if k else "",
                                note=a["unmapped_reason"]))
    d6 = pd.DataFrame(d6_rows)
    dg.add("D6", "ungated_rows_sites_variant_ids", f"{len(ung)} rows, {ung.site.nunique()} sites, {ung.variant_id.nunique()} variant_ids, "
           f"{ung.study.nunique()} studies; max susie_pip/raw_cs_pip_sum {ung.pip_d6.max():.6g}")
    dg.add("D6", "ungated_sites_not_in_gated", sorted(set(ung.site) - set(gated.site)))
    if not debug:
        prim_genes = set(agg)
        d6_genes = set(d6.ensembl_id[d6.ensembl_id != ""])
        only = sorted(d6_genes - prim_genes)
        dg.add("D6", "genes_only_without_PIP_sum_gate",
               [dict(ensembl_id=g, symbol=d6[d6.ensembl_id == g].symbol.iloc[0], max_pip_d6=float(d6[d6.ensembl_id == g].pip_d6.max()),
                     sites=sorted(d6[d6.ensembl_id == g].site.unique())) for g in only])
        dg.add("D6", "pip_ge_0.1_and_0.5_site_sets_identical_with_and_without_gate", bool((ung.pip_d6 < 0.1).all()))
        dg.add("D6", "unmapped_ungated_variant_ids", sorted(set(d6.variant_id[d6.ensembl_id == ""])))

    # ---- D5 / summary --------------------------------------------------------------------------------
    gcol = None
    summary = {}
    if not debug:
        gdf = pd.DataFrame.from_dict(agg, orient="index")
        gdf.index.name = "ensembl_id"
        gdf = gdf.reset_index()
        gdf["in_catalog"] = gdf.ensembl_id.isin(set(cat.ensembl_id))
        gdf = gdf.merge(cat[["ensembl_id", "primary_evidence_class"]], on="ensembl_id", how="left")
        classes = ["unresolved", "genetically_anchored", "established_state_associated", "untested", "concordant"]
        for lab, thr in (("any", -1.0), ("ge_0.1", THRESH_COUNT), ("ge_0.5", THRESH_FINE)):
            s = gdf[gdf.coding_cs_max_pip >= thr]
            sc = s[s.in_catalog]
            counts = {c: int((sc.primary_evidence_class == c).sum()) for c in classes}
            summary[lab] = dict(genes=len(s), in_catalog=len(sc), class_counts=counts,
                                not_in_catalog=sorted(s[~s.in_catalog].coding_cs_symbol_v49 + "(" + s[~s.in_catalog].ensembl_id + ")"))
            dg.add("D5", f"genes_{lab}", summary[lab])
        s05 = gdf[(gdf.coding_cs_max_pip >= THRESH_FINE) & gdf.in_catalog]
        s05_clean = s05[~s05.ensembl_id.isin(s5_genes)]
        summary["ge_0.5_without_S5_flag"] = dict(in_catalog=len(s05_clean),
                                                 class_counts={c: int((s05_clean.primary_evidence_class == c).sum()) for c in classes})
        dg.add("D5", "genes_ge_0.5_in_catalog_without_S5_build_flag", summary["ge_0.5_without_S5_flag"])
        for lab, thr in (("ge_0.5", THRESH_FINE), ("ge_0.1", THRESH_COUNT)):
            s = gdf[(gdf.coding_cs_max_pip >= thr) & gdf.in_catalog]
            summary[f"unresolved_or_untested_gain_{lab}"] = dict(
                unresolved=int((s.primary_evidence_class == "unresolved").sum()),
                untested=int((s.primary_evidence_class == "untested").sum()))
            dg.add("D5", f"unresolved_or_untested_genes_gaining_entry_{lab}", summary[f"unresolved_or_untested_gain_{lab}"])
        ct = pd.crosstab(gdf[gdf.in_catalog].primary_evidence_class, gdf[gdf.in_catalog].coding_cs_state)
        dg.add("D5", "crosstab_class_by_state", ct.to_dict())
        state["summary"] = summary

        # D7: production against gate
        gate = pd.read_csv(C.INPUTS["gate_gene_summary"][0], sep="\t", dtype={"ensembl_id": str})
        gate_ids = set(gate.ensembl_id)
        gate_num = set(gate.ensembl_id[gate.max_pip.notna()])
        # gate site-gene pairs: liftover member map plus gate VEP JSON (merge_vep.py logic)
        gmap = pd.read_csv(C.INPUTS["gate_member_map"][0], sep="\t", dtype=str, keep_default_na=False)
        gpairs = defaultdict(set)
        for r_ in gmap.itertuples():
            if r_.ensembl_id != "NA":
                gpairs[r_.ensembl_id].add(r_.site)
        gate_pa = {"missense_variant", "stop_gained", "stop_lost", "start_lost", "frameshift_variant", "inframe_insertion",
                   "inframe_deletion", "protein_altering_variant", "splice_donor_variant", "splice_acceptor_variant"}
        for f in sorted((C.ROUND4 / "gate/R15-01/tmp/vep").glob("*.json")):
            ch, ps = f.name[:-5].split("_")
            for it in load_json(f):
                for t in it.get("transcript_consequences", []) or []:
                    if t.get("biotype") == "protein_coding" and gate_pa & set(t.get("consequence_terms", [])):
                        gpairs[t["gene_id"]].add(f"{ch}:{ps}")
        prod_ids = set(gdf.ensembl_id)
        lost, gained = [], []
        site_assign = defaultdict(list)
        for vid in prim_ids:
            site_assign[":".join(vid.split(":")[:2])].append(assign[vid])
        ung_sites = set(ung.site) - set(gated.site)
        for gid in sorted(gate_ids - prod_ids):
            sym = gate[gate.ensembl_id == gid].symbol.iloc[0]
            sites = sorted(gpairs.get(gid, set()))
            reasons = []
            if gid not in gate_num or (sites and set(sites) <= ung_sites):
                reasons.append("only ungated members (PIP-sum gate)")
            for s in sites:
                for a in site_assign.get(s, []):
                    if any(d.startswith(gid) or f"->{gid}(" in d for d in a["dropped"]):
                        reasons.append(f"{s}: dropped by A4 ({[d for d in a['dropped'] if gid in d][0].split(':')[-1]})")
                    elif any(gid in k["gene_ids37"] for k in a["kept"]):
                        reasons.append(f"{s}: GRCh37 ID remapped by A3 to {[k['v49_id'] for k in a['kept'] if gid in k['gene_ids37']]}")
                    elif gid in a["a2_genes"]:
                        reasons.append(f"{s}: removed by A5 canonical-transcript rule or no v49 match")
                    else:
                        g37 = [t for vid in prim_ids if vid.startswith(s + ":") for t in vep[vid].get("transcript_consequences", []) or []
                               if t.get("gene_symbol") == sym or t.get("gene_id") == gid]
                        if g37:
                            reasons.append(f"{s}: GRCh37 VEP returns {sym} without a T term on a protein_coding transcript "
                                           f"({sorted({x for t in g37 for x in t.get('consequence_terms', [])})})")
                        else:
                            reasons.append(f"{s}: {sym} not in the GRCh37 VEP response (gate used the GENCODE v49 liftover CDS/exon+/-2 path)")
            lost.append(dict(ensembl_id=gid, symbol=sym, gate_max_pip=gate[gate.ensembl_id == gid].max_pip.iloc[0],
                             gate_sites=sites, reasons=sorted(set(reasons)) or ["no gate site found"]))
        gate_sites_all = set().union(*gpairs.values()) if gpairs else set()
        for gid in sorted(prod_ids - gate_ids):
            g = mm[mm.ensembl_id == gid]
            sites = sorted(g.site.unique())
            reasons = []
            for s in sites:
                if s not in gate_sites_all:
                    reasons.append(f"{s}: site not mapped to any gene by the gate (liftover found no v49 CDS/exon+/-2 hit at the anchor base, or VEP path not run)")
                else:
                    other = sorted(k for k, v in gpairs.items() if s in v)
                    reasons.append(f"{s}: gate assigned {other} at this site")
            if "id_remapped_by_symbol" in ";".join(g.mapping_path):
                reasons.append("A3 symbol remap from a GRCh37 ID")
            gained.append(dict(ensembl_id=gid, symbol=g.symbol.iloc[0], max_pip=float(g.normalized_susie_pip.max()),
                               sites=sites, reasons=sorted(set(reasons))))
        dg.add("D7", "gate_gene_rows", f"{len(gate)} (numeric max PIP {len(gate_num)})")
        dg.add("D7", "production_genes_any_pip", len(prod_ids))
        dg.add("D7", "lost_against_gate", lost)
        dg.add("D7", "gained_against_gate", gained)
        state["d7"] = (lost, gained)
        common = sorted(gate_ids & prod_ids)
        moved = []
        for gid in common:
            gv = gate[gate.ensembl_id == gid].max_pip.iloc[0]
            pv = agg[gid]["coding_cs_max_pip"]
            if not (pd.isna(gv) and pd.isna(pv)) and (pd.isna(gv) or abs(gv - pv) > 1e-6):
                moved.append(dict(ensembl_id=gid, symbol=agg[gid]["coding_cs_symbol_v49"], gate=gv, production=pv))
        dg.add("D7", "max_pip_changed_against_gate", moved)
        gcol = gdf

    # ---- D8: Figure 6C rows ------------------------------------------------------------------------
    d8 = {}
    for sym in FIG6C:
        gid, r = gene_row(sym)
        cls = cat[cat.symbol == sym].primary_evidence_class.tolist()
        d8[sym] = dict(ensembl_id=gid, primary_evidence_class=cls[0] if cls else None,
                       **(r if r else dict(coding_cs_state="no_coding_cs_member")))
    dg.add("D8", "fig6c_rows", d8)
    state["d8"] = d8

    # ---- D9: per gene at PIP >= 0.5, study list --------------------------------------------------
    d9 = []
    for gid, r in sorted(agg.items(), key=lambda kv: -kv[1]["coding_cs_max_pip"]):
        if r["coding_cs_max_pip"] < THRESH_FINE:
            continue
        g = mm_used[(mm_used.ensembl_id == gid) & (mm_used.normalized_susie_pip >= THRESH_COUNT)]
        studies = (g.groupby(["study", "trait_scope", "ancestry"]).normalized_susie_pip.max().reset_index()
                   .sort_values(["trait_scope", "study"]))
        d9.append(dict(ensembl_id=gid, symbol=r["coding_cs_symbol_v49"], max_pip=r["coding_cs_max_pip"],
                       studies=[f"{s.study}|{s.trait_scope}|{s.ancestry}|{s.normalized_susie_pip:.4g}" for s in studies.itertuples()]))
    dg.add("D9", "studies_per_gene_ge_0.5", d9)

    # ---- Stop rules evaluated before any column is written -------------------------------------
    if s3_fail and not args.null_debug:
        raise Stop("S3", "; ".join(s3_fail))
    if len(s5_genes) > 2 and not args.null_debug:
        raise Stop("S5", f"{len(s5_genes)} genes at max PIP >= 0.5 are build-discordant: {s5_genes}")

    # ---- Write outputs ----------------------------------------------------------------------------------
    prefix = "DEBUG_" if debug else ("NULL_" if args.null_debug else "")
    mm_out = mm.copy()
    mm_out["normalized_susie_pip"] = mm_out.normalized_susie_pip.map(fmt)
    mm_out.sort_values(["ensembl_id", "site", "study", "variant_id"]).to_csv(
        out / f"{prefix}h87_coding_member_gene_map.tsv", sep="\t", index=False)
    pairs.to_csv(out / f"{prefix}h87_d3_site_gene_pairs.tsv", sep="\t", index=False)
    d6.to_csv(out / f"{prefix}h87_d6_ungated_member_map.tsv", sep="\t", index=False)
    if debug:
        pd.DataFrame.from_dict(agg, orient="index").to_csv(out / "DEBUG_control_gene_rows.tsv", sep="\t")
        return
    cols = ["coding_cs_max_pip", "coding_cs_n_studies", "coding_cs_n_sites", "coding_cs_n_studies_any", "coding_cs_n_sites_any",
            "coding_cs_top_variant", "coding_cs_top_study", "coding_cs_top_consequence", "coding_cs_top_trait_tier",
            "coding_cs_top_pos_grch38", "coding_cs_trait_tiers", "coding_cs_state", "coding_cs_flags"]
    full = cat[["ensembl_id", "symbol", "primary_evidence_class"]].copy()
    add = pd.DataFrame.from_dict(agg, orient="index")[cols]
    add.index.name = "ensembl_id"
    full = full.merge(add.reset_index(), on="ensembl_id", how="left")
    full["coding_cs_state"] = full.coding_cs_state.fillna("no_coding_cs_member")
    for c in ("coding_cs_n_studies", "coding_cs_n_sites", "coding_cs_n_studies_any", "coding_cs_n_sites_any"):
        full[c] = full[c].astype("Int64")
    full["coding_cs_max_pip"] = full.coding_cs_max_pip.map(fmt)
    full = full.sort_values("ensembl_id")
    assert len(full) == 28233 and full.ensembl_id.is_unique
    full.to_csv(out / f"{prefix}h87_coding_cs_gene_column.tsv", sep="\t", index=False)
    gcol[~gcol.in_catalog].to_csv(out / f"{prefix}h87_coding_genes_not_in_catalog.tsv", sep="\t", index=False)
    dg.add("OUT", "gene_column_rows", len(full))
    dg.add("OUT", "gene_column_state_counts", full.coding_cs_state.value_counts().to_dict())
    state["written"] = True


# ----------------------------------------------------------------------------------------
# RESULT
# ----------------------------------------------------------------------------------------
def verdict(state):
    if state.get("stop"):
        rule = state["stop"][0]
        if rule == "S4":
            return "UNTESTABLE_OR_STOPPED_S4"
        return f"STOPPED_{rule}"
    apoe_serp_miss = [s for s in ("APOE", "SERPINA1") if not state["pc2"][s]["pass_"]]
    n5 = len(state.get("s5", []))
    if n5 in (1, 2) or apoe_serp_miss:
        return "PARTIAL"
    return "PASS"


def write_result(out: Path, state: dict, dg: Diag, args):
    v = verdict(state)
    L = []
    L.append("# RESULT H87: per-gene coding fine-mapping column beside the eQTL-colocalization class")
    L.append("")
    L.append(f"Run directory: `{out}`. Prespec `PRESPEC_H87.md` sha256 `{C.PRESPEC_SHA}`. Seed {args.seed}. "
             f"Written {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}.")
    L.append("")
    L.append(f"**Verdict (computed by the prespec Reading rule): {v}**")
    L.append("")
    L.append("H87 is a Resource annotation. It computes no test and no p value, so there is no testing family and no "
             "multiple-testing correction. Unit: the gene as a Catalog row; studies are counted, not replicates.")
    L.append("")
    if not state.get("stop") and "summary" not in state:
        L.append("DEBUG run restricted to control genes: no Catalog-wide numbers were computed.")
    elif state.get("stop"):
        L.append(f"Stop rule fired: {state['stop'][0]}: {state['stop'][1]}")
        L.append("")
        L.append("No column and no new number. The header change \"Genetics\" to \"eQTL colocalization\" can still go ahead "
                 "because it rests on the confirmed review finding, not on H87.")
    else:
        s = state["summary"]
        L.append("## Key numbers")
        L.append("")
        L.append("| Set | Genes with an entry | In Catalog | unresolved | genetically_anchored | established_state_associated | untested | concordant |")
        L.append("|---|---|---|---|---|---|---|---|")
        for lab, name in (("any", "any PIP"), ("ge_0.1", "PIP >= 0.1"), ("ge_0.5", "PIP >= 0.5")):
            c = s[lab]["class_counts"]
            L.append(f"| {name} | {s[lab]['genes']} | {s[lab]['in_catalog']} | {c['unresolved']} | {c['genetically_anchored']} | "
                     f"{c['established_state_associated']} | {c['untested']} | {c['concordant']} |")
        c = s["ge_0.5_without_S5_flag"]
        L.append(f"| PIP >= 0.5 without an S5 build flag | | {c['in_catalog']} | {c['class_counts']['unresolved']} | "
                 f"{c['class_counts']['genetically_anchored']} | {c['class_counts']['established_state_associated']} | "
                 f"{c['class_counts']['untested']} | {c['class_counts']['concordant']} |")
        L.append("")
        L.append(f"Genes not in the Catalog: any PIP {s['any']['not_in_catalog']}; PIP >= 0.5 {s['ge_0.5']['not_in_catalog']}.")
        L.append("")
        for lab in ("ge_0.5", "ge_0.1"):
            g = s[f"unresolved_or_untested_gain_{lab}"]
            L.append(f"- Catalog genes classed `unresolved` or `untested` that gain a visible coding entry at PIP {lab[3:]} or above: "
                     f"{g['unresolved']} unresolved, {g['untested']} untested.")
        L.append(f"- S5 build-discordant genes at PIP >= 0.5: {len(state['s5'])} {state['s5']}.")
        L.append(f"- S6 agreement: {state['s6_agree']} of 147 variant_ids carry a term from T in the new GRCh37 response (stop below 140).")
        L.append("")
    L.append("## Controls (predicted before the run, observed here)")
    L.append("")
    if "pc3a" in state:
        L.append(f"- PC-3a fixed-input A3/A4 check: {'pass' if state['pc3a'][0] else 'FAIL'} ({json.dumps(state['pc3a'][1])}).")
    if "pc3b" in state:
        pc3b, viol, eppk = state["pc3b"]
        L.append(f"- PC-3b real response: SPECC1L-ADORA2A appears={pc3b['ENSG00000258555']['appears']}, "
                 f"ENSG00000286231 appears={pc3b['ENSG00000286231']['appears']}, present without removal={viol}; "
                 f"EPPK1 ENSG00000227184 appears={pc3b['ENSG00000227184']['appears']}, remapped to ENSG00000261150 at 8:144942477={eppk}.")
    if "pc1" in state:
        pc1, s3 = state["pc1"]
        for sym, r in pc1.items():
            L.append(f"- PC-1 {sym}: observed max PIP {r.get('observed')!r}, predicted {r['expected']!r}, within 1e-9 = {r['ok']}; "
                     f"top {r.get('top_variant')} in {r.get('top_study')}.")
        L.append(f"- S3 failures: {s3 or 'none'}.")
    if "pc2" in state:
        for sym, r in state["pc2"].items():
            L.append(f"- PC-2 {sym}: observed {r['observed']}, predicted {r['predicted']}, state {r['state']}, at >= 0.5 = {r['pass_']}, "
                     f"equals gate value to 3 dp = {r['value_matches_gate_to_3dp']}"
                     + (f"; {r.get('note')}" if r.get("note") else "") + ".")
    if "nc1" in state:
        for sym, r in state["nc1"].items():
            L.append(f"- NC-1 {sym}: {r if isinstance(r, str) else json.dumps(r, default=str)}.")
    if "nc2" in state:
        L.append(f"- NC-2 HSD17B13: {json.dumps(state['nc2'], default=str)} (predicted max PIP {NC2['max_pip']}, state {NC2['state']}, "
                 f"1 site, 2 studies).")
    if "nc3" in state:
        n = state["nc3"]
        L.append(f"- NC-3: {n['n_with_assigned_gene']} of {n['n_variant_ids']} synonymous_or_utr variant_ids at PIP >= 0.5 get an "
                 f"assigned gene (predicted <= 1; pass={n['pass_']}); {n['assigned']}.")
    L.append("")
    if "d8" in state:
        L.append("## Figure 6C rows (D8)")
        L.append("")
        for sym, r in state["d8"].items():
            if r.get("coding_cs_max_pip") is not None:
                L.append(f"- {sym} ({r['primary_evidence_class']}): max PIP {r['coding_cs_max_pip']:.4g}, {r['coding_cs_top_consequence']}, "
                         f"top {r['coding_cs_top_variant']} in {r['coding_cs_top_study']}, studies at PIP >= 0.1 {r['coding_cs_n_studies']}, "
                         f"state {r['coding_cs_state']}, flags '{r['coding_cs_flags']}'.")
            else:
                L.append(f"- {sym} ({r['primary_evidence_class']}): none (no_coding_cs_member).")
        L.append("")
    if "d7" in state:
        lost, gained = state["d7"]
        L.append("## Production against the gate (D7)")
        L.append("")
        L.append(f"Lost ({len(lost)}):")
        for r in lost:
            L.append(f"- {r['symbol']} ({r['ensembl_id']}): {'; '.join(r['reasons'])}")
        L.append(f"Gained ({len(gained)}):")
        for r in gained:
            L.append(f"- {r['symbol']} ({r['ensembl_id']}, max PIP {r['max_pip']:.4g}): {'; '.join(r['reasons'])}")
        L.append("")
    if not state.get("stop") and "summary" in state:
        s = state["summary"]
        c = s["ge_0.5_without_S5_flag"]
        cc = c["class_counts"]
        mt = state["d8"].get("MTARC1", {})
        mt_pip = f"{mt.get('coding_cs_max_pip'):.2f}" if mt.get("coding_cs_max_pip") is not None else "none"
        L.append("## Proposed text changes (returned as data; not installed)")
        L.append("")
        L.append(f"- RESULTS / Figure 6 legend sentence: \"{c['in_catalog']} Catalog genes carry a fine-mapped protein-altering or "
                 f"canonical-splice credible-set member at normalized SuSiE PIP >= 0.5 ({cc['unresolved']} unresolved, "
                 f"{cc['genetically_anchored']} genetically anchored, {cc['established_state_associated']} established-state associated, "
                 f"{cc['untested']} untested, {cc['concordant']} concordant). This field reports fine-mapped coding variants; it does not "
                 f"change the eQTL-colocalization class and does not identify a causal gene.\"")
        L.append(f"- Figure 6C: header \"Genetics\" to \"eQTL colocalization\"; new column \"Coding fine-mapping\"; MTARC1 prints "
                 f"\"PIP {mt_pip}, {(mt.get('coding_cs_top_consequence') or '').replace('_variant', '')}\".")
        L.append(f"- docs/manuscript/05_figure_legends.md:437-441: MTARC1 \"with a fine-mapped missense credible-set member at "
                 f"normalized PIP {mt_pip}\".")
        L.append(f"- docs/RESULTS.md:2746-2748 addition: \"({s['ge_0.5']['in_catalog']} genes at normalized PIP >= 0.5)\" per the prespec text.")
        L.append("- These edits wait for the corrected COLOC rerun to be adopted (prespec, Dependence on releases).")
        L.append("")
    L.append("## Files")
    L.append("")
    for f in sorted(out.iterdir()):
        if f.is_file():
            L.append(f"- `{f.name}`")
    L.append("")
    L.append("All diagnostics D1-D9 are in `diagnostics.tsv` and `build.log`.")
    text = "\n".join(L) + "\n"
    name = ("DEBUG_RESULT_H87.md" if args.debug_controls_only else
            "NULL_RESULT_H87.md" if args.null_debug else "RESULT_H87.md")
    with open(out / name, "x") as handle:
        handle.write(text)
    dg.log(f"VERDICT {v}")
    return v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ensembl", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--seed", type=int, default=C.SEED)
    ap.add_argument("--debug-controls-only", action="store_true")
    ap.add_argument("--null-debug", action="store_true",
                    help="exercise the production code on random gene IDs and permuted PIPs (no real outcome)")
    args = ap.parse_args()
    if args.seed != C.SEED:
        raise SystemExit(f"seed must be {C.SEED}")
    dg = Diag(args.out)
    state = {}
    try:
        run(args, dg, state)
    except Stop as err:
        state["stop"] = (err.rule, err.message)
        dg.add("STOP", err.rule, err.message)
    except Exception:  # an unexpected error is a failure, not a stop rule
        dg.log(traceback.format_exc())
        dg.write("diagnostics_partial.tsv")
        raise
    dg.write()
    write_result(args.out, state, dg, args)


if __name__ == "__main__":
    main()
