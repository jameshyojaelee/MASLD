#!/usr/bin/env python3
"""EXTERNAL arbiter for the re-oriented high-PIP palindromes.
gnomAD (via Ensembl REST) is INDEPENDENT of the 1kg panel that triggered each
re-orientation. Per variant: does gnomAD side with the PANEL (study AF was
inverted -> flip was right) or with the STUDY (panel was the outlier -> a PIP=1
anchor carries a wrongly inverted beta)?
1000GENOMES:phase_3:* is reported too but is NOT independent -- it is the same
source as the local sidecar."""
import csv, json, subprocess, time, sys, collections
BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/tmp_diag/reorient_audit"
GN={"EUR":"gnomADg:nfe","AFR":"gnomADg:afr","AMR":"gnomADg:amr","EAS":"gnomADg:eas","SAS":"gnomADg:sas"}
KG={"EUR":"1000GENOMES:phase_3:EUR","AFR":"1000GENOMES:phase_3:AFR","AMR":"1000GENOMES:phase_3:AMR",
    "EAS":"1000GENOMES:phase_3:EAS","SAS":"1000GENOMES:phase_3:SAS"}
def fetch(rs):
    url=f"https://rest.ensembl.org/variation/human/{rs}?content-type=application/json;pops=1"
    for _ in range(3):
        p=subprocess.run(["wget","-q","-O","-","--tries=2","--timeout=30",url],capture_output=True,text=True)
        if p.returncode==0 and p.stdout.strip():
            try: return json.loads(p.stdout)
            except Exception: pass
        time.sleep(2)
    return None
rows=list(csv.DictReader(open(f"{BASE}/reoriented_highpip.tsv"),delimiter="\t"))
cache={}
out=[]
for r in rows:
    rs=r["rsid"]; anc=r["ancestry"]; ea=r["effect_allele"]
    if rs not in cache:
        cache[rs]=fetch(rs); time.sleep(0.35)
    d=cache[rs]
    alleles=None; gn=None; kg=None; gn_pop=GN.get(anc); kg_pop=KG.get(anc)
    if d:
        m=d.get("mappings") or [{}]
        alleles=m[0].get("allele_string")
        for p in d.get("populations",[]):
            if p["population"]==gn_pop and p["allele"]==ea: gn=p["frequency"]
            if p["population"]==kg_pop and p["allele"]==ea: kg=p["frequency"]
    study=float(r["af_on_effect"]); panel=float(r["panel_af_a1"])
    verdict="NO_GNOMAD_DATA"
    if gn is not None:
        ds=abs(gn-study); dp=abs(gn-panel)
        if dp < ds-0.10: verdict="FLIP_CORRECT (gnomAD sides with panel)"
        elif ds < dp-0.10: verdict="*** FLIP_WRONG (gnomAD sides with study) ***"
        else: verdict="ambiguous"
    out.append(dict(rsid=rs,study=r["study_name"],anc=anc,chr=r["chromosome"],pos=r["position"],
        ea=ea,oa=r["other_allele"],alleles=alleles,multiallelic=(alleles and len(alleles.split("/"))>2),
        study_af=round(study,4),panel_af=round(panel,4),
        gnomad=None if gn is None else round(gn,4),kg=None if kg is None else round(kg,4),
        margin=round(float(r["margin"]),4),pip=round(float(r["pip"]),3),
        primary=r["primary_eligible"],verdict=verdict))
w=csv.DictWriter(open(f"{BASE}/ensembl_verdicts.tsv","w"),fieldnames=list(out[0].keys()),delimiter="\t")
w.writeheader(); [w.writerow(o) for o in out]
print(f"{'rsid':13s} {'study':14s} {'an':4s} {'ea':2s} {'alleles':10s} {'studyAF':>8s} {'panelAF':>8s} {'gnomAD':>7s} {'1kg':>6s} {'pip':>5s} {'pri':>5s}  verdict")
for o in out:
    print(f"{o['rsid']:13s} {o['study'][:14]:14s} {o['anc']:4s} {o['ea']:2s} {str(o['alleles'])[:10]:10s} "
          f"{o['study_af']:8.4f} {o['panel_af']:8.4f} {str(o['gnomad']):>7s} {str(o['kg']):>6s} {o['pip']:5.2f} "
          f"{o['primary'][:5]:>5s}  {o['verdict']}")
print()
for k,v in collections.Counter(o["verdict"] for o in out).most_common(): print(f"  {v:3d}  {k}")
print(f"\n  multi-allelic in dbSNP: {sum(1 for o in out if o['multiallelic'])} of {len(out)}")
