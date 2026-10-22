import json, collections, sys
p="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/catalog-r8-completion-20261014T141645Z_22236527/molecular_objects/molecular_object_sidecar.jsonl"
recs=[json.loads(l) for l in open(p)]
print("n_records",len(recs))
keys=set()
for r in recs: keys|=set(r)
print("keys",sorted(keys))
def cnt(f): return collections.Counter(str(f(r)) for r in recs)
print("candidate_object",cnt(lambda r:r["candidate_object"]))
print("noncoding_dna_context",cnt(lambda r:r["noncoding_dna_context"]))
print("cs_link",cnt(lambda r:r["credible_set_to_gene_link_status"]))
print("cs_id",cnt(lambda r:r["credible_set_id"]=="" ))
print("coding_pip",cnt(lambda r:r["credible_set_coding_pip_mass"]))
print("noncoding_pip",cnt(lambda r:r["credible_set_noncoding_pip_mass"]))
print("target_link_basis",cnt(lambda r:r["target_link_basis"]))
aa=[json.loads(r["assay_applicability"]) for r in recs]
for k in ["bulk_rna","chromatin","proteomics","spatial"]:
    print("applic",k,collections.Counter(a.get(k) for a in aa))
print("rna_call",cnt(lambda r:r["rna_call"]))
print("genetics_call",cnt(lambda r:r["genetics_call"]))
print("rna_call x bulk_rna", collections.Counter((r["rna_call"],json.loads(r["assay_applicability"])["bulk_rna"]) for r in recs))
print("biotype",cnt(lambda r:r["gene_biotype"]).most_common(8))
print("lncrna class",cnt(lambda r:r["lncrna_genomic_class"]))
for g in ["GNMT","HKDC1","THRB","MTARC1","GLP1R","PNPLA3"]:
    for r in recs:
        if r["symbol"]==g:
            print(g, {k:r[k] for k in ["candidate_object","primary_evidence_class","genetics_call","rna_call","assay_applicability","noncoding_dna_context","credible_set_id","credible_set_to_gene_link_status","credible_set_coding_pip_mass","next_experiment_rule_id"]})
# any field mentioning regulatory/element
reg=[r for r in recs if any(isinstance(v,str) and ("regulatory" in v.lower() or "element" in v.lower() or "enhancer" in v.lower()) for v in r.values())]
print("records with regulatory/element/enhancer in any string value",len(reg))
vals=collections.Counter()
for r in reg:
    for k,v in r.items():
        if isinstance(v,str) and ("regulatory" in v.lower() or "element" in v.lower() or "enhancer" in v.lower()): vals[(k,v[:120])]+=1
print(vals.most_common(10))
