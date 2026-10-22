#!/bin/bash
# H87 n-by-group counts on the gated coding credible-set members (no gene mapping, no outcome).
F=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/credible_set_members_annotated.tsv
V=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/candidates/noncoding-dna-vep-2026-08-12/run_19782908/credible_set_variant_consequences_complete.tsv
awk -F'\t' 'NR>1{rows[$17]++} END{for(k in rows) print "all member rows gate="k, rows[k]}' $F
for g in true false any; do
awk -F'\t' -v G=$g 'NR>1 && (G=="any"||$17==G) && ($20=="protein_altering"||$20=="canonical_splice"){
 r[$20]++; s[$11":"$12]=1; v[$10]=1; st[$3]=1; tsite[$5" "$11":"$12]=1;
 if($18>=0.1){a[$11":"$12]=1; as[$3]=1}
 if($18>=0.5){h[$11":"$12]=1; hs[$3]=1; hv[$10]=1; ht[$5" "$11":"$12]=1; hts[$5" "$3]=1}}
 END{printf "gate=%s rows_protein_altering=%d rows_canonical_splice=%d sites=%d variant_ids=%d studies=%d sites_pip>=0.1=%d studies_pip>=0.1=%d sites_pip>=0.5=%d variant_ids_pip>=0.5=%d studies_pip>=0.5=%d\n",G,r["protein_altering"],r["canonical_splice"],length(s),length(v),length(st),length(a),length(as),length(h),length(hv),length(hs);
 for(k in ht){split(k,x," ");n[x[1]]++} for(k in hts){split(k,x," ");m[x[1]]++} for(k in n) printf "  gate=%s pip>=0.5 %s sites=%d studies=%d\n",G,k,n[k],m[k]}' $F
done
# coding variant_ids present in the repo GRCh37 VEP table (oriented alleles)
awk -F'\t' 'NR==FNR{if(FNR>1 && $17=="true" && ($20=="protein_altering"||$20=="canonical_splice")) c[$10]=1; next} FNR>1 && ($1 in c){f++} END{print "gated coding variant_ids",length(c),"found in repo VEP table",f}' $F $V
