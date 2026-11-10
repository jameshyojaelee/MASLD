#!/usr/bin/env python3
from pathlib import Path
import json
import numpy as np
import pandas as pd
import pysam

RUN=Path(__file__).resolve().parents[1];ROOT=RUN.parents[2];OUT=RUN/'gram_inputs'
ref=ROOT/'data/external/allelic_refs/kgp_highcov'
ped=pd.read_csv(ref/'20130606_g1k_3202_samples_ped_population.txt',sep=r'\s+',dtype=str)
samples=ped[(ped.Superpopulation=='EUR') & (ped.FatherID=='0') & (ped.MotherID=='0')].sort_values('SampleID').drop_duplicates('FamilyID')
assert len(samples)==388
samples.to_csv(OUT/'reference_participants.tsv',sep='\t',index=False)
vcf=pysam.VariantFile(ref/'1kGP_high_coverage_Illumina.chr17.filtered.SNV_INDEL_SV_phased_panel.vcf.gz')
assert set(samples.SampleID).issubset(vcf.header.samples)
vcf.subset_samples(samples.SampleID.tolist())
d=pd.read_csv(OUT/'mapped_reference_events.tsv',sep='\t')
target_positions=set(d.loc[d.reciprocal_unique,'position_hg38'].astype(int))
candidates={p:[] for p in target_positions}
for record in vcf.fetch('chr17',min(target_positions)-1,max(target_positions)):
    if record.pos in candidates and len(record.alts or [])==1:
        candidates[record.pos].append(record.copy())
columns=[];selected=[];coverage=[]
for r in d.itertuples(index=False):
    row=dict(SNP=r.SNP,position=r.position,position_hg38=r.position_hg38)
    hits=[]
    if r.reciprocal_unique:
        for rec in candidates[int(r.position_hg38)]:
            direct=(rec.ref,rec.alts[0])==(r.allele1,r.allele2)
            reverse=(len(r.allele1)==len(r.allele2)==1 and
                     (rec.ref,rec.alts[0])==(r.allele2,r.allele1))
            if direct or reverse: hits.append((rec,0 if direct else 1))
    if len(hits)!=1:
        row['status']='nonunique_liftover' if not r.reciprocal_unique else 'absent_event' if not hits else 'multiple_VCF_records'
        coverage.append(row);continue
    rec,index=hits[0]
    dosage=[]
    for name in samples.SampleID:
        gt=rec.samples[name]['GT']
        dosage.append(sum(a==index for a in gt) if len(gt)==2 and all(a in (0,1) for a in gt) else np.nan)
    x=np.asarray(dosage,dtype=float)
    row.update(missing_fraction=float(np.isnan(x).mean()),effect_allele1_frequency=float(np.nanmean(x)/2),vcf_ref=rec.ref,vcf_alt=rec.alts[0])
    if row['missing_fraction']>.05:
        row['status']='excess_missing_genotypes';coverage.append(row);continue
    x[np.isnan(x)]=np.nanmean(x)
    if np.std(x)==0:
        row['status']='monomorphic';coverage.append(row);continue
    row['status']='retained';coverage.append(row)
    columns.append((x-x.mean())/np.linalg.norm(x-x.mean()));selected.append(r.SNP)
X=np.column_stack(columns);R=X.T@X
assert np.isfinite(R).all() and np.allclose(np.diag(R),1,atol=1e-12)
R.ravel(order='F').tofile(OUT/'genotype_LD_float64.bin')
pd.DataFrame({'SNP':selected}).to_csv(OUT/'genotype_LD_order.tsv',sep='\t',index=False)
pd.DataFrame(coverage).to_csv(OUT/'variant_coverage.tsv',sep='\t',index=False)
np.savez_compressed(OUT/'standardized_dosages.npz',X=X)
summary=dict(reference_participants=len(samples),input_events=len(d),retained_events=len(selected),
             status_counts=pd.DataFrame(coverage).status.value_counts().to_dict(),source_matched_to_MVP=False)
(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary,flush=True)
