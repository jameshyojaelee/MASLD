import pandas as pd
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
mem=pd.read_csv(B+'GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/credible_set_members_annotated.tsv',sep='\t',low_memory=False)
print([c for c in mem.columns][:60])
# known variants: PNPLA3 rs738409 hg19 22:44324727 / hg38 22:43928847; TM6SF2 rs58542926 hg19 19:19379549 / hg38 19:19268740
for c,p19,p38,lab in [(22,44324727,43928847,'PNPLA3 rs738409'),(19,19379549,19268740,'TM6SF2 rs58542926')]:
    for p,b in [(p19,'hg19'),(p38,'hg38')]:
        h=mem[(mem.chromosome.astype(str)==str(c))&(mem.position==p)] if True else None
        if h is not None:
            print(lab,b,p,h.groupby('ancestry').size().to_dict(), h[['study','consequence_category','lineage_accessibility_status']].drop_duplicates().head(6).to_string())
# liftover failure by study for non-EUR
nc=mem[mem.pip_sum_gate_passed.astype(str).str.lower()=='true']
t=nc.groupby('study').lineage_accessibility_status.apply(lambda s:(s=='unresolved_liftover').mean()).sort_values()
print(t.round(3).to_string())
print(nc[nc.lineage_accessibility_status=='unresolved_liftover'].variant_id.head(8).tolist())
print(nc[(nc.lineage_accessibility_status=='resolved')&(nc.ancestry!='EUR')].variant_id.head(8).tolist())
