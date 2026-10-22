# Codex independent check: T1 tag orientation and founder selection, 2026-09-29

Scope: public GTEx v8 liver eGene table, public 1kGP high-coverage phased panel and pedigree, and T1's public tag tables. No MASLD participant data, allelic outcomes, or sealed set were read. This is a candidate scientific check, not a Resource result.

## Answer

- **[Verified] Lead alleles:** among 5,127 autosomal GTEx liver eGenes with `qval <= 0.05` and a one-base REF/ALT lead, 5,065 leads match a biallelic 1kGP panel record exactly. None has REF/ALT swapped. Fifty-eight have no panel record at that position; four have a record there with different alleles (one is an indel). T1's exact `(position, REF, ALT)` join excludes those 62; the 62 cases are in `codex_t1_audit_20260929/nonidentical_leads.tsv`. No multi-allelic panel record was the sole explanation in this set.
- **[Verified] Stored LD orientation:** a seeded sample of 150 of the 4,390 T1 tag rows (seed 20260929; 50 sampled from negative-orientation rows, then 100 from the remainder, yielding 57 negative rows) matches independently recalculated EUR and EAS phased-haplotype r and MAF. All 150 meet both r² and MAF thresholds, both populations have the stored sign, and the largest absolute r discrepancy is `9.3e-15`. The GTEx-ALT to tag-allele to `a_o` chain agrees in the ten examples below. GTEx defines the eQTL effect allele as ALT rather than the minor allele: https://sites.google.com/broadinstitute.org/gtex-faqs/home (A10/A16).
- **[Verified] Founder-selection problem:** `t1_tags.py:22-26` retains 311 EUR and 371 EAS samples, but the filter removes every person listed as a parent while also requiring that their own parent fields be zero. It retains only 8/179 CEU and 7/157 IBS samples while retaining 99/99 FIN; in EAS it retains 13/163 CHS but 103/103 CHB. Five EUR and six EAS pairs with the same `FamilyID` still remain. The pedigree excludes recorded parent-child pairs among those retained, but does not establish that these 11 same-family pairs are genetically unrelated.
- **[Verified sensitivity; not a replacement panel]** Selecting one sample per `FamilyID` from the same public pedigree, preferring one with no recorded parents, gives 389 EUR and 437 EAS samples. Recalculating the same 150 tags preserves every LD sign; 143 still meet T1's two-population threshold. Seven sampled PPFIBP1 tag rows move just below the EUR r² cutoff: `-0.896230` (r² `0.803229`) with T1's set versus `-0.893038` (r² `0.797516`) in this alternative. They are seven rows from one gene, and PPFIBP1 has other retained tags; this does not establish a lost eGene. Maximum absolute r change is `0.01924`.

**[Inference]** The selected tag set depends at least locally on the pedigree filter and its strong subpopulation shift. Before freezing T1, choose and document a defensible unrelated reference set, then rerun T1/T3 and compare tag and eGene memberships. A one-per-family rule is only a sensitivity here: `FamilyID` alone cannot exclude relatedness across families, and its population mix is also uneven. The 150-row check does not quantify full-set changes or new tags that a revised reference might gain.

## Ten sign traces

`+` or `-` in the GTEx column means the effect of the **lead ALT** on expression. Positive LD puts the tag ALT on the lead-ALT haplotype; negative LD puts the tag REF there. The last column is the count used as `a_o`, the lead-ALT count. This checks orientation, not whether the in-biopsy effect matches GTEx.

| Gene | GTEx lead | Tag | GTEx slope sign | LD sign | Tag allele on lead-ALT | `a_o` source |
|---|---|---|---:|---:|---|---|
| PPFIBP1 | chr12:27641042 C>T | chr12:27636337 C>T | - | - | C | `ref_reads` |
| L3HYPDH | chr14:59537557 C>T | chr14:59474860 A>G | + | - | A | `ref_reads` |
| STXBP5 | chr6:147199829 T>C | chr6:147314277 A>G | + | - | A | `ref_reads` |
| ZNF175 | chr19:51571113 T>G | chr19:51573707 T>C | + | - | T | `ref_reads` |
| PNKD | chr2:218317758 G>A | chr2:218323178 C>T | - | - | C | `ref_reads` |
| ERCC6L2 | chr9:96024975 C>G | chr9:95966906 C>G | + | + | G | `alt_reads` |
| TIGD1 | chr2:232549950 G>A | chr2:232546284 C>T | + | + | T | `alt_reads` |
| RP11-499P20.2 | chr10:18530940 C>G | chr10:18531919 G>A | + | + | A | `alt_reads` |
| TMEM161B-AS1 | chr5:88292209 A>G | chr5:88295130 G>A | + | + | A | `alt_reads` |
| CAPN8 | chr1:223627873 G>A | chr1:223628875 T>C | - | + | C | `alt_reads` |

## Reproduction and limits

Input files: `data/external/allelic_refs/gtex_v8/GTEx_Analysis_v8_eQTL/Liver.v8.egenes.txt.gz`, `data/external/allelic_refs/kgp_highcov/1kGP_high_coverage_Illumina.chrN.filtered.SNV_INDEL_SV_phased_panel.vcf.gz`, `data/external/allelic_refs/kgp_highcov/20130606_g1k_3202_samples_ped_population.txt`, and `Analysis/MASLD_Model_Benchmark/executions/in-biopsy-allelic-tags-20260925T194248Z/tags/chr*.tags.tsv`. Producing code: `scripts/analysis/in_biopsy_allelic/tags/t1_tags.py:22-26,29-40,85-105`; consumer: `scripts/analysis/in_biopsy_allelic/tags/t5a_allelic_table.py:312-315`.

Scripts and sampled row-level outputs are in `docs/technical/agent_exchange/codex_t1_audit_20260929/`. The scripts run with `PYTHONNOUSERSITE=1` under the `rnaseq` Python environment and `bcftools/1.21 htslib/1.21`. The first audit queries panel alleles at every eligible lead and recomputes r by haplotype covariance. `bcftools query -S` emits genotypes in the sample-file order; `query -l -S` was not used to align genotype columns because it returned a different order for this panel. The second script tests the one-per-family sensitivity on precisely the first script's 150 sampled rows. No formal kinship matrix was available in the public inputs reviewed here, so cryptic relatedness is unresolved. The scripts and row-level outputs do not contain MASLD participant data.
