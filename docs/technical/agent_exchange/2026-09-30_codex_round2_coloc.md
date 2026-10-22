# Codex round-2 COLOC direction review, 2026-09-30

**Verdict: revise.** The v5 counts reproduce, but its marginal GWAS sign has not been checked
against the colocalized GWAS component. This limits the interpretation of the direction column.
No builder file or release authority was edited.

## Independent counts

CPU job 21985809 (1536, 1 CPU, 2 GB) completed successfully. The recount reads public direction
tables only, with no MASLD allelic outcome. Output:
[round2_summaries_21985809.json](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/round2_summaries_21985809.json).

| Release | Trait rows | Genes | Increase | Decrease | Not directional | Tag r² ≥ 0.8 genes |
|---|---:|---:|---:|---:|---:|---:|
| adopted_2026-08-17 | 628 | 462 | 217 | 220 | 191 | 34 |
| candidate_r2_20260924T053244Z | 589 | 436 | 212 | 213 | 164 | 34 |

There are zero duplicate gene×trait keys within either release and zero v3→v5 sign flips.
All 437/425 directional rows pass both p < 1e-3 floors; 407/396 also pass p < 1e-5.
Every reported directional row has an available eQTL fit and agreement between its eQTL
component and marginal signs. Gene sign patterns reproduce 160/152/15/135 and 156/147/14/119
(increase only/decrease only/both/none). Pair counts are 925/840; 85/0 are ineligible under the
corrected rule. Those historical adopted pairs are explicitly labelled. Counts remain separated
by release. The primary tag column is `t3prime_v2_kept`, as the lead requested.

## Review items

### C10-G1, major: the marginal floor cannot establish GWAS component direction

Evidence: [coloc_direction.py:206](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/coloc_direction/coloc_direction.py:206)
tests marginal p; line 258 multiplies marginal GWAS sign, LD sign and marginal eQTL sign;
line 286 checks only the eQTL component. No GWAS component sign is exported.
SuSiE COLOC compares pairs of signals at loci with multiple effects
([official COLOC documentation](https://chr1swallace.github.io/coloc/articles/a06_SuSiE.html)).

A direct counterexample uses two standardized SNPs with r = 0.8 and conditional effects
(+1, −2). The first marginal effect is 1 + 0.8×(−2) = −0.6. Its p-value can be arbitrarily
small with increasing n, while its conditional effect remains positive. A stricter marginal
floor therefore does not resolve the missing component check. [Unverified] The frequency of
such disagreements in these outputs is unknown; this review does not claim an observed reversal.

Requested fix: obtain the matching GWAS component effects and report their agreement with the
oriented marginal signs. If those fits remain unavailable, label the current estimate explicitly
as a marginal trait-associated haplotype direction, mark conditional GWAS direction unverified,
and exclude it from signal-specific direction claims. The p floor is useful evidence of a
marginal association, but is not a substitute for this check.

### C10-G2, minor: N4BP2L2 passes the written rule but is sensitive to the MAF cap

Evidence: N4BP2L2 ALT, UKBB_ALT, 13:33198176, T/A, has panel frequency 0.39961759082217974,
no source frequency, and `gwas_allele1_follows_file_convention=True`. Its GWAS and eQTL p-values
are 4.706e-10 and 2.047322e-124, and direction is −1 in both releases.
[The palindrome rule at line 165](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/coloc_direction/coloc_direction.py:165)
permits the empirical allele-order convention when source frequency is absent. This row passes
the prespecified 0.4 cap on the selected 523-founder panel. The builder reports that the
525-founder frequency was 0.4010, on the other side of the cap.

Requested fix: retain this disclosed sensitivity with the row and avoid using N4BP2L2 as an
unqualified direction example. Do not change the cap after inspecting this result. The allele
order supports the coding convention; the near-0.4 frequency alone is weak strand evidence.

### C10-G3, minor: call the allele-order denominator pair records

Evidence: the recount finds 1,239/1,239 non-palindromic SNV pair records following the convention.
There are 165 indel pair records, with all ten convention disagreements among indels.
Records can repeat a study/lead across signals or releases.

Requested fix: describe 1,239 as pair records and additionally report distinct study×lead
counts if claiming the number of independently checked leads. No directional SNV convention
failure was found.

## Checks and limits

Python compilation and launcher `bash -n` passed. Two earlier recount jobs (21985807/21985808)
stopped because the reviewer initially treated an absent development-cohort verdict as an
ordinary G1 verdict; the corrected recount recognizes `not_computed` with zero sites.
No raw GWAS refit, LD-panel reconstruction or independent allele-frequency recalculation was
run in this round. The public tables and source code support the checks above. The corrected
genetic release remains a candidate.

**Verdict: revise**, principally for C10-G1.
