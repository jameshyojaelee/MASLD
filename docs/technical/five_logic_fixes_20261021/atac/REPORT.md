# ATAC indel normalization correction, 2026-10-21

Changed source: `Analysis/Multimodal_Program_Projection/scripts/atac_context_v3/genetics_context.py`.
Added regression tests: adjacent `tests/test_indel_normalization.py`.

The old left-alignment loop required both trailing allele bases to equal the preceding reference base. A valid rotation requires shared trailing allele bases; its preceding anchor may differ. The corrected loop therefore permits movement from a homopolymer/tandem repeat onto the preceding different base. Chromosome-start anchors remain nonempty. Nonpositive input positions are rejected before a negative reference fetch. Existing first-matching allele orientation is preserved and its ambiguity is now documented. No source-reference match or SNP/MNV normalization rule otherwise changed.

## Exact current-result impact

**Zero current numerical outputs change.** Complete streaming enumeration read 22,356,697 posterior rows from the original raw liftover input (558,299,548 compressed bytes). All 767 current r2 retained primary signal pairs were observed, representing 3,475,640 posterior rows and 3,835 lineage rows. Of these, 198,730 have length-changing alleles. 198,727 uniquely lifted, standard-chromosome indel rows were normalized using both the before-snapshot implementation and corrected implementation against the exact frozen hg38 FASTA. All normalized (position,ref,alt,orientation) tuples are identical. The remaining three indel rows retain the same upstream mapping exclusion. Retained inputs contain zero nonpositive hg38 positions. Equal-length allele behavior is identical by inspection of the narrow patch.

Consequently, all retained variant peak-membership flags, mapped/lost/shared/source-specific/union posterior masses, joint masses, primary-pair selection, evidence states and Figure5D numerical inputs are identical. This conclusion follows from equality of every applicable normalized variant tuple, rather than recomputing an unchanged downstream sum. No named release, current table, figure or result file was written.

Before/after tuple-stream SHA256 (both): `d5673044bb641f265134a8cfecebde803528d7845839d425009a3210e178ac43`.
`retained_indel_before_after.tsv.gz` records every compared indel occurrence, its input coordinates/alleles and before/after tuples.
`retained_normalization_changes.tsv` contains only its header because no retained variant changes.
`retained_indel_changes_summary.json`, `variant_scan_summary.json` record complete counts.
`fixture_before_after.tsv` shows the intended changed behavior on previously failing examples: reference CAAAAG, input pos3 AA>A, changes pos2 AA>A to leftmost pos1 CA>C; GATATATC pos4 ATA>A changes pos4 ATA>A to pos1 GAT>G. These equivalent representations preserve the edited haplotype.

## Verification executed

- New `test_indel_normalization.py`:6tests passed (homopolymer/tandem insertion/deletion with haplotype equivalence, chromosome start, SNP/MNV trimming, identifiable allele swap/reverse complement, unmatched alleles, invalid positions).
- Existing `test_contracts.py`:8tests passed.
- Independent comparison with installed bcftools 1.21: 6,723 exhaustively generated substitutions and indels across five small reference sequences, zero representation disagreements. VCF/FASTA, normalized output, comparison script and JSON are in this directory. Tested ALT lengths1-3, REF lengths1-3 at every feasible sequence coordinate, including chromosome start.
- `git diff --check` passed for the two changed files.
- Complete raw-input streaming scan53.5s; before/after retained-indel comparison11.8s on repeat,1CPU; cached normalization bounded65,536entries. No production pipeline, jobs, installation or figure rendering run.

Commands: `PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s Analysis/Multimodal_Program_Projection/scripts/atac_context_v3/tests -p test_indel_normalization.py -v`; analogous `test_contracts.py`. Scientific comparisons use existing `/gpfs/commons/home/jameslee/micromamba/envs/allelic_qtl/bin/python`,1BLAS thread,PYTHONNOUSERSITE=1.

## Interpretation limit

An unordered homopolymer allele pair can match the reference in both orders; ref-versus-alt direction is not identifiable from these arguments alone. This correction preserves the prior first-match convention rather than inventing a reference assignment. Swap tests cover cases where the reference orientation is identifiable. Current retained tuples are unchanged under that convention.
