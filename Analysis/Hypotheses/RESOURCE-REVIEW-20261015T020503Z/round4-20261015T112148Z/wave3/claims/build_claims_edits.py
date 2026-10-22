"""Build claims_edits.json for round-4 claim-level text (title, abstract, intro, Discussion, PAPER, README, STATUS).

Every old_text is taken from ST/stage0 and checked to occur exactly once.
"""
import json
import re
from pathlib import Path

ST = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/"
          "RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/staged-20261015T184532Z/stage0")
OUT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/"
           "RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/wave3/claims/claims_edits.json")

H = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses"
R4 = f"{H}/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z"
SRC = {
    "H57": f"{H}/H57-20261015T112421Z/run-20261015T113107Z/RESULT_H57.md",
    "H53": f"{H}/H53-20261015T112655Z/RESULT_H53.md",
    "H55": f"{H}/H55-20261015T112558Z/RESULT_H55.md",
    "H58": f"{H}/H58-20261015T112423Z/jobs/5800/H58-20261015T113519Z_22252009/RESULT_H58.md",
    "H59": f"{H}/H59-20261015T112432Z/RESULT_H59.md",
    "H60": f"{H}/H60-20261015T112434Z/out-final-20261015T120955Z/RESULT_H60.md",
    "H71": f"{H}/H71-20261015T124239Z/run-20261015T125815Z/RESULT_H71.md",
    "C7": f"{R4}/wave3/c7/20261015T185127Z/c7_results.tsv",
    "R9C": f"{R4}/catalog_r9b_candidate-20261015T185230Z/build_out/catalog_r9b_candidate_gene_table.parquet (fields carried to r9c)",
    "RESULTS": "ST/stage0 docs/RESULTS.md:55-58 (774-unit refit) and :1853-1860 (Stromal ECM fibrosis table)",
    "PAPER": "ST/stage0 docs/PAPER.md:316-320 (COLOC r2 436/996, 50 studies), :546-590 (H58 values, 19 direct clusters, 423/13)",
    "TITLE": f"{R4}/TITLE_CLAUSE_DECISION.md",
    "NOV": f"{R4}/wave3/novelty/novelty_paragraph.md",
    "CLAIMS": "ST/Claims.json R13-02 check (15.81% / 3.26% credible-set weighted)",
    "PROGCTX": f"{R4}/catalog_r9b_candidate-20261015T185230Z/build_out/sidecar_program_context_r9.tsv",
}

edits = []


def text(rel):
    return (ST / rel).read_text()


def lines(rel, a, b):
    """Exact text of lines a..b (1-based, inclusive) without the final newline."""
    ls = text(rel).split("\n")
    return "\n".join(ls[a - 1:b])


def add(rel, old, new, claim, source, rationale):
    n = text(rel).count(old)
    assert n == 1, (rel, n, old[:80])
    assert "—" not in new, ("em dash", rel, new[:80])
    assert "–" not in new, ("en dash", rel, new[:80])
    assert not re.search(r"\[(Inference|Unverified|Speculation|Estimate)", new), rel
    edits.append(dict(file=rel, old_text=old, new_text=new, claim_level=claim,
                      source=source, rationale=rationale))


TITLE_LINE = ("A multimodal human MASLD Resource records liver-trait genetic evidence and\n"
              "disease-state liver expression side by side for each gene")

# ---------------------------------------------------------------- title.md
T = "docs/manuscript/draft/title.md"
historical = lines(T, 27, 38)
assert historical.startswith("The historical analysis reported 2,426")
new_title = f"""## Binding title

**{TITLE_LINE}**

Adopted 2026-10-15 with the user's approval, replacing "A multimodal human
MASLD Resource resolves the architecture of inherited risk and the
coordination of multicellular remodeling" (adopted 2026-09-07). The
coordination clause was removed under the joint H52/H54 rule (row: H52
untestable / H54 replicated); the user agreed on 2026-10-15. "Resolves the
architecture" was dropped because H58 found the direct-trait versus
liver-enzyme coding contrast confounded with GWAS ancestry and fine-mapping
coverage, so no trait-class call is made. The title states what the Resource
records. It does not assert a finding about either layer or a link between
them.

"Liver-trait genetic evidence" covers multi-signal and single-signal liver
eQTL colocalization, fine-mapped protein-altering or canonical-splice
credible-set members, and colocalized liver chromatin-accessibility regions.
The manuscript names which of these applies to a gene. None of them alone
identifies a causal gene, and the strongest colocalization evidence for 423 of
the 436 primary multi-signal genes comes from a liver-enzyme trait.

"Disease-state liver expression" means the pooled cohort-adjusted
disease-versus-control calls from 774 genetic-identity units in five cohorts,
and the cross-sectional fibrosis-stage contrasts. The units are genetic
clusters, not verified donors, and stage contrasts compare different people.

Where the text describes how fine-mapped posterior probability at liver-trait
loci divides between protein-altering, splice and noncoding variants, the
direct-trait versus liver-enzyme contrast is descriptive, rests on 19
direct-trait regions and is pooled across ancestries. It is confounded with
GWAS ancestry and fine-mapping coverage and receives no trait-class call in
any ancestry stratum. It does not claim a resolved causal variant at any
locus.

**Coordination is not a title claim.** The clause was removed under the joint
H52/H54 rule (row: H52 untestable / H54 replicated). H52 could not be read in
GSE213621 because its program-group permutation control failed. Program-pair
covariance is reported in Results as a discovery-cohort result, with the H54
replication in GSE281797 and GSE225740 stated beside it and never pooled. Bulk
program scores mix cell abundance with program activity, so no claim of
coordination between cell types follows. The paragraph below is historical.

{historical}

The predictive-model result is not in the title. It remains Figure 7."""
add(T, lines(T, 8, 41), new_title, True,
    f"{SRC['TITLE']}; {SRC['H58']}; {SRC['PAPER']}; {SRC['RESULTS']}; ST/stage0 fig5_discussion.md:228-232 (GSE281797, GSE225740)",
    "Resolves conflict 0 and overlap groups G08, G25, G26 in one edit. One title, with neither 'resolves the architecture' (H58 no call) "
    "nor the coordination clause (joint rule, H52 untestable). The Claims title 'describes the coding and noncoding architecture ... maps "
    "cross-sectional multicellular remodeling' was not used: it leads with the H58 contrast that now receives no call, and 'multicellular "
    "remodeling' rests on bulk scores and a single-cell layer with 0 of 117 stage calls (H59). E189 placeholders filled. 'keeps its "
    "abstract sentence' removed: the stage0 abstract has no Figure 7 sentence. The 74.9%/93.7% no-context figures of W0359 are not used "
    "because H58 attributes the unresolved part to failed lookup. Title change needs user approval at install (title.md lines 3-4).")

# ---------------------------------------------------------------- abstract.md
A = "docs/manuscript/draft/abstract.md"
abs_old = lines(A, 8, 69)
assert abs_old.startswith("Inherited genetic variation") and abs_old.endswith("resolve\nmechanism.")
abs_new = """Inherited genetic variation and the molecular state of established disease
have each been measured in metabolic dysfunction-associated steatotic liver
disease (MASLD), but gene-level evidence from the two is rarely kept side by
side with a record of what each assay could test. We assembled a human MASLD
Resource that combines liver bulk RNA sequencing, analyzed as 774
genetic-identity units (138 control, 636 disease) from five cohorts in a
pooled cohort-adjusted model, with liver-trait GWAS from 50 studies,
fine-mapping and multi-signal colocalization with liver expression QTL. The
bulk units are conservative genetic clusters, not verified people. The
disease-versus-control fit gave 1,854 canonical disease-state calls among
24,115 tested genes, and colocalization identified 436 primary multi-signal
genes and 996 broader genetic candidates. For 423 of the 436 genes the
strongest colocalization evidence came from a liver-enzyme trait, and for 13
from a direct disease or liver-fat trait. Among 3,769 genes evaluated by
multi-signal colocalization, 26 were both colocalized and canonical
disease-state calls, against 21.2 expected for genes of the same expression
level (odds ratio 1.29, 95% CI 0.84 to 2.01). Expression was not fully
controlled, so neither depletion nor enrichment was detected and no direction
is stated. Most genes called only in the disease state (935 of 1,129) had no
multi-signal colocalization result; they are recorded as untested, not as
genetic negatives. Reanalysis of a 138-donor liver chromatin-accessibility QTL
map (Currin et al., 2025) placed a colocalized accessible region at 191 of 581
liver GWAS signals, 45 of them reported by the source; at 108 Catalog genes
the same region also colocalized with the gene's expression QTL. In 222 GTEx
postmortem livers, the 1,002 measurable disease-state genes changed with
image-measured histologic fat in the same direction as in the biopsy
disease-versus-control contrast (Spearman rho 0.638), and the agreement held
in the 138 donors who did not die on a ventilator; the cause of the fat is not
recorded. Donor-aware single-cell analysis placed 117 fixed gene programs in
five cell lineages; none passed the 117-program family test for a
cross-sectional stage association. The MASLD Gene Catalog gives 28,233 gene
identities separate evidence states for eQTL colocalization, fine-mapped
protein-altering variants, colocalized accessible regions, disease-state RNA
and postmortem fat response, with liver and plasma protein detection,
single-cell lineage detection and healthy-liver zonation. It distinguishes a
layer that was untested from one that was tested without support, and assigns
each gene the next discriminating experiment. This Resource records
inherited-variation and disease-state evidence side by side for each gene,
with the experiment needed to connect them."""
add(A, abs_old, abs_new, True,
    f"{SRC['RESULTS']}; {SRC['PAPER']}; {SRC['H57']} (T1, D6 935/194 of 1,129); {SRC['H53']} (Catalog deliverable 191/581, 45 Currin, M1 108); "
    f"PRESPEC_H53 (138 donors); {SRC['H55']} (rho 0.638, 222 donors, 1,002 genes, Q1 n 138); {SRC['H59']} (0 of 117); {SRC['R9C']} (28,233 rows; fields)",
    "Full rewrite of the abstract body. Resolves G06 and G23 and conflict 27 (last sentence). Keeps only the results the round-4 rule allows: "
    "774-unit calls, COLOC r2 counts, H57 with the S6 label, H53, H55, Catalog fields; the H59 null is stated as a null. Removed: the "
    "pooled 6.03 coding ratio and 57-comparison menu (H58 no trait-class call), the hepatocyte-specificity comparison (failed prediction, "
    "post hoc only; moved to intro and Discussion with the conflict-7 rerun), the GNMT/MAT1A/CYP2C19 sentence (not a round-4 passed "
    "result; it also said SuSiE was 'not run' for GNMT, which TP-F3cdg-00 shows is wrong), the continuum sentence, the Visium "
    "ECM/IGFBP7 sentence, and 'reusable map from inherited variation to multicellular MASLD biology'. 'Regulatory DNA elements' is no "
    "longer listed as a Catalog object (R07-01).")

# ---------------------------------------------------------------- 01_intro.md
I = "docs/manuscript/draft/01_intro.md"
add(I,
    """None of these combines fine-mapped multi-ancestry
colocalization with multicohort disease-state expression and donor-aware cell
and tissue context. What remains missing is a human MASLD Resource that
connects inherited variation to multicohort disease remodeling and then places
both in donor-aware cell and tissue context. This connection must preserve the
distinct biological meaning, coverage, provenance, and inferential unit of each
assay rather than collapsing them into one modality count.""",
    """None of these keeps fine-mapped multi-ancestry
colocalization, multicohort disease-state expression and donor-aware cell and
tissue context side by side for each gene. What remains missing is a per-gene
record of both kinds of evidence that also states, for each layer, whether the
gene was untested or tested without support. Such a record must preserve the
distinct biological meaning, coverage, source, and inferential unit of each
assay rather than collapsing them into one modality count.""",
    True, f"{SRC['NOV']} (not_novel list: Piras 2024, Sveinbjornsson 2022)",
    "The gap statement now names a per-gene side-by-side record, not a 'connection' the corrected results do not show (H57 overlap "
    "near expectation, H59 null). 'provenance' replaced with 'source'.")

add(I, lines(I, 58, 71),
    """We assembled a multimodal human MASLD Resource to keep these two kinds of
evidence side by side. The transcriptomic arm includes bulk RNA sequencing
across nine cohorts, with 774 selected genetic-identity clusters from five
control-bearing cohorts entering a pooled, cohort-adjusted disease-versus-control
analysis. These clusters are not verified independent people. The genetic arm
integrates prespecified GWAS strata across ancestry groups with liver eQTL
colocalization, fine-mapped coding and noncoding context, and a reanalysis of a
published liver chromatin-accessibility QTL map (Currin et al., 2025). Direct
disease and imaging phenotypes remain distinct from liver-enzyme proxy traits.
Donor-aware single-cell analyses score prespecified gene programs in each cell
lineage, while protein, chromatin-accessibility, and spatial assays provide
assay-native molecular and physical context. We retain fibrosis stage in the
primary severity analyses and add an externally anchored transcriptomic
continuum as a complementary layer for describing molecular variation among
participants assigned the same stage.""",
    True, f"{SRC['RESULTS']}; ST/stage0 fig1.md:16 (nine cohorts); {SRC['H53']}",
    "'to make this connection' becomes 'to keep these two kinds of evidence side by side'; adds the H53 caQTL reanalysis with source "
    "credit; 'testing coordinated molecular variation' becomes 'describing molecular variation' (coordination clause removed).")

add(I, lines(I, 73, 87),
    """The two maps can be compared gene by gene. Among 3,769 genes evaluated by
multi-signal colocalization, the 426 colocalized and 220 disease-state genes
share 26, against 21.2 expected for genes of the same expression level
(Mantel-Haenszel odds ratio 1.29, 95% CI 0.84 to 2.01). Expression was not
fully controlled, so no direction is stated, and the overlap is not itself
informative: the interval includes 1, odds ratios between about 0.46 and 1.77
had less than 80% power to be detected, and two sets of these sizes are mostly
disjoint whether or not the underlying biology differs. Most genes called only
in the disease state (935 of 1,129) had no multi-signal colocalization result;
they are untested, not genetic negatives. Matched for expression, colocalized
genes are not more hepatocyte-specific than disease-state genes, so the
prespecified directional prediction failed. In a post hoc comparison the mean
genetic-minus-state difference ran the other way (−0.265515 log2 units in 281
matched pairs; 126 positive, 155 negative), and it kept that direction when
restricted to state-only genes that colocalization had tested (−0.342 log2
units in 129 matched pairs). No difference is detected in liver protein
response to fibrosis stage (73 pairs, nominal sign-flip P = 0.487465); this is
not equivalence. These comparisons do not establish formal biological class
separation.""",
    True, f"{SRC['H57']} (T1, G3 power, S6, D6, D9); {SRC['C7']} (arms A, C)",
    "Adds the H57 S6 label (TP-F1-00) and resolves conflict 7 with the completed rerun (job 22257837): the restricted comparison keeps "
    "the direction. The post hoc P values stay out of the introduction, as the PAPER rule and R14-10 require; they are in the Discussion.")

add(I, lines(I, 89, 101),
    """Fine-mapping the trait-associated regions describes how posterior
probability divides between variant classes. Averaged over one-Mb regions and
pooled across ancestries, credible sets for direct MASLD and liver-fat
phenotypes carry 23.34% of their posterior probability in protein-altering or
canonical-splice variants, against 3.87% for the liver enzymes ALT, AST and
GGT (15.81% and 3.26% when each credible set is weighted equally). The pooled
contrast mixes trait class with GWAS ancestry and fine-mapping coverage. Within
European studies the fractions are 29.88% and 6.92%, within AMR 6.73% and
4.54%, and no stratum supports a calibrated trait-class call. Separately, the
strongest evidence for 423 of 436 colocalized genes comes from a liver-enzyme
trait and for 13 from a direct disease or liver-fat trait, and 22 genes have
any direct-trait colocalization. Far more enzyme than direct-trait
gene-by-GWAS pairs were tested by multi-signal colocalization (12,121 against
771). Some protein-altering variants may act through protein function without
changing steady-state RNA, while others change RNA abundance through
nonsense-mediated decay or splicing (Rivas et al., 2015). The Catalog
therefore keeps fine-mapped protein-altering variants in a field of their own,
beside expression colocalization.""",
    True, f"{SRC['H58']}; {SRC['RESULTS']} lines 1-8 (EUR/AMR); {SRC['CLAIMS']}; ST/stage0 abstract.md:33-37 and fig5_discussion.md:131-135 (22 genes; 12,121 vs 771); {SRC['R9C']} (coding_cs_state)",
    "Resolves G07: merges W0036 (H58 stratification, no call) and W0358 (region-weighted vs credible-set-weighted values). Drops the "
    "ratio with its bootstrap interval, 'largely a map of liver-enzyme regulation', and the hypothesis sentence that leaned on the "
    "trait-class coding difference (TP-F2d-00).")

add(I,
    """ECM/IGFBP7 and ductular-injury/BICC1 test whether a disease-associated program,
its cell-of-origin assignment, and its intact-tissue organization support the
same conclusion.""",
    """ECM/IGFBP7 and ductular-injury/BICC1 test whether a program's disease
association, its cell-of-origin assignment, and its intact-tissue organization
support the same conclusion. For ECM/IGFBP7 the bulk fibrosis association does
not reach BH q<0.05 in the 364 units of the three 774-unit cohorts with staged
disease (q=0.214; q=0.421 after composition adjustment); it reaches q=0.045
only when 93 GSE174478 source records are added, and q=0.230 after composition
adjustment in that arm.""",
    False, "ST/stage0 docs/RESULTS.md:1853-1860 (Stromal ECM rows: 364 q 0.214/0.421; 457 q 0.045/0.230)",
    "Resolves G12 (W0150, W0210). Uses only values printed in the RESULTS table; W0150's '224 program-by-axis tests' and '0.23 to 0.25 "
    "depending on the fraction set' were not checked against a file and are left out.")

add(I, lines(I, 118, 127),
    """We therefore placed these dimensions in donor-resolved cellular programs,
liver protein, accessible chromatin, and intact-tissue organization. Each
assay retains its native effect unit and biological unit. Coverage,
applicability and donor identity determine whether an apparent absence is
interpretable. The resulting MASLD Gene Catalog gives each protein-coding gene,
lncRNA transcript and unresolved locus separate evidence states for inherited
variation and disease-state RNA, adds gene-level fields for colocalized
accessible regions, fine-mapped protein-altering variants, and protein and
single-cell detection, links genes to cell programs by membership, and assigns
the next discriminating experiment; spatial context is attached at program
level. This per-gene record, which keeps the two kinds of evidence apart and
separates untested from unsupported, is the standalone Resource contribution.""",
    True, f"{SRC['R9C']} (caqtl_*, coding_cs_*, liver_/plasma_proteome_*, sc_detected); {SRC['PROGCTX']}",
    "Resolves G24 (W0290, W0371): no regulatory DNA element objects, no program-member gene count (3,422 in W0290 vs 3,498 in the "
    "R07-01 note disagree). Gene-level fields follow the r9c candidate so the introduction matches the abstract; 'provenance' removed.")

# ---------------------------------------------------------------- fig5_discussion.md
D = "docs/manuscript/draft/fig5_discussion.md"
add(D, lines(D, 90, 93),
    """Across Figures 1 to 7, the Resource records inherited genetic evidence and
disease-state remodeling for each gene and places both in cell and tissue
context where an assay can observe them. The Gene Catalog is the reusable
endpoint: it preserves the biology, evidence boundary, and next experiment for
each nominated molecular object.""",
    True, "wording only; no new number",
    "Same overreach as conflict 27: the corrected results do not show a connection between the two maps.")

add(D, lines(D, 97, 100),
    """This Resource records liver-trait genetic evidence and cross-sectional
disease-state remodeling in established MASLD for each gene, and keeps their
assay boundaries separate. Sparse membership overlap does not establish
biological separation, equivalence or a privileged intersection of the two
maps.""",
    True, "wording only; no new number", "Opening claim matches the new title.")

add(D,
    """One hypothesis about which genes appear in
both sets is that expression-based colocalization detects the coding component
of direct-trait risk poorly; the present comparison does not test it, because
it does not hold gene eligibility, phenotype, study power and expression
observability fixed.""",
    """One hypothesis about which genes appear in
both sets is that expression-based colocalization misses risk variants that act
through protein function rather than RNA abundance; the present comparison does
not test it, because it does not hold gene eligibility, phenotype, study power
and expression observability fixed.""",
    True, "wording only; no new number",
    "TP-F2d-00: the old hypothesis leaned on a trait-class coding difference that H58 leaves without a call.")

anchor = """This is a selection-conditioned
descriptive difference, not confirmatory class separation or independent
biological replication."""
add(D, anchor, anchor + """ Most state-only genes (935 of 1,129) had no
multi-signal colocalization result, so the state-only pool is mostly genes that
colocalization could not test; the 193 single-cell-observable state-only genes
it did test have mean AveExpr 1.3405, against 0.4164 for all 1,124
single-cell-observable state-only genes. Restricted post hoc to those 193
genes, the difference kept its direction: in 129 expression-matched pairs
(expression SMD −0.0004) the mean genetic-minus-state difference was −0.342
log2 units (median −0.183, paired d −0.217; nominal sign-flip P=0.0154,
100,000 draws, seed 20260923; Wilcoxon P=0.0314; 57 positive and 72 negative
pairs, sign test P=0.218). The design had 80% power for d=0.247. Within the
281 published pairs, the difference comes mainly from the 56 pairs whose state
partner was tested (mean −0.748, nominal P=0.00148) rather than the 225 with an
untested partner (−0.145, P=0.233). Genes and pairs are the units here, the
sign-flip null does not account for genes sharing loci, and none of these P
values is calibrated.""",
    True, f"{SRC['C7']} (arms C, D); {SRC['H57']} D6/D7",
    "Resolves conflict 7 with the completed restricted rerun. All values read from c7_results.tsv (job 22257837).")

add(D, "(Gobeil et al., 2026b)",
    "(Gobeil et al., 2026b; medRxiv doi:10.64898/2026.06.04.26354903)",
    False, "bioRxiv/medRxiv API get_preprint 10.64898/2026.06.04.26354903 (checked 2026-10-15: Gobeil E et al., 504 liver biopsies, proteogenomic, v1 2026-06-06, not published)",
    "Novelty note 4 asked whether Gobeil 2026b exists. It does, as a medRxiv preprint; the DOI is added so the citation can be checked.")

add(D,
    """Program-pair covariance seen in
the discovery cohorts recurred in two further bulk cohorts (GSE281797 and
GSE225740); their participants were not checked by genotype against the
discovery cohorts, and bulk covariance mixes cell abundance with program
activity. Neither""",
    """Program-pair covariance seen in
the discovery cohorts recurred in two further bulk cohorts (GSE281797 and
GSE225740) and could not be tested in GSE213621, where the program-group
permutation control failed. The replication participants were not checked by
genotype against the discovery cohorts, and bulk covariance mixes cell
abundance with program activity, so it does not show coordination between cell
types. Neither""",
    True, f"{SRC['TITLE']} (H52 untestable, NC2 failed; H54 replicated)",
    "States H52 beside H54 and the cell-type limit, so the Discussion does not imply the removed title clause.")

nov_par = Path(SRC["NOV"]).read_text().strip()
assert nov_par.startswith("Each component of this Resource has precedent.")
add(D, lines(D, 237, 256),
    "\n".join([
        nov_par,
        "",
        "The 117 programs and 43 systems used here were fixed before any outcome was",
        "read and are never reweighted by network position, so a sparsely annotated",
        "gene keeps its assigned weight and an inadequately observed one is recorded",
        "as untestable rather than ranked low.",
    ]),
    True, f"{SRC['NOV']} (paragraph and citations as given); {SRC['H57']}; {SRC['H71']}; {SRC['H53']}",
    "Replaces the old positioning paragraph as the novelty note suggests. Removes the statement that the earlier meta-analyses and the "
    "continuum framework 'carry no fine-mapped inherited-variation layer' (not_novel list: Piras 2024 integrated GWAS) and the "
    "unchecked claim that Kumak hub lists are almost all previously published genes. The fixed-program design point is kept.")

anchor = ("No accepted spatial branch provides a MASLD-versus-healthy disease-effect "
          "estimand, and this paper contains no new experimental perturbation.")
add(D, anchor, anchor + (
    " The GTEx comparison measures histologic fat of unrecorded cause in postmortem liver; its"
    " slide tissue-area control came close to its stop rule (Spearman rho 0.354, P=0.0564), and"
    " ischemic time also tracked the same genes (rho −0.350), although less strongly than fat."
    " The caQTL reanalysis rests on one 138-donor map; its chr1 and chr19 positive controls were"
    " seen in a pilot and are not blind, 45 of the 191 colocalized signals were already reported"
    " by its authors, and overlap with published three-way eQTL, caQTL and GWAS colocalizations"
    " (Broadaway et al., 2024) was not checked."),
    True, f"{SRC['H55']} (NC1, NC2, Q2); {SRC['H53']}; PRESPEC_H53 (pilot controls, 138 donors); {SRC['NOV']} (Broadaway PMID 39173627)",
    "Limits for the two new abstract-level results.")

add(D,
    "Preserving these distinctions—and releasing the provenance needed to audit them—is more useful than interpreting a single aggregate rank as target truth.",
    "Keeping these distinctions, and releasing the source tables needed to check them, is more useful than reading a single aggregate rank as target truth.",
    False, "wording only", "Removes two em dashes and 'provenance ... audit'.")

# ---------------------------------------------------------------- 05_figure_legends.md (G22)
L = "docs/manuscript/05_figure_legends.md"
add(L, "Protein, chromatin and spatial context retain their assay or program unit.",
    """Chromatin and spatial context are attached at program level for two programs
(Stromal ECM and Ductular injury), not as gene-level calls; liver protein
(PXD051911) was untestable for both. No record is a regulatory DNA element
object.""",
    False, f"{SRC['PROGCTX']} (14 rows; PXD051911 liver_DIA_MS untestable for both programs in r8 and r9)",
    "Resolves G22 (W0288, W0373). The protein state was checked in the program-context table.")

# ---------------------------------------------------------------- PAPER.md
P = "docs/PAPER.md"
add(P, lines(P, 66, 67), f"**{TITLE_LINE}**", True, SRC["TITLE"],
    "Resolves G27 (W0353, W0354) with the same title as title.md; needs user approval at install.")

add(P,
    """Their overlap (26 genes against 21.2
   expected within expression deciles, among 3,769 genes evaluated by
   multi-signal COLOC) does not establish biological complementarity or
   antagonism.""",
    """Their overlap (26 genes against 21.2
   expected for genes of the same expression level, among 3,769 genes
   evaluated by multi-signal COLOC; expression not fully controlled, so no
   direction is stated) does not establish biological complementarity or
   antagonism. Most state-only genes (935 of 1,129) have no multi-signal COLOC
   result and are untested, not genetic negatives.""",
    True, f"{SRC['H57']} (T1, S6, D6)", "TP-F1-00 label and the conflict-7 caveat.")

add(P, "cell lineages; no single-cell program has a family-level stage association.",
    "cell lineages; no single-cell program passes the 117-program family test for\n   a stage association (H59).",
    True, SRC["H59"], "TP-F4abcd-01 wording.")

add(P,
    """The paper connects **inherited genetic variation** to **multicellular disease
remodeling** in human MASLD.""",
    """The paper records **liver-trait genetic evidence** and **disease-state
remodeling** in human MASLD side by side for each gene.""",
    True, "wording only", "Matches the new title; the corrected results do not show a connection.")

add(P,
    """headline. The impact comes from connecting inherited variation to
multicellular remodeling and making both maps usable through cell and tissue
context, a shared set of evidence rules, and the Gene Catalog.""",
    """headline. The contribution is a per-gene record that keeps inherited
genetic evidence and disease-state remodeling side by side, states whether each
layer was untested or tested without support, and makes both usable through
cell and tissue context, a shared set of evidence rules, and the Gene Catalog.""",
    True, "wording only", "Same reason.")

add(P, "where SuSiE was not run for\nGNMT)",
    "where SuSiE was attempted for\nGNMT but returned no credible-set result)",
    True, "ST/text_problems.json NR-F3cdg-R14-09 (method abf_fallback, lambda_s 0.00123)",
    "TP-F3cdg-00: SuSiE was attempted, not 'not run'.")

add(P,
    """- The Resource maps trait-associated inherited genetic variation and
  multicellular remodeling in established MASLD, then connects both to cell and
  tissue context.""",
    """- The Resource records liver-trait genetic evidence and cross-sectional
  disease-state remodeling in established MASLD side by side for each gene, and
  places both in cell and tissue context where an assay can observe them.""",
    True, "wording only", "First claim line matches the title.")

anchor = "  The earlier 272-pair, −0.44, p=2×10⁻⁵ statement is historical."
add(P, anchor, anchor + """
  Most state-only genes (935 of 1,129) have no multi-signal COLOC result.
  Restricted post hoc to the 193 single-cell-observable state-only genes that
  COLOC tested, the difference keeps its direction (129 expression-matched
  pairs, −0.341885 log2 units, nominal sign-flip p=0.01543, seed 20260923;
  sign test 57/72, p=0.218; MDE d 0.247). Within the 281 published pairs, the
  56 with a tested state partner carry most of the difference (−0.748,
  p=0.00148) and the 225 with an untested partner do not (−0.145, p=0.233).
  These P values are nominal and do not account for genes sharing loci; none is
  printed in the abstract or introduction (conflict-7 rerun, job 22257837,
  `Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/wave3/c7/20261015T185127Z/c7_results.tsv`).""",
    True, SRC["C7"], "Resolves conflict 7 in the claim ledger.")

add(P, "stage by resolving coordinated cross-sectional transcriptional remodeling",
    "stage by describing cross-sectional transcriptional variation",
    True, "wording only", "Drops 'coordinated' after the title clause was removed.")

anchor = """- The Resource can prospectively define hypotheses for a separate functional
  screen."""
add(P, anchor, """- Reanalysis of a published 138-donor liver caQTL map (Currin et al. 2025,
  PMID 40393811) places a colocalized accessible region at 191 of 581 liver
  GWAS study-lead signals; 45 of the 191 were already reported by the source,
  so this is a credited reanalysis, not a discovery. At 108 Catalog genes the
  same region also colocalizes with the gene's eQTL. Among loci selected for
  liver-eQTL colocalization, liver-only loci share a liver caQTL more often than
  platelet-only loci at equal GWAS and caQTL strength (adjusted OR 3.350, 95%
  CI 1.622 to 7.397; 172 against 78 loci), which agrees with the source's
  trait-level pattern (H53).
- In 222 GTEx postmortem livers, the 1,002 measurable canonical disease-state
  genes change with image-measured histologic fat in the same direction as
  their biopsy disease-versus-control logFC (Spearman rho 0.638, Freedman-Lane
  p 0.0002 at the permutation floor), also in the 138 donors outside the
  ventilator death class (rho 0.632), and more strongly for fat than for
  ischemic time (H55). The claim is limited to histologic steatosis of
  unrecorded cause in postmortem liver; the tissue-area control came close to
  its stop rule (rho 0.354, p 0.0564).
- The Gene Catalog keeps, for each gene, whether a layer was untested or tested
  without support. For the Kumak et al. (2026) lists, 661 of 902 NAS genes and
  714 of 943 fibrosis genes have no multi-signal COLOC result, and among
  evaluated genes their COLOC rate does not differ detectably from that of genes
  of the same expression level (H71).
""" + anchor,
    True, f"{SRC['H53']}; PRESPEC_H53 (138 donors); {SRC['H55']}; {SRC['H71']}",
    "Adds the two round-4 results allowed at abstract level and the H71 untested-versus-negative point.")

anchor = "| epigenetic causation or mechanism from accessibility | assay-native chromatin context, with peak-level evidence state |"
add(P, anchor, anchor + """
| "coordination of multicellular remodeling" or coordination between cell types | bulk program-pair covariance reported in Results with its replication cohorts; bulk scores mix cell abundance with program activity |
| a trait-class difference in coding or splice architecture | pooled descriptive fraction, confounded with GWAS ancestry and fine-mapping coverage; no trait-class call (H58) |
| a stage-associated single-cell program | 0 of 117 programs pass the family rule (H59) |
| direction agreement between ALT TWAS and disease-state expression | genome-wide signed agreement bounded inside ±0.03 without composition adjustment; no direction (H60) |
| liver caQTL colocalization at GWAS signals as this paper's finding | credited reanalysis of Currin et al. 2025; 45 of 191 signals already reported |""",
    True, f"{SRC['TITLE']}; {SRC['H58']}; {SRC['H59']}; {SRC['H60']}; {SRC['H53']}",
    "Records the round-4 negative and untestable outcomes as claims we cannot make.")

anchor = "- a gene-level Catalog that records context, limitations, and suitable follow-up."
add(P, anchor, anchor + """

Each component has precedent, listed in the Discussion: staged bulk
transcriptomes and their meta-analysis with gene-level GWAS (Piras 2024),
liver caQTL colocalization with GWAS signals (Currin 2025), three-way eQTL,
caQTL and GWAS colocalization (Broadaway 2024), coding plus expression gene
nomination (Sveinbjornsson 2022), and assignment of MASH risk loci to cell
types (Elison 2025). None of these is claimed as new here.""",
    True, SRC["NOV"], "Positioning follows the novelty search.")

# ---------------------------------------------------------------- README.md / STATUS.md
RM = "docs/README.md"
add(RM,
    """**Title:** *A multimodal human MASLD Resource describes the coding and
noncoding architecture of liver-trait genetic signals and maps
cross-sectional multicellular remodeling.* (Adopted
2026-09-07; the coordination evidence is conditional on resolving donor identity
in the already-used GSE193066 evaluation cohort.)""",
    f"""**Title:** *{TITLE_LINE}.* (Adopted 2026-10-15 with the user's approval. The
2026-09-07 coordination clause was removed under the joint H52/H54 rule, and
"resolves the architecture" was dropped after H58.)""",
    True, SRC["TITLE"], "Title line matches title.md and PAPER.md.")

SM = "docs/STATUS.md"
add(SM,
    """The binding title and framing are **A multimodal human MASLD Resource describes
the coding and noncoding architecture of liver-trait genetic signals and maps
cross-sectional multicellular remodeling**. The coordination clause was removed
under the joint H52/H54 rule.""",
    f"""The binding title and framing are **{TITLE_LINE}**. The coordination clause was removed under the joint H52/H54
rule (H52 untestable, H54 replicated), and "resolves the architecture" was
dropped after H58.""",
    True, SRC["TITLE"], "Title line matches title.md and PAPER.md.")

# sequential application check: every edit must still apply after the earlier ones in the same file
by_file = {}
for e in edits:
    by_file.setdefault(e["file"], []).append(e)
for rel, es in by_file.items():
    t = text(rel)
    for e in es:
        assert t.count(e["old_text"]) == 1, ("sequential", rel, e["old_text"][:60])
        t = t.replace(e["old_text"], e["new_text"])
    for bad in ["reusable map from inherited", "resolves the architecture of inherited\nrisk and the coordination"]:
        if rel.endswith(("abstract.md", "PAPER.md")):
            assert bad not in t, (rel, bad)

OUT.write_text(json.dumps(edits, indent=1, ensure_ascii=False))
print(len(edits), "edits written to", OUT)
for rel, es in by_file.items():
    print(rel, len(es))
