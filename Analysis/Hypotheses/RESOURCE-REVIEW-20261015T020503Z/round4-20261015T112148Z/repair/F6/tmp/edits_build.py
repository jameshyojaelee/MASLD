#!/usr/bin/env python3
"""Build the F6 repair text edits, check each old_text against the current file, write JSON."""
import json
import re
from pathlib import Path

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
OUT = Path(__file__).resolve().parent / "f6_repair_text_edits.json"
LEG = "docs/manuscript/05_figure_legends.md"
PAP = "docs/PAPER.md"
RES = "docs/RESULTS.md"
MET = "docs/manuscript/METHODS.md"
STA = "docs/STATUS.md"
DIS = "docs/manuscript/draft/fig5_discussion.md"
ABS = "docs/manuscript/draft/abstract.md"
INT = "docs/manuscript/draft/01_intro.md"
RDM = "figures/main/fig6_gene_catalog/README.md"

E = {}

E["R07-02"] = [
    (LEG, False,
     """**(B)** Genetics-call by RNA-call matrix over 28,233 selected r8 gene identities.
Supported, indeterminate and untestable states remain separate; zero strict
`tested_negative` calls occur. The 708 interface-untested genes with one
supported domain (698 RNA and 10 genetics) demonstrate why missing eligibility
is not negative evidence.""",
     """**(B)** Genetics-call by RNA-call matrix over 28,233 selected r8 gene identities.
Supported, indeterminate and untestable states remain separate. No negative-call
rule was defined for either domain, so `tested_negative` was not assessed and the
matrix has no tested-negative row or column. RNA `indeterminate` includes 11,301
genes with BH q < 0.05 whose absolute log2FC is at most 0.5. Genetics `untestable`
means no Tier-1/2 COLOC record (9,312 genes) or SuSiE-untestable for insufficient
shared posterior in one Tier-1/2 study (7); RNA `untestable` means not detected
(4,232). The Catalog holds 1,853 of the 1,854 canonical RNA calls: ENSG00000278732
(Metazoa_SRP) is one of 114 tested Ensembl IDs without a Catalog record, 103 of
which share a symbol with another ID. The 708 interface-untested genes with one
supported domain (698 RNA and 10 genetics) demonstrate why missing eligibility
is not negative evidence."""),
    (PAP, True,
     """The scale of assay limitation is an important Resource result. CosMx adequately
observes only 12 of 117 programs; GeoMx and GSE287826 cannot test any complete
program under the accepted rules. In single-cell ATAC, 38,066 of 39,914 jointly
testable peaks are `indeterminate` and 1,844 are `source_dependent`. The current
Catalog has zero strict `tested_negative` calls. These are not failed
validations. They show why modality-vote counting converts absent coverage,
unresolved donor identity, and assay mismatch into false confidence.""",
     """Observability limits are an important Resource result, and they have three
different causes. Panel content limits CosMx, which adequately observes 12 of
117 programs. Unresolved processing and metadata limit GeoMx (no gene axis
loaded; direction and normalization unresolved) and GSE287826 (no gene axis
loaded; no AOI-to-donor key): each currently tests 0 of 117 programs, and their
observability is not yet determined. Donor number limits single-cell ATAC (5
normal and 9 MASH donors in GSE244832, 6 and 6 in GSE281367): 38,066 of 39,914
jointly testable peaks are `indeterminate`, and 1,844 are `source_dependent`
(q < 0.05 in one cohort only; 1,505 of them have the same sign in the other
cohort). Against these, the four Visium/ST datasets each observe 114 of 117
programs. The current Catalog defines no negative-call rule, so
`tested_negative` is not assessed. These are not failed validations. They show
why modality-vote counting converts absent coverage, unresolved donor identity,
and assay mismatch into false confidence."""),
    (RES, False,
     "| Catalog strict `tested_negative` calls | 0 | No current row satisfies every adequacy requirement for a negative claim |",
     "| Catalog `tested_negative` calls | not assessed | No negative-call rule was defined for either domain; the negative-rule fields are empty in all 56,466 evidence rows, so a count of zero would be fixed by construction |"),
    (RES, False,
     "| Strict `tested_negative` calls | 0 |",
     "| `tested_negative` calls | not assessed (no negative-call rule defined) |"),
    (DIS, False,
     """These states have
different scientific consequences and remain distinct in the downloadable
tables and source graph. In the corrected r8 Catalog, zero rows satisfy the
strict `tested_negative` rule in the corrected Catalog build; that absence is a
valid Resource result, not a software failure.""",
     """These states have
different scientific consequences, but the r8 Catalog tables do not yet keep all
of them apart: `not_applicable` is used in no row, `source_dependent` is true in
every row, and the 9,312 genes with no Tier-1/2 COLOC record carry the label
`underpowered_source` although no power assessment was made. The r8 Catalog
defines no negative-call rule, so `tested_negative` was not assessed. Among the
456 genes with an eligible SuSiE signal pair in a Tier-1/2 study, 20 have PP.H3
above 0.8 and best PP.H4 at most 0.5, the pattern expected when the trait and
expression signals have distinct causal variants; 18 of them have the genetics
call indeterminate, the same label given to genes with no evidence for either a
shared or a distinct causal variant."""),
    (MET, False,
     """check, or an assay that cannot observe the proposed mechanism remains
indeterminate or untestable.""",
     """check, or an assay that cannot observe the proposed mechanism remains
indeterminate or untestable. The r8 Gene Catalog defines no such rule for either
domain, so no gene can be `tested_negative` there."""),
    (MET, False,
     """Complete per-assay coverage is in the Supplement. In summary: CosMx adequately
observes 12 of 117 programs and GeoMx and GSE287826 test no complete program
under the accepted rules; in single-cell ATAC, 38,066 of 39,914 jointly testable
peaks are `indeterminate` and 1,844 are `source_dependent`; and the Gene Catalog
contains zero strict `tested_negative` calls. These are limits on what the
assembled assays can currently resolve, and they are the reason modality-vote
counting is not used anywhere in this work.""",
     """Complete per-assay coverage is in the Supplement. In summary: CosMx adequately
observes 12 of 117 programs because of its targeted panel; GeoMx (no gene axis
loaded; direction and normalization unresolved) and GSE287826 (no gene axis
loaded; no AOI-to-donor key) test no complete program in the current data
processing, so their observability is not yet determined; and the four
Visium/ST datasets each observe 114 of 117 programs. In single-cell ATAC, with
5/9 and 6/6 normal/MASH donors per cohort, 38,066 of 39,914 jointly testable
peaks are `indeterminate` and 1,844 are `source_dependent` (q < 0.05 in one
cohort only; 1,505 share sign with the other cohort). The Gene Catalog defines
no negative-call rule, so `tested_negative` is not assessed. These are limits on
what the assembled assays and data can currently resolve, and they are the
reason modality-vote counting is not used anywhere in this work."""),
    (STA, False,
     "programs with 6,633 memberships, zero tested-negative rows.",
     "programs with 6,633 memberships; no negative-call rule is defined, so tested-negative is not assessed."),
    (RDM, False,
     """Zero
strict `tested_negative` calls is a central result, not an empty-output defect.""",
     """No
negative-call rule is defined, so `tested_negative` is not assessed; the zero is not a result."""),
]

E["R07-04"] = [
    (RES, False,
     "| `indeterminate` | Available data do not resolve the question |",
     "| `indeterminate` | No call under the prespecified rule. For RNA this includes 11,301 genes with BH q < 0.05 whose absolute log2FC is at most 0.5: the direction is resolved but the change is below the effect floor |"),
]

E["R07-05"] = [
    (MET, False,
     "Genes absent from an assay universe, and Catalog identities excluded by annotation or assay qualification, are labeled `not_jointly_testable` rather than neither.",
     "Genes absent from an assay universe are labeled `not_jointly_testable` rather than neither. The Catalog keeps one record per adjudicated identity in the union of the bulk fit and the Tier-1/2 COLOC tables. Identities quarantined during adjudication (one symbol mapping to more than one Ensembl ID, or an ID absent from GENCODE v49), second Ensembl IDs of a symbol already kept under another ID, and program-member genes outside both assay universes have no record and are not counted in any class; RESULTS lists them. The Catalog therefore holds 1,853 of the 1,854 canonical RNA calls."),
    (MET, False,
     "- `not_jointly_testable`: absent from an assay universe or excluded from the Catalog joint universe by annotation or assay qualification.",
     "- `not_jointly_testable`: absent from an assay universe, or untestable in one domain; quarantined identities have no Catalog record and are not counted in any class."),
    (RES, False,
     """One of the other 28 quarantines has
unilateral RNA support, so quarantine is not equivalent to absence of evidence.""",
     """One of the other 28 quarantines has
unilateral RNA support, so quarantine is not equivalent to absence of evidence.
Over the bulk fit, 114 tested Ensembl IDs have no Catalog record: 109 are
multi-Ensembl quarantines and 5 are second IDs of a symbol kept under another
ID (two POLR1HASP, two ABCC6P2 and one CES1P1), with no quarantine record. 35
of the 114 have BH q < 0.05 and one, ENSG00000278732 (Metazoa_SRP, log2FC
+0.535), is a canonical call, so the Catalog holds 1,853 of the 1,854 canonical
calls. Of the 61 Tier-1/2 COLOC IDs without a record, 32 are multi-Ensembl
quarantines, 22 are absent from GENCODE v49, six are the genetic IDs of the six
genes described below, and one (ENSG00000223797) is a second ENTPD3-AS1 ID whose
symbol is kept under ENSG00000280739. The 76 program-member IDs without a
record (82 of 6,633 memberships) are outside the bulk fit and every COLOC
table."""),
]

E["R07-06"] = [
    (LEG, False,
     """**(C)** The four fixed worked entries show their
separate calls, interface class and proposed gene-level and molecular-object
routes: *GLP1R*, untested because RNA is untestable; *HKDC1*, concordant;
*MTARC1*, unresolved; and *THRB*, genetically anchored. Gene-level genetic
support does not establish a credible-set target link or a causal mechanism.
**(D)** Counts routed by the five source classes: untested 13,551, unresolved
13,127, established-state associated 1,129, genetically anchored 400 and
concordant 26.""",
     """**(C)** The four fixed worked entries show their
separate calls with native estimates, interface class and proposed gene-level
and molecular-object routes. *GLP1R* is untested: RNA not detected and no
Tier-1/2 SuSiE result (best ABF PP.H4 0.022, MVP NAFLD EAS). *HKDC1* is
supported in both domains: SuSiE PP.H4 0.988 (UKBB AST) and log2FC +1.07 (BH q
8.0e-7). *MTARC1* is unresolved: no Tier-1/2 SuSiE result (best ABF PP.H4
0.046, FinnGen NASH) and log2FC -0.15 (q 1.2e-4). *THRB* is genetically
anchored: SuSiE PP.H4 > 0.999 (UKBB GGT) and log2FC -0.20 (q 4.2e-10). The RNA
calls of *MTARC1* and *THRB* are indeterminate because the change is below the
0.5 log2FC floor, not because the direction is unresolved. RNA estimates come
from 774 genetic-identity clusters (138 control, 636 disease). Gene-level genetic
support does not establish a credible-set target link or a causal mechanism.
**(D)** Counts routed by the five source classes: untested 13,551 (4,232 with
no RNA call because the gene was not detected; 9,319 with no genetic call, of
which 9,312 have no Tier-1/2 COLOC record and 7 are SuSiE-untestable in one
Tier-1/2 study), unresolved 13,127, established-state associated 1,129,
genetically anchored 400 and supported in both domains 26 (Catalog class
`concordant`). That class requires a supported call in each domain and does not
compare the direction of the genetic and RNA effects."""),
    (MET, False,
     "- `convergent`: both SuSiE `PP.H4 > 0.5` and canonical DEG.",
     "- `convergent` (Catalog class `concordant`): both SuSiE `PP.H4 > 0.5` and canonical DEG. The direction of the genetic effect is not compared with the RNA direction; the Catalog genetics rows carry no direction. The Catalog names for `genetic_only`, `disease_state_only`, `neither` and `not_jointly_testable` are genetically anchored, established-state associated, unresolved and untested."),
    (PAP, False,
     "anchored, 26 concordant, 1,129 established-state associated, 13,127 unresolved",
     "anchored, 26 supported in both domains (class `concordant`; effect directions not compared), 1,129 established-state associated, 13,127 unresolved"),
    (STA, False,
     """Primary classes:
genetically anchored 400, concordant 26, established-state associated 1,129,""",
     """Primary classes:
genetically anchored 400, concordant (both supported; direction not compared) 26, established-state associated 1,129,"""),
    (STA, False,
     "(6,633 memberships); classes genetically anchored 400, concordant 26, established-state associated 1,129",
     "(6,633 memberships); classes genetically anchored 400, concordant (both supported; direction not compared) 26, established-state associated 1,129"),
    (STA, False,
     """400 genetically anchored, 26 concordant); the working 6B and 6D were on the""",
     """400 genetically anchored, 26 concordant, meaning both supported with direction
not compared); the working 6B and 6D were on the"""),
    (STA, False,
     "full classes 400 genetically anchored, 26 concordant, 1,129 state-associated",
     "full classes 400 genetically anchored, 26 concordant (both supported; direction not compared), 1,129 state-associated"),
    (RES, False,
     "| Concordant | 26 |",
     "| Concordant (both domains supported; direction not compared) | 26 |"),
    (RES, False,
     """genetically anchored 400,
concordant 26. The working 6B and 6D were still""",
     """genetically anchored 400,
concordant 26 (both supported; direction not compared). The working 6B and 6D were still"""),
]

E["R07-07"] = [
    (PAP, True,
     """- The MASLD Gene Catalog identifies what is supported, what is untested or
  source-dependent, and which experiment can discriminate the leading
  explanations.""",
     """- The MASLD Gene Catalog identifies what is supported, what is untested
  because one domain has no call (RNA not detected for 4,232 genes; no Tier-1/2
  COLOC record for 9,312; SuSiE-untestable in one Tier-1/2 study for 7), what
  was tested without a call, and which experiment can discriminate the leading
  explanations."""),
    (DIS, False,
     """Calls distinguish supported, discordant, indeterminate, untestable, not
applicable, source-dependent, and prespecified adequate-negative states.""",
     """The vocabulary distinguishes supported, discordant, indeterminate, untestable, not
applicable, source-dependent, and prespecified adequate-negative states; the r8
gene calls use three of them (supported, indeterminate and untestable)."""),
    (RES, False,
     "| `not_applicable` | The assay does not address that biological question |",
     "| `not_applicable` | The assay does not address that biological question; not used in the r8 gene calls (0 of 56,466 rows) |"),
    (RES, False,
     "| `source_dependent` | Interpretation depends materially on one source or mapping assumption |",
     "| `source_dependent` | Interpretation depends materially on one source or mapping assumption; the r8 field is true in all 56,466 rows and does not separate genes |"),
]

E["R15-02"] = [
    (RES, False,
     """(3) 555 indeterminate genetics rows carry a PP.H4 above 0.5 without a supported
SuSiE call, which is the ABF value where SuSiE is absent or at most 0.5.""",
     """(3) the r8 PP.H4 column shows the best SuSiE PP.H4 for the 436 supported genes
and the best ABF PP.H4 for 15,231 of the 18,478 indeterminate genes, with no
field naming the method. 555 indeterminate rows show a value above 0.5 (152 at
or above 0.8); 222 of these have a finite SuSiE PP.H4 of at most 0.5 that the
ABF value replaces (GNMT shows ABF 0.706; its SuSiE value is 0.128). The
genetics call uses SuSiE only. A candidate per-gene table, not adopted, carries
the best Tier-1/2 SuSiE and ABF PP.H4 in separate columns, each with its GWAS
study and that study's ancestry, and for RNA the log2FC, standard error, t, p
and BH q from the 774-unit fit (138 control and 636 disease genetic-identity
clusters)."""),
]

E["R15-03"] = [
    (MET, False,
     """Each gene or program retains a stable identifier, evidence
class, testing universe, assay-specific call, biological unit, source dataset,
source-publication relationship, and parent-file graph.""",
     """Each gene or program retains a stable identifier, evidence
class, testing universe, assay-specific call, source dataset,
source-publication relationship, and parent-file graph. The r8 evidence rows
carry the native estimate and, for RNA, the BH q, but no standard error,
interval, p value or unit count."""),
]

E["R15-05"] = [
    (RES, False,
     "contains 28,233 records from the same r8 classes.",
     """contains 28,233 records from the same r8 classes. Its assay_applicability field
reads `indeterminate` for bulk RNA in all 1,853 RNA-supported records and for
chromatin and spatial context in every record; proteomics reads `indeterminate`
in 21,305 records and `not_applicable` in 6,928, and no record carries a protein
measurement. A candidate correction, not adopted, sets bulk RNA to the gene's
RNA call, proteomics to quantified or not quantified in PXD051911 (6,042 genes
quantified: 5,952 by direct symbol match and 90 after symbol remapping), and
chromatin and spatial to not joined at gene level. Protein significance is
reported in Supplementary Table 1 and the Figure 5B legend, not in the Catalog."""),
    (LEG, False,
     "Protein, chromatin and spatial context retain their assay or program unit.",
     """Chromatin and spatial context are attached for two programs (Stromal ECM and
Ductular injury); liver protein was untestable for both. No gene-level protein,
chromatin or spatial call is joined to the r8 Catalog, and no record is a
regulatory DNA element object."""),
    (ABS, True,
     """The MASLD Gene Catalog connects protein-coding genes, lncRNA
transcripts, regulatory DNA elements, and unresolved loci to their genetic,
disease-state, cell, and tissue contexts and to the next discriminating
experiment.""",
     """The MASLD Gene Catalog records separate genetic and disease-state RNA
evidence states for 28,233 gene identities, including protein-coding genes,
lncRNA transcripts and unresolved loci, and assigns each the next
discriminating experiment. It links its 3,422 program-member genes to the 117
cell programs; chromatin and spatial context are attached for two programs."""),
    (INT, False,
     """The resulting MASLD Gene Catalog links
protein-coding genes, lncRNA transcripts, regulatory DNA elements, and
unresolved loci to inherited variation, multicellular remodeling, cell and
tissue context, evidentiary boundaries, and the next discriminating
experiment.""",
     """The resulting MASLD Gene Catalog links
protein-coding genes, lncRNA transcripts and unresolved loci to inherited
variation and disease-state RNA evidence, links its 3,422 program-member genes
to the 117 cell programs, records evidentiary boundaries, and assigns the next
discriminating experiment; chromatin and spatial context are attached for two
programs."""),
]

E["R15-09"] = [
    (RES, False,
     """is nevertheless `source_dependent` because the ten arrays cannot be joined to
biological donors.""",
     """is nevertheless `source_dependent` because the ten arrays cannot be joined to
biological donors. The r8 Catalog program-context row for this result has call
state `supported` with the source-dependence flag set; its call state should
read `source_dependent`."""),
    (RES, False,
     """It remains separate from r8's fourteen historical
multimodal program contexts.""",
     """It remains separate from r8's fourteen historical
multimodal program contexts, which are the only program rows in the Catalog:
two programs (Stromal ECM and Ductular injury) in seven datasets from the
2026-08-07 candidate. No 117-program stage result or HMSMA spatial row is in
the Catalog."""),
]

E["R10-03"] = [
    (RES, False,
     """The 28 historical Myojin class-assay rows still test earlier classes and are
retained as historical, nonconfirmatory records.""",
     """The 28 historical Myojin class-assay rows still test earlier classes and are
retained as historical, nonconfirmatory records. Rerun with the r8 classes and
the frozen Plan 40 code on the 7,678 screen genes in the four prespecified
classes (84 untested genes excluded; OLS with covariates log1p HLF TPM, minimum
guide count and HLF Chronos), the prespecified omnibus partial F for the three
class indicators is 0.395, Freedman-Lane permutation p = 0.761 (100,000 draws
within 25 guide-count and expression strata). Standardized effects against
unresolved genes are -0.044 for genetically anchored (n = 207), 0.035 for
established-state associated (n = 348) and 0.181 for concordant (n = 11; HC3
95% CI -0.218 to 0.365 LFC). Counting the 84 untested genes in the reference
gives F 0.399, p = 0.760. No contrast reaches the prespecified standardized
effect of 0.20, so the nonconfirmatory verdict is unchanged. The unit is the
gene in a pooled CRISPR screen in one cell line, not a donor."""),
]

E["R07-11"] = [
    (RES, False,
     """Catalog-versus-native jointly eligible IDs differ by six Catalog-only unique
GENCODE-v49 identity resolutions, all unresolved, and 28 native-only IDs (21
quarantined and seven genetically untestable).""",
     """Catalog-versus-native jointly eligible IDs differ by six Catalog-only unique
GENCODE-v49 identity resolutions, all unresolved, and 28 native-only IDs (21
quarantined and seven genetically untestable). The six (CTSLP2, ITFG2-AS1,
MATR3, PAXX, POLR2J4 and SUGT1P3) are keyed on their bulk Ensembl ID, but
their genetics row is the COLOC record of a different ID that is absent from
GENCODE v49; the bulk IDs have no COLOC record. Under a same-ID rule they would
be untested (unresolved 13,121, untested 13,557)."""),
]

texts = {}
out = {}
bad = []
for fid, edits in E.items():
    out[fid] = []
    for f, claim, old, new in edits:
        if f not in texts:
            texts[f] = (ROOT / f).read_text()
        t = texts[f]
        n = t.count(old)
        if n != 1:
            bad.append((fid, f, n, old[:60]))
            line = "NOT_FOUND"
        else:
            start = t[: t.index(old)].count("\n") + 1
            end = start + old.count("\n")
            line = f"{start}" if start == end else f"{start}-{end}"
        for ch in ["—", "–"]:
            if ch in new:
                bad.append((fid, f, "dash", new[:60]))
        if re.search(r"\[(Inference|Unverified|Speculation|Estimate)\]", new):
            bad.append((fid, f, "tag", new[:60]))
        out[fid].append(dict(file=f, line=line, claim_level=claim, old_text=old, new_text=new))
# overlap check across all edits in the same file (spans must be disjoint)
spans = {}
for fid, eds in out.items():
    for e in eds:
        if e["line"] == "NOT_FOUND":
            continue
        t = texts[e["file"]]
        s = t.index(e["old_text"])
        spans.setdefault(e["file"], []).append((s, s + len(e["old_text"]), fid))
for f, sp in spans.items():
    sp.sort()
    for a, b in zip(sp, sp[1:]):
        if b[0] < a[1]:
            bad.append(("overlap", f, a[2], b[2]))
json.dump(out, open(OUT, "w"), indent=1)
print("problems:", bad)
for fid, eds in out.items():
    print(fid, [(e["file"].split("/")[-1], e["line"]) for e in eds])
