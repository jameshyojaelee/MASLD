> **HISTORICAL CAMPAIGN BRIEF — SUPERSEDED 2026-08-10.** Retained for
> provenance only. It cannot set the current thesis, figures, journals,
> required analyses, hero genes, convergence framing, or downstream-library
> membership. Use `docs/RESOURCE_PAPER_SCOPE.md` and the Plan 00/50/60 chain.

# Competitor-Positioning Brief — LOCKED FOR ITS 2026-07-13 CAMPAIGN ONLY

> **Purpose.** One source of truth for how the five figure-writing teams frame every external
> competitor. Paste or paraphrase the **LOCKED SENTENCE** blocks so no two figures contradict each
> other. This brief is **polish + consistency only** — it invents no new novelty claim and stays
> inside `docs/paper_outline.md` "Novelty Positioning" (L236–246) and "Major Limitations" (L248–261),
> which remain authoritative. If a figure needs to say more than this brief allows, escalate to the
> team lead; do not extend the framing locally.
>
> **Numbers authority:** `docs/manuscript/NUMBERS.md` (Broadaway comparison §2b, L216–228, L436–438).
> **Register authority:** `docs/competitor_paper_details.md` (#1–#16).

---

## 0. Master posture (applies to every comparator)

The framing is **complementary, never superiority**. Our contribution is the explicit population-scale
comparison of liver-trait regulatory genetics against the established disease-state transcriptome, and
the independent characterization of the resulting evidence classes. Competitors either (a) supply a
layer we do not generate (pQTL, cell-type chromatin QTL, MPRA function) that **independently reproduces
our separation** at a different modality, or (b) share our data/cohorts and are differentiated on
**method, not discovery**.

Three invariants that every team must hold:

1. **Corroboration, not convergence input.** External layers (Gobeil pQTL, Elison chromatin QTL,
   Zhu/Hu MPRA) are cited annotation joined onto the **frozen** evidence-class table. They are
   **never** inputs to the evidence-weighted heuristic (convergence) score. Do not revive
   "convergence-as-truth."
2. **Same symbol ≠ same signal.** A shared gene symbol is a descriptive cross-reference, not
   established shared causal support — especially where their tables are embargoed (Elison) or the
   colocalizing trait differs (enzyme vs direct-MASLD).
3. **Report the direct-MASLD value.** When we claim an external study corroborates one of our genetic
   hits, cite our **direct-MASLD** colocalization strength, not an off-trait / liver-enzyme-anchored
   maximum. Several shared-locus overlaps are enzyme-anchored on our side (see per-competitor notes).

Use the outline's required vocabulary: *independent corroboration*, *liver-trait regulatory genetics*,
*disease-state-associated*, *complementary*, *colocalized candidate effector*, *evidence-weighted
heuristic*. See the DO-NOT-SAY list (§8) for banned framings.

---

## 1. Broadaway 2024 — the reused liver-eQTL panel

**Register #(source panel); NUMBERS.md §2b.** Broadaway et al. 2024 (AJHG, PMID 39173627) is the
**source of our liver-eQTL panel (N = 1,183 APEX meta-analysis)** — the eQTL side is identical, so we
must **acknowledge reuse up front** (outline Limitation L250: "the liver eQTL source is not newly
generated and is largely European"). Their published colocalization (29 cardiometabolic GWAS,
`coloc.abf` PP.H4 ≥ 0.7, GCTA-cojo isolation) yields **747 effector eGenes**, lipid-dominated
(TC/LDL/nonHDL/logTG/HDL) with **no AST and no per-variant PIP**. The differentiation is **added
genetic content on the GWAS / fine-mapping side**, not a better eQTL panel.

> **LOCKED SENTENCE (differentiation).** "Our colocalization reuses the Broadaway et al. 2024 liver
> eQTL panel (N = 1,183); the added value is MASLD/liver-trait tiering, expanded and multi-ancestry
> GWAS integration, and comparison with diseased human tissue. Of our 751-gene colocalization set
> (PP.H4 > 0.5, either method), only 269 (36%) overlap their 747 lipid-dominated cardiometabolic
> eGenes; 482 are new to us, of which 328 (68%) are reachable only through data their design lacked
> (European AST, non-European ancestry, and dedicated NAFLD/PDFF/NASH GWAS). We additionally generate
> 5,172 SuSiE/CARMA credible-set variants with per-variant fine-mapping posterior probability (1,430
> PIP > 0.5; 952 PIP > 0.9; 198 loci, four ancestries) that their GCTA-cojo design never produced."

**Pre-empts:** the reviewer objection "the genetics is not novel — you reused Broadaway." Answer with
the **482 novel-to-us genes + 5,172 credible-set variants**, framed as expanded integration, not as
beating Broadaway. **Scope note:** the **751** count is the *Broadaway-comparison denominator*
(PP.H4 > 0.5, either method, full portfolio) — do **not** use it as the Fig 2 headline colocalization
count (headline = 473 SuSiE / 1,031 Tier-1/2 union; see NUMBERS.md). Keep 751 inside the Broadaway
paragraph only.

---

## 2. Elison / Gaulton 2025 — single-cell chromatin QTL (NEUTRAL prominence)

**Register #16; medRxiv 2025-05-11, DOI 10.1101/2025.05.09.25327043.** 86 deceased-donor livers,
five single-cell modalities (snRNA+snATAC multiome, Paired-Tag, scHi-C, Visium HD), per-cell-type
caQTL / eQTL / hQTL. Their **own** key result: caQTL/hQTL and eQTL at the same locus are **largely
distinct**, and most colocalizing loci act through **cell-type chromatin, not expression** — i.e. our
genetic-versus-disease-state separation **recurs one modality down**. This is the layer we lack
(single-cell chromatin / per-cell-type regulatory QTL); we cite it, we do not compete on it.

> **LOCKED SENTENCE (corroboration).** "An independent single-cell multiomic liver atlas (Elison
> et al. 2025) reports that regulatory-accessibility QTL and expression QTL colocalize with largely
> distinct risk loci and that most colocalizing loci act through cell-type chromatin rather than
> expression, reproducing the genetic-versus-disease-state separation one modality down; it
> corroborates our genetic candidates (RORA, PPP1R3B) and disease-state-associated genes (RORA, CUX2,
> RELB, EFHD1, KRT8) and resolves the cell of action for several of our genetic-only loci, while the
> canonical lipid-genetics loci absent from our expression colocalization are also absent from their
> chromatin atlas."

> **LOCKED SENTENCE (limitation-complement, for Discussion/Limitations).** "Our causal layer is
> tissue-level rather than cell-type-resolved and does not assign the cell of action or regulatory
> element at a risk locus; an independent single-cell multiomic atlas maps per-cell-type
> chromatin-accessibility QTL that supply that resolution (Elison et al. 2025), a complementary layer
> this resource does not generate."

**INTERNAL handling (do not print):** Bing Ren is institution-adjacent to our corresponding author
(NYGC) — **neutral prominence, no head-on chromatin competition, and do not present their orthogonality
quote as a claim about our results.** Their genome-wide tables are **embargoed** (as of 2026-07-11): use **only per-locus calls
stated in the preprint**, and state the cross-reference is a descriptive prototype (same-symbol
overlap does not establish same-signal support), extensible to full tables on release. **Caveat to
carry:** several shared-locus headline colocalizations (KRT8, EFHD1, TRPS1, FCGR2B) are
**liver-enzyme-anchored on our side** and weak on direct-MASLD — report the direct-MASLD value.

---

## 3. Gobeil / Arsenault 2026 — liver proteogenomics (pQTL) at the protein layer

**Register #15; medRxiv 2026-06-06, DOI 10.64898/2026.06.04.26354903.** Québec Obesity Biobank,
504 liver biopsies, DIA-MS proteomics + RNA-seq + genotyping in the **same** biopsies → in-house liver
**pQTL and eQTL**. Their headline — most differentially abundant liver proteins are **reactive
biomarkers, not causal drivers** — is our transcriptome thesis restated at the **protein layer**:
**independent corroboration of the biomarker-versus-driver separation.** They supply the liver-tissue
pQTL arm we do not generate (ours is expression-QTL only). Of 14 MR-significant proteins, **5
colocalized (PPH4 ≥ 0.80): ERLIN1, MTARC1, GCKR, LGALS1, MTTP.**

> **LOCKED SENTENCE (corroboration).** "An independent liver proteogenomic study (Gobeil et al. 2026,
> 504 biopsies) reproduces the biomarker-versus-driver separation at the protein layer and supplies a
> complementary protein-QTL causal map our expression-only arm does not generate; of their five
> colocalization-prioritized targets we independently flag four (ERLIN1, GCKR, LGALS1, MTTP) and miss
> exactly one, MTARC1 — a protein-altering effector (their protein-QTL PPH4 = 1.00 vs eQTL PPH4 =
> 0.003) that sits in the neither class of our expression map (our COLOC = 0.045) because its risk
> variant acts on protein, not transcript."

> **LOCKED SENTENCE (limitation-complement, for Discussion/Limitations).** "Our causal layer is
> expression-QTL only and does not resolve protein-altering or protein-level mechanisms (for example
> MTARC1); an independent liver proteogenomic study recovers such effectors through protein-QTL
> colocalization (Gobeil et al. 2026), a complementary causal layer this resource does not generate."

**Carry two caveats:** (1) **GCKR is direction-discordant** — their pQTL MR is protective, the eQTL
side (which our expression map captures) is risk-increasing; state that our colocalization captures
the expression side and the protein-level causal sign is opposite. (2) Their pQTL is **liver-tissue**
(the correct compartment) — the archived "plasma pQTL is the wrong compartment" defense does **not**
apply here. **Never** feed their calls into the convergence score; ingested as external cited
annotation only.

---

## 4. Zhu / Hu 2026 — MPRA + chromatin of MASLD risk variants (INDEPENDENT-MODALITY AGREEMENT)

**Nat Genet 2026;58(6):1280–1294, DOI 10.1038/s41588-026-02617-8** (Zhu, He, Xiao et al.; Guangzhou
National Lab / UPenn). Massively-parallel-reporter-assay + chromatin accessibility on MASLD GWAS risk
variants. **Not our snATAC source** (they reuse GSE244832 = Kim/Rosenthal 2025 and GSE189600 = Payen
2021 as public data, same accessions we use). Their MPRA-supported target genes cross-checked against
our COLOC gate and Cas13 library:

- **SLC22A3** — corroborates our genetic gate (COLOC ABF 0.83 / SuSiE 0.86); already Tier-1 and in the
  Cas13 library (coloc tier). *Note: direct-MASLD sub-gate 0.05–0.21; the coloc is enzyme/ALT-anchored
  — report the direct-MASLD value.*
- **LPL** — corroborated at the disease-state layer (DEG logFC +2.79, FDR 6.6e-32); already Tier-1 and
  in the library (core tier).
- **APOA5** — honest **gap**: MPRA-supported, our COLOC ~null (0.015–0.03), **not in the library**.
- **ANGPTL3** — **coincidental** overlap: in the library as a **positive control** (known
  lipid-lowering mechanism), not because it cleared our evidence gate.

> **LOCKED SENTENCE (agreement).** "A massively-parallel-reporter-assay study of MASLD risk variants
> (Zhu et al. 2026) independently supports two of our colocalization-gated candidates, SLC22A3 and LPL
> (both already Tier-1 and in our functional library) — independent-modality agreement between
> statistical colocalization and direct variant-level reporter function, not independent discovery.
> APOA5 is an honest gap (reporter-supported but null in our expression colocalization and absent from
> our library)."

**Framing rule:** this is **COLOC + MPRA agreeing**, **not** us discovering SLC22A3/LPL — do not claim
priority. Do **not** add APOA5 to the library on this basis (would require re-derived genetic support,
not an external MPRA hit).

---

## 5. Feng / MegaMASLD 2026 — differentiate on discordance, not another binary score

**Register #4 & #10; bioRxiv, DOI 10.64898/2026.02.19.706502.** 29 datasets, 2,640 samples (mixed
RNA-seq + microarray), 6-layer **binary** target score; EFHD1 ranked #1; MLIP validated in vitro.
**Reference obliquely** as "recent multi-omic meta-analyses" — **do not name confrontationally.** Our
differentiator is **not adding another binary score** but **testing and preserving evidence-class
discordance** (outline L243).

> **LOCKED SENTENCE (differentiation).** "Unlike recent multi-omic target scores that collapse
> evidence into a single binary or additive rank, this resource preserves discordance between evidence
> classes: a gene supported by liver-enzyme genetics is not equated with one supported by direct
> disease-state expression. For example, a top-ranked candidate from binary multi-omic scoring
> colocalizes in our atlas with ALT/GGT rather than with NAFLD or NASH (COLOC-Strong for liver-enzyme
> traits, not the disease phenotype) — an evidence-class distinction a binary score discards."

**Anchor fact (MLIP):** liver-enzyme-pleiotropic in our atlas (COLOC-Strong for ALT/GGT, not
NAFLD/NASH). **Anchor fact (EFHD1):** their #1; mediocre in our heuristic (rank ~top 22%). Use these
as illustrations of preserved discordance, not as attacks.

---

## 6. Open Targets / Nelson / Minikel — we do not rediscover genetic support

We do **not** re-derive the established result that genetic support improves drug-development success.
We use it to motivate why genetic support and tissue dysregulation are **not interchangeable** MASLD
target filters (outline L244).

> **LOCKED SENTENCE.** "That human genetic support improves drug-target success is established (Nelson
> et al.; Minikel et al.; Open Targets); this resource does not rediscover it. Instead it shows why
> genetic support and tissue dysregulation cannot be treated as interchangeable MASLD target filters —
> the two evidence classes are largely non-overlapping and identify mechanistically distinct
> candidates."

---

## 7. Figure → competitor routing map

| Figure | Competitors invoked | Role |
|---|---|---|
| **Fig 1** (atlas overview / framing) | Open Targets / Nelson / Minikel; (Feng obliquely) | Motivate the genetic-vs-disease-state framing; assert complementarity, not a leaderboard. |
| **Fig 2** (liver-trait colocalization) | **Broadaway** (reuse + 482/5,172 added content); **Gobeil** (boundary / MTARC1 panel E beside PNPLA3); **Zhu/Hu** (SLC22A3, LPL agreement) | Panel E: MTARC1 = pQTL-recovers / eQTL-misses exemplar. Broadaway overlap Venn in supplement (figS04). |
| **Fig 3** (disease-state stage-ordered remodeling) | none primary | Bulk/disease-state map; no head-to-head external comparator needed. |
| **Fig 4** (independent tissue modalities separate evidence classes) | **Elison** (chromatin QTL corroboration, neutral); modality corroboration (spatial/proteomic) | Elison as cited corroboration one modality down; figS2Q external panel. |
| **Fig 5** (mechanism-specific uses & limits) | **Open Targets / Nelson** (Discussion); **Gobeil** + **Elison** as complementary-layer limitations | THRB genetics-only; boundary framing; no convergence-superiority. |
| **Discussion / Limitations** | **Gobeil** (protein-QTL gap), **Elison** (cell-type-chromatin gap), **Zhu/Hu** (APOA5 gap), **Feng** (discordance) | Each named as a complementary layer this resource does not generate. |

---

## 8. DO-NOT-SAY list (all teams)

- ❌ **No priority/discovery language:** never "scooped," "first," "novel" as a priority claim,
  "we discovered [gene]" where an external study also reports it. Say **"independent corroboration"** /
  **"independent-modality agreement."**
- ❌ **No convergence-superiority revival.** Never present the convergence/heuristic score as a
  calibrated target probability or a universally superior class; the convergent primary class is small
  (outline L255). External pQTL/chromatin/MPRA layers are **annotation, never convergence inputs.**
- ❌ **No same-symbol = same-signal overclaim.** Shared gene symbols are descriptive cross-references;
  do not assert shared causal support, especially for Elison (embargoed tables) or enzyme-anchored
  overlaps.
- ❌ **No off-trait inflation.** Do not cite an enzyme/off-trait COLOC maximum as MASLD corroboration;
  report the **direct-MASLD** value (SLC22A3, KRT8, EFHD1, TRPS1, FCGR2B are enzyme-anchored on our
  side).
- ❌ **Elison: no head-on chromatin competition; neutral prominence** (institution-adjacent to our
  corresponding author). Do not present their orthogonality quote as evidence about our results.
- ❌ **Feng: no confrontational naming** — "recent multi-omic meta-analyses."
- ❌ **Do not claim we generate** pQTL, single-cell chromatin QTL, or MPRA data — we cite others' calls.
- ❌ **Outline-banned terms** (paper_outline.md L223–234): *causal gene based on COLOC alone; acquired
  risk; tissue-real; validated target (for observational evidence); therapeutic window; multicellular
  cascade (for cross-sectional data); prospective (absent a preregistered prediction); full
  drug-development gradient; recovers both approved targets; convergence is required.*

---

## 9. Locked numbers (cite these exact figures; source = NUMBERS.md)

| Fact | Value | Source |
|---|---|---|
| Broadaway liver eQTL panel (reused) | N = 1,183 (APEX meta) | NUMBERS.md §2b L218 |
| Broadaway effector eGenes | 747 (lipid-dominated, EUR, no PIP, no AST) | NUMBERS.md L218 |
| Our COLOC set for the comparison | 751 genes (PP.H4 > 0.5, either method) | NUMBERS.md L220 |
| Overlap with Broadaway | 269 / 751 = **36%** | NUMBERS.md L224, L436 |
| Novel to us | **482** (478 theirs-only) | NUMBERS.md L224 |
| Novel genes needing data Broadaway lacked | **328 / 482 = 68%** (AST 141 / non-EUR 120 / MASLD GWAS 93) | NUMBERS.md L228, L437 |
| Credible-set variants (PIP) we add | **5,172** (1,430 PIP>0.5; 952 PIP>0.9; 198 loci, 4 ancestries) | NUMBERS.md L228, L438 |
| Gobeil biopsies / proteins | 504 biopsies; DIA-MS 2,744 proteins | Register #15 |
| Gobeil coloc-prioritized set | 5 (ERLIN1, MTARC1, GCKR, LGALS1, MTTP); we flag 4/5 | Register #15 |
| MTARC1 (the miss) | their pQTL PPH4 = 1.00 / eQTL 0.003; our COLOC 0.045 | Register #15 |
| Elison donors / modalities | 86 livers; snRNA+snATAC / Paired-Tag / scHi-C / Visium HD | Register #16 |
| Elison corroborated (genetic / disease-state) | RORA, PPP1R3B / RORA, CUX2, RELB, EFHD1, KRT8 | Register #16 |
| Zhu/Hu agreement / gap | SLC22A3, LPL (agree, in library) / APOA5 (gap, not in library) | Zhu memo 2026-07-12 |
| Feng score / MLIP | 6-layer binary, 39 genes; MLIP = ALT/GGT COLOC not NAFLD/NASH | Register #4/#10 |

> **Reconciliation note.** The 5,172 credible-set count is the current reconciled NUMBERS.md figure
> (§2b, L228/L438) with full provenance; it supersedes the 2026-06-13 megareview "unverified" flag.
> If a team needs a variant count for a caption, pull it from NUMBERS.md at write time rather than
> hard-coding a differently-scoped number (e.g. the 6,912 GWAS-ATAC credible set is a *different*
> denominator).
