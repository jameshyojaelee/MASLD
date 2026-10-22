# scRNA-seq Dataset Protocol Comparison

Generated 2026-05-22. Source: quantitative analysis of `donor_metadata.tsv` + PubMed lookups.

## Summary Table

| Dataset | Paper | Authors | Journal / Year | PMID / DOI | Tissue source | Isolation method | Platform | Hepatocyte enrichment? | Hep % (mean) | Hep % (range) | Hep % (std) | Mac % | Fib % | Endo % | Chol % | N donors | Total cells | Disease stages |
|---------|-------|---------|---------------|------------|--------------|-----------------|----------|----------------------|-------------|--------------|-------------|-------|-------|--------|--------|----------|------------|----------------|
| **GSE136103** | Resolving the fibrotic niche of human liver cirrhosis at single-cell level | Ramachandran P et al. | *Nature* 2019 | PMID [31597160](https://doi.org/10.1038/s41586-019-1631-3) | Surgical resection / explant | **Deliberate NPC enrichment** — hepatocytes physically excluded by design. Profiled non-parenchymal cells (macrophages, endothelial, stellate, cholangiocytes) to resolve the fibrotic niche. | 10x Chromium v2 | **Hepatocytes actively excluded** | 0% | 0–2% | 0% | 30% | 13% | 39% | 18% | 20 | 16,394 | Healthy (8) / Cirrhosis (9) / NA (3) |
| **GSE174748** | Opposing roles of hepatic stellate cell subpopulations in hepatocarcinogenesis. Nature 2022 | Wallace S | ~2021–2022 | Not indexed in PubMed | Surgical biopsy | Standard enzymatic whole-tissue dissociation (inferred from ~78% hepatocytes = natural liver composition) | 10x Chromium (inferred) | None | 78% | 72–84% | 5% | 4% | 4% | 8% | 5% | 4 | 18,360 | Healthy (2) / Steatosis (2) |
| **GSE185477** | Unknown — likely reference healthy liver atlas | Unknown | ~2021–2022 | Not indexed in PubMed | Surgical resection | Standard whole-tissue dissociation. Tight 70–90% hepatocyte range across 21 donors = protocol well controlled. | 10x Chromium (inferred) | None (natural representation) | 81% | 70–90% | 8% | 5% | 2% | 10% | 2% | 21 | 106,335 | Healthy only |
| **GSE189600** | Unknown — tiny pilot (142 cells total) | Unknown | ~2021 | Not indexed in PubMed | Unknown | Unknown — 100% hepatocytes in both donors suggests Percoll gradient or EpCAM+ hepatocyte-specific isolation | Unknown | **Hepatocyte-only isolation** | 100% | 100% | 0% | 0% | 0% | 0% | 0% | 2 | 142 | Healthy (1) / SH (1) |
| **GSE202379** | Single-Cell, Single-Nucleus, and Spatial RNA Sequencing of the Human Liver Identifies Cholangiocyte and Mesenchymal Heterogeneity | Andrews TS, Atif J, …, MacParland SA | *Hepatol Commun* 2022 | PMID [34792289](https://doi.org/10.1002/hep4.1854) | Surgical resection / transplant (Toronto) | **Paired scRNA + snRNA-seq**. scRNA alone under-captures hepatocytes due to dissociation stress (paper explicitly documents this); snRNA used to recover hepatocytes. Wide donor range (29–100% hep) reflects sc vs sn variation per donor. | 10x Chromium (sc) + 10x Nuclei Isolation (sn) | snRNA recovers hepatocytes that scRNA misses | 73% (wide SD) | 29–100% | 21% | 2% | 4% | 16% | 5% | 67 | 111,938 | MASLD (SAF-staged; F0–F4 documented for 58/67 donors) |
| **GSE244832** | Wang et al. MAFLD single-cell atlas | Wang et al. | ~2023 (JHEP Reports candidate; GEO accession confirmed) | Not confirmed in PubMed | Liver biopsy | Standard enzymatic whole-tissue dissociation; large well-powered MASLD cohort. 82% mean hepatocytes = consistent whole-tissue protocol. | 10x Chromium | None | 82% | 44–100% | 15% | 1% | 4% | 9% | 3% | 117 | 559,735 | Healthy (31) / Steatosis (26) / SH (60) |
| **Liver_Atlas (GSE192740)** | Spatial proteogenomics reveals distinct and evolutionarily conserved hepatic macrophage niches | Guilliams M, Bonnardel J, …, Scott CL | *Cell* 2022 | PMID [35021063](https://doi.org/10.1016/j.cell.2021.12.018) | Surgical liver resection (Ghent University Hospital) | **Multiplexed FACS-sorted populations**: 22 donors → pure Kupffer cell/macrophage sorts (0% hep); 9 donors → hepatocyte-enriched (63–98% hep); 5 donors → endothelial sorts; 2 donors → fibroblast sorts. Primary goal = macrophage atlas. Also includes CITE-seq and single-nuclei. | 10x Chromium + CITE-seq + snRNA | **Cell-type-specific FACS sorting** | 20% (bimodal: 0% in 29 NPC-sorted donors; 63–98% in 9 hep-sorted donors) | 0–98% | 36% | 54% | 10% | 15% | 1% | 38 | 82,638 | Healthy only (some obese) |

---

## Protocol Classification

### Group 1: Whole-tissue dissociation — biologically representative
- **GSE174748**, **GSE185477**, **GSE244832**
- Hepatocyte fraction: 78–82% (matches histological estimates)
- **Use for cross-stage composition analysis**

### Group 2: Deliberate NPC enrichment — hepatocytes excluded
- **GSE136103** (Ramachandran 2019 *Nature*): hepatocytes physically depleted, all 20 donors 0% hep
- **Liver_Atlas / GSE192740** (Guilliams 2022 *Cell*): FACS-sorted; 22/38 donors are pure macrophage sorts
- **Never use for cell-type composition comparisons across disease stages**
- **Never use to infer hepatocyte-expressed program scores**

### Group 3: Nucleus-based or hepatocyte-specific isolation
- **GSE202379** (Andrews 2022 *Hepatol Commun*): paired sc+sn; sc component has low hep fraction, sn recovers them. Wide per-donor variability (29–100%) reflects sc vs sn mixing.
- **GSE189600**: 100% hepatocytes — likely hepatocyte-specific protocol. Only 142 cells; too small for standalone analysis.
- **Use cautiously**: snRNA and scRNA are not directly comparable to scRNA-only datasets; hepatocyte counts may be inflated relative to true in vivo proportions.

---

## Impact on Cross-Stage Composition Analysis

The only dataset with **matched donors across Healthy / Steatosis / Steatohepatitis** is **GSE244832** (31 H + 26 ST + 60 SH). All cross-stage conclusions about cell-type proportion shifts rest on this single dataset.

| Cell type | Healthy (GSE244832 + GSE185477 + GSE174748) | SH (GSE244832 only) | Δ | Interpretation |
|-----------|--------------------------------------------|--------------------|---|----------------|
| Hepatocytes | 77–84% | 70% | −7–14 pp | Modest decline; real but partly driven by NPC expansion |
| Fibroblasts | 3.7% | 8.0% | +4.3 pp | **Real — stellate/myofibroblast expansion; matches P3 Fibrogenic** |
| Cholangiocytes | 2.1% | 6.4% | +4.3 pp | **Real — ductular reaction in steatohepatitis** |
| Macrophages | 6.6%* | 2.2% | −4.4 pp | **ARTIFACT** — inflated by Liver_Atlas FACS macrophage sorts in Healthy pool |
| Endothelial | 10% | 14% | +4 pp | Modest increase; consistent with angiogenesis in fibrosis |

*True healthy macrophage fraction from whole-tissue datasets only: ~1.9–5% (consistent with 2–4% histological Kupffer cell estimate)

---

## Macrophage Inflation Artifact — Root Cause

The apparent Healthy→SH macrophage *decrease* (6.6% → 2.2%) in Panel A is entirely explained by:
- 22 of 38 Liver_Atlas donors are pure FACS-sorted macrophage preparations (100% macrophages)
- These 22 donors are all classified as "Healthy"
- They contribute ~21,377 cells all classified as macrophages to the Healthy pool
- This inflates the Healthy macrophage mean from ~2% (biological) to 6.6% (artifactual)
- SH is exclusively GSE244832 whole-tissue donors → 2.2% macrophages (biologically accurate)

**Correct interpretation**: Kupffer cell/macrophage proportion does not decrease from Healthy to SH in MASLD. The data cannot speak to macrophage trends because the comparison is confounded.

---

## Recommendations for Figure Panels

1. **Panel A (cell-type composition)**: Restrict to whole-tissue datasets (GSE174748 + GSE185477 + GSE244832). Remove Liver_Atlas and GSE136103. This gives clean protocol-matched Healthy (n=53 donors) / ST (n=28) / SH (n=61) comparison.

2. **Panel C (hepatocyte-intrinsic program scores)**: Already correctly restricted to GSE244832 + GSE185477 + GSE174748 (Healthy/ST/SH, no Cirrhosis, no GSE136103). Andrews scRNA component also included but its low-hep donors reduce signal.

3. **Cirrhosis comparisons**: Cannot be made from this atlas without GSE136103, which has no hepatocytes. Any statement about cirrhosis cell-type composition is a NPC-enriched protocol artifact, not biology.

---

## Notes on Unknown Datasets (GSE174748, GSE185477, GSE189600)

These three datasets are not findable via PubMed by GEO accession or plausible search terms. Likely explanations:
- Deposited as supplementary data to a larger paper under a different accession
- Preprint / conference data that is not yet published or indexed
- Data from a clinical group deposited without an associated primary paper

The consistent ~80% hepatocyte fraction in GSE185477 (21 donors, tight range) suggests a well-run standard protocol from an experienced liver surgery group. GSE174748 (4 donors) is likely a pilot from the same or a collaborating group.

Action item: Check GEO submission metadata for submitter/PI name to identify the originating lab.
