/**
 * Canonical headline numbers for the MASLD Gene Catalog portal.
 *
 * Every number that appears in portal copy MUST be imported from here so a
 * results refresh is a one-file edit. Sources (as of 2026-07-08):
 *   - DEG counts: canonical_deg_results.csv, effect-size-aware interval-null
 *     FDR gate (fdr < 0.05 at lfc = 0.25). 1,918 DEGs (1,419 up / 499 down).
 *   - GWAS portfolio: 50 GWAS across 5 ancestries.
 *   - COLOC effector genes: 473 SuSiE / 1,031 union (gene_level_coloc.csv).
 *   - Convergence: Tier-1 = 677 (convergence_evidence.csv).
 *   - Drug pipeline: 2 approved / 208 clinical / 1,504 preclinical.
 *   - Atlas dimensions: 27,187 genes.
 *   - Single-cell atlas: 1,232,318 cells.
 *
 * NOTE: the DEG gate is named ONLY by DEG_GATE_LABEL below. Never spell the
 * name of the underlying test function in copy (build guard enforces this).
 */

// --- Differential expression -------------------------------------------------
export const DEG_COUNT = 1918;
export const DEG_UP = 1419;
export const DEG_DOWN = 499;

/** LFC floor folded into the interval-null gate. */
export const DEG_LFC = 0.25;
/** FDR threshold of the interval-null gate. */
export const DEG_FDR = 0.05;

/**
 * The exact, approved phrasing for the DEG-calling gate. This is the ONLY
 * string permitted to describe the method in UI copy.
 */
export const DEG_GATE_LABEL =
  "effect-size-aware interval-null FDR gate (fdr < 0.05 at lfc = 0.25)";

// --- Genetics / COLOC --------------------------------------------------------
export const GWAS_COUNT = 50;
export const COLOC_SUSIE = 473;
export const COLOC_UNION = 1031;

// --- Convergence -------------------------------------------------------------
export const CONVERGENCE_TIER1 = 677;

// --- Drug pipeline -----------------------------------------------------------
export const DRUGS_APPROVED = 2;
export const DRUGS_CLINICAL = 208;
export const DRUGS_PRECLINICAL = 1504;

// --- Atlas / single-cell dimensions -----------------------------------------
export const ATLAS_GENES = 27187;
export const SC_CELLS = 1232318;

/** en-US thousands-separated formatting for display. */
export const fmt = (n: number): string => n.toLocaleString("en-US");
