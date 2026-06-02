/** Fisher's exact test engine for client-side gene set enrichment.
 *
 * Uses the hypergeometric distribution to compute over-representation
 * p-values, with Benjamini-Hochberg FDR correction.
 */

// ---------------------------------------------------------------------------
// Log-factorial with lazy cache (safe for universe sizes up to ~40 000)
// ---------------------------------------------------------------------------
const LOG_FACT_CACHE: number[] = [0];

function logFact(n: number): number {
  if (n <= 0) return 0;
  while (LOG_FACT_CACHE.length <= n) {
    LOG_FACT_CACHE.push(
      LOG_FACT_CACHE[LOG_FACT_CACHE.length - 1] + Math.log(LOG_FACT_CACHE.length)
    );
  }
  return LOG_FACT_CACHE[n];
}

function logHypergeometric(
  N: number,
  K: number,
  n: number,
  k: number
): number {
  return (
    logFact(K) +
    logFact(N - K) +
    logFact(n) +
    logFact(N - n) -
    logFact(k) -
    logFact(K - k) -
    logFact(n - k) -
    logFact(N - n - K + k) -
    logFact(N)
  );
}

// ---------------------------------------------------------------------------
// One-sided Fisher exact test (over-representation / enrichment)
// ---------------------------------------------------------------------------

/** Compute one-sided (upper tail) Fisher exact p-value.
 *
 * @param overlap  k – genes in both user list and pathway
 * @param listSize n – user list size
 * @param pathwaySize K – pathway gene set size
 * @param universeSize N – total genes in universe
 */
export function fisherExactOneSided(
  overlap: number,
  listSize: number,
  pathwaySize: number,
  universeSize: number
): number {
  let pval = 0;
  const maxK = Math.min(listSize, pathwaySize);
  for (let i = overlap; i <= maxK; i++) {
    pval += Math.exp(
      logHypergeometric(universeSize, pathwaySize, listSize, i)
    );
  }
  return Math.min(1, pval);
}

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

export interface EnrichmentResult {
  pathway_id: string;
  pathway_name: string;
  overlap: number;
  pathway_size: number;
  pvalue: number;
  padj: number;
  overlap_genes: string[];
}

export interface PathwayGeneSet {
  id: string;
  name: string;
  genes: string[];
  size: number;
}

// ---------------------------------------------------------------------------
// Main enrichment runner
// ---------------------------------------------------------------------------

/** Run Fisher enrichment on a user gene list against a collection of gene sets.
 *
 * Returns results sorted by ascending p-value with BH-adjusted p-values.
 */
export function runEnrichment(
  userGenes: string[],
  geneSets: PathwayGeneSet[],
  universeSize: number
): EnrichmentResult[] {
  const userSet = new Set(userGenes.map((g) => g.toUpperCase()));
  const n = userSet.size;

  const results: EnrichmentResult[] = [];

  for (const gs of geneSets) {
    const gsUpper = new Set(gs.genes.map((g) => g.toUpperCase()));
    const overlap_genes = [...userSet].filter((g) => gsUpper.has(g));
    const k = overlap_genes.length;
    if (k === 0) continue;

    const pval = fisherExactOneSided(k, n, gs.size, universeSize);
    results.push({
      pathway_id: gs.id,
      pathway_name: gs.name,
      overlap: k,
      pathway_size: gs.size,
      pvalue: pval,
      padj: 0, // filled after BH correction
      overlap_genes,
    });
  }

  // Benjamini-Hochberg FDR correction
  results.sort((a, b) => a.pvalue - b.pvalue);
  const m = results.length;
  for (let i = 0; i < m; i++) {
    results[i].padj = Math.min(1, (results[i].pvalue * m) / (i + 1));
  }
  // Enforce monotonicity (step-up)
  for (let i = m - 2; i >= 0; i--) {
    results[i].padj = Math.min(results[i].padj, results[i + 1].padj);
  }

  return results;
}
