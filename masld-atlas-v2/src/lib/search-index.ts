import Fuse from "fuse.js";
import type { GeneIndexEntry } from "./types";

let fuseInstance: Fuse<GeneIndexEntry> | null = null;
let geneData: GeneIndexEntry[] = [];
let loadPromise: Promise<void> | null = null;

async function ensureLoaded(): Promise<void> {
  if (fuseInstance) return;
  if (loadPromise) return loadPromise;

  loadPromise = (async () => {
    const response = await fetch("/data/gene_index.json");
    geneData = await response.json();
    fuseInstance = new Fuse(geneData, {
      keys: [
        { name: "symbol", weight: 2 },
        { name: "ensembl_id", weight: 1 },
      ],
      threshold: 0.3,
      includeScore: true,
      minMatchCharLength: 2,
    });
  })();

  return loadPromise;
}

export async function searchGenes(
  query: string,
  limit = 20
): Promise<GeneIndexEntry[]> {
  if (!query || query.length < 2) return [];
  await ensureLoaded();
  if (!fuseInstance) return [];
  const results = fuseInstance.search(query, { limit });
  return results.map((r) => r.item);
}

export async function getGeneBySymbol(
  symbol: string
): Promise<GeneIndexEntry | null> {
  await ensureLoaded();
  return geneData.find((g) => g.symbol === symbol) ?? null;
}

export async function getGeneIndex(): Promise<GeneIndexEntry[]> {
  await ensureLoaded();
  return geneData;
}
