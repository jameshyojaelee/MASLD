/**
 * Hero constellation data.
 *
 * Builds a small, decorative subgraph from the already-shipped
 * `public/data/knowledge_graph.json` (272 nodes / 200 edges): the top-N genes
 * by `evidence_score`, plus the pathway/drug "connector" nodes that link at
 * least two of those genes (so the constellation actually has edges — the KG
 * carries no direct gene–gene edges). Purely CHROME: nodes are sized by
 * evidence and colored from brand tokens, with NO diverging logFC scale and NO
 * quantitative legend. The result is memoized per `topN`.
 */

import { dataUrl } from "./data-base";

// ---------------------------------------------------------------------------
// Raw knowledge-graph shape (only the fields the hero reads)
// ---------------------------------------------------------------------------

interface RawKgNode {
  id: string;
  type: string; // "gene" | "pathway" | "drug"
  label: string;
  evidence_score?: number;
  composite_score?: number;
}

interface RawKgEdge {
  source: string;
  target: string;
  type?: string;
}

interface RawKg {
  nodes: RawKgNode[];
  edges: RawKgEdge[];
}

// ---------------------------------------------------------------------------
// Hero subgraph shape
// ---------------------------------------------------------------------------

export interface HeroNode {
  id: string;
  label: string;
  /** "gene" star (brand-colored, evidence-sized) vs faint "connector" hub. */
  kind: "gene" | "connector";
  /** Normalized size weight in [0, 1] (evidence for genes; small for hubs). */
  weight: number;
}

export interface HeroLink {
  source: string;
  target: string;
}

export interface HeroGraph {
  nodes: HeroNode[];
  links: HeroLink[];
}

// ---------------------------------------------------------------------------
// Memoized builder
// ---------------------------------------------------------------------------

const cache = new Map<number, HeroGraph>();
const inflight = new Map<number, Promise<HeroGraph>>();

function buildHeroGraph(kg: RawKg, topN: number): HeroGraph {
  const nodes = kg.nodes ?? [];
  const edges = kg.edges ?? [];
  const byId = new Map<string, RawKgNode>();
  for (const n of nodes) byId.set(n.id, n);

  const genes = nodes
    .filter((n) => n.type === "gene")
    .sort(
      (a, b) =>
        (b.evidence_score ?? 0) - (a.evidence_score ?? 0) ||
        (b.composite_score ?? 0) - (a.composite_score ?? 0)
    )
    .slice(0, topN);

  const geneIds = new Set(genes.map((g) => g.id));

  // Connector hubs = non-gene nodes adjacent to >= 2 of the selected genes.
  const hubHits = new Map<string, number>();
  for (const e of edges) {
    const s = e.source;
    const t = e.target;
    if (geneIds.has(s) && byId.get(t)?.type !== "gene") {
      hubHits.set(t, (hubHits.get(t) ?? 0) + 1);
    } else if (geneIds.has(t) && byId.get(s)?.type !== "gene") {
      hubHits.set(s, (hubHits.get(s) ?? 0) + 1);
    }
  }
  const hubIds = new Set(
    [...hubHits.entries()].filter(([, c]) => c >= 2).map(([id]) => id)
  );

  const keep = new Set<string>([...geneIds, ...hubIds]);

  // Normalize gene evidence to [0, 1] for sizing.
  const scores = genes.map((g) => g.evidence_score ?? 0);
  const min = scores.length ? Math.min(...scores) : 0;
  const max = scores.length ? Math.max(...scores) : 1;
  const span = max - min || 1;

  const heroNodes: HeroNode[] = [
    ...genes.map((g) => ({
      id: g.id,
      label: g.label,
      kind: "gene" as const,
      weight: ((g.evidence_score ?? 0) - min) / span,
    })),
    ...[...hubIds].map((id) => ({
      id,
      label: byId.get(id)?.label ?? id,
      kind: "connector" as const,
      weight: 0.12,
    })),
  ];

  const heroLinks: HeroLink[] = edges
    .filter((e) => keep.has(e.source) && keep.has(e.target))
    .map((e) => ({ source: e.source, target: e.target }));

  return { nodes: heroNodes, links: heroLinks };
}

/**
 * Fetch + build the memoized hero subgraph. Never throws to the caller for a
 * decorative element — a failed fetch resolves to an empty graph so the hero
 * simply keeps its static fallback.
 */
export function loadHeroGraph(topN = 40): Promise<HeroGraph> {
  const cached = cache.get(topN);
  if (cached) return Promise.resolve(cached);
  const pending = inflight.get(topN);
  if (pending) return pending;

  const p = fetch(dataUrl("knowledge_graph.json"))
    .then((r) => {
      if (!r.ok) throw new Error(`knowledge_graph.json ${r.status}`);
      return r.json() as Promise<RawKg>;
    })
    .then((kg) => {
      const graph = buildHeroGraph(kg, topN);
      cache.set(topN, graph);
      inflight.delete(topN);
      return graph;
    })
    .catch(() => {
      const empty: HeroGraph = { nodes: [], links: [] };
      inflight.delete(topN);
      return empty;
    });

  inflight.set(topN, p);
  return p;
}
