/**
 * Network filter execution.
 *
 * Uses pre-built inverted indexes under
 * `/data/network/portal_export_v2/filter_indexes/*.json` when available
 * (built by the Python portal-export pipeline). When not deployed,
 * gracefully falls back to client-side scanning of the v2 search index
 * (which already carries `is_deg`, `is_coloc`, `is_drug`, and community
 * membership fields).
 *
 * The returned set is the intersection of all active filter predicates.
 */

import { loadSearchIndexV2 } from "./network-data";
import type { NetworkSearchEntryV2 } from "./network-types";
import { dataUrl } from "./data-base";

const BASE_V2 = dataUrl("network/portal_export_v2");

export type FStageFilter =
  | "all"
  | "F0-F1"
  | "F2"
  | "F3-F4"
  | "F2_emerging"
  | "F2_dissolving";

export interface NetworkFilters {
  fstage: FStageFilter;
  hasColoc: boolean;
  druggable: boolean;
  druggablePair: boolean;
  cohortSupport: number; // 0..10 minimum
  communityF01: number | null;
  communityF34: number | null;
}

export const DEFAULT_FILTERS: NetworkFilters = {
  fstage: "all",
  hasColoc: false,
  druggable: false,
  druggablePair: false,
  cohortSupport: 0,
  communityF01: null,
  communityF34: null,
};

export function filtersActive(f: NetworkFilters): boolean {
  return (
    f.fstage !== "all" ||
    f.hasColoc ||
    f.druggable ||
    f.druggablePair ||
    f.cohortSupport > 0 ||
    f.communityF01 !== null ||
    f.communityF34 !== null
  );
}

// ---------------------------------------------------------------------------
// Inverted-index cache (best-effort; absent indexes fall back to scan)
// ---------------------------------------------------------------------------

type IndexName =
  | "by_f_stage"
  | "by_coloc"
  | "by_druggable"
  | "by_druggable_pair"
  | "by_community_F01"
  | "by_community_F34"
  | "by_cohort_support";

const indexCache = new Map<IndexName, unknown>();
const indexMisses = new Set<IndexName>();

export async function loadFilterIndex<T = unknown>(
  name: IndexName
): Promise<T | null> {
  if (indexMisses.has(name)) return null;
  const cached = indexCache.get(name);
  if (cached !== undefined) return cached as T;
  try {
    const res = await fetch(`${BASE_V2}/filter_indexes/${name}.json`);
    if (!res.ok) {
      indexMisses.add(name);
      return null;
    }
    const data = (await res.json()) as T;
    indexCache.set(name, data);
    return data;
  } catch {
    indexMisses.add(name);
    return null;
  }
}

function intersect(a: Set<string> | null, b: Set<string>): Set<string> {
  if (a === null) return b;
  const out = new Set<string>();
  for (const x of a) if (b.has(x)) out.add(x);
  return out;
}

interface FStageIndexPayload {
  // { F2_emerging: ["GENE1", ...], F2_dissolving: [...], F0-F1: [...], ... }
  [stage: string]: string[];
}

interface NumericIndexPayload {
  // Either list-of-genes keyed by value, or map gene -> value
  [key: string]: string[] | number;
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/**
 * Executes the active filter set and returns matching gene symbols.
 * Uses inverted indexes when present, otherwise falls back to the v2
 * search index (which has is_deg/is_coloc/is_drug/community flags).
 *
 * For criteria that cannot be resolved from fallbacks (e.g. F2_emerging
 * or cohort_support without indexes), the predicate is silently skipped
 * and `unresolved` is populated so the UI can warn the user.
 */
export async function executeFilter(
  filters: NetworkFilters
): Promise<{ genes: string[]; unresolved: string[] }> {
  const unresolved: string[] = [];
  let acc: Set<string> | null = null;

  const searchIdx = await loadSearchIndexV2().catch(() => [] as NetworkSearchEntryV2[]);
  const allGenes = new Set<string>(searchIdx.map((g) => g.s));

  // --- COLOC ---
  if (filters.hasColoc) {
    const idx = await loadFilterIndex<FStageIndexPayload>("by_coloc");
    let set: Set<string>;
    if (idx && Array.isArray(idx["genes"])) {
      set = new Set(idx["genes"] as string[]);
    } else {
      // Fallback: search index is_coloc flag
      set = new Set(searchIdx.filter((g) => g.is_coloc).map((g) => g.s));
    }
    acc = intersect(acc, set);
  }

  // --- Druggable (single endpoint) ---
  if (filters.druggable) {
    const idx = await loadFilterIndex<FStageIndexPayload>("by_druggable");
    let set: Set<string>;
    if (idx && Array.isArray(idx["genes"])) {
      set = new Set(idx["genes"] as string[]);
    } else {
      set = new Set(searchIdx.filter((g) => g.is_drug).map((g) => g.s));
    }
    acc = intersect(acc, set);
  }

  // --- Druggable pair (both endpoints druggable on at least one edge) ---
  if (filters.druggablePair) {
    const idx = await loadFilterIndex<FStageIndexPayload>("by_druggable_pair");
    if (idx && Array.isArray(idx["genes"])) {
      acc = intersect(acc, new Set(idx["genes"] as string[]));
    } else {
      // No reliable fallback — skip but surface warning
      unresolved.push("druggable_pair");
    }
  }

  // --- F-stage ---
  if (filters.fstage !== "all") {
    const idx = await loadFilterIndex<FStageIndexPayload>("by_f_stage");
    if (idx && Array.isArray(idx[filters.fstage])) {
      acc = intersect(acc, new Set(idx[filters.fstage] as string[]));
    } else {
      unresolved.push(`f_stage=${filters.fstage}`);
    }
  }

  // --- Cohort support ---
  if (filters.cohortSupport > 0) {
    const idx = await loadFilterIndex<NumericIndexPayload>("by_cohort_support");
    if (idx) {
      // Expect either per-gene map { GENE: n } or buckets { "5": [...], "6": [...] }
      const matched = new Set<string>();
      // Detect shape
      const firstVal = Object.values(idx)[0];
      if (typeof firstVal === "number") {
        for (const [gene, n] of Object.entries(idx as Record<string, number>)) {
          if (n >= filters.cohortSupport) matched.add(gene);
        }
      } else {
        for (const [k, v] of Object.entries(idx)) {
          const n = parseInt(k, 10);
          if (!isNaN(n) && n >= filters.cohortSupport && Array.isArray(v)) {
            for (const g of v) matched.add(g);
          }
        }
      }
      acc = intersect(acc, matched);
    } else {
      unresolved.push(`cohort_support>=${filters.cohortSupport}`);
    }
  }

  // --- Community F01 / F34 (client-side resolvable via search index) ---
  if (filters.communityF01 !== null) {
    const target = filters.communityF01;
    const set = new Set(
      searchIdx.filter((g) => g.c_macro_F01 === target).map((g) => g.s)
    );
    acc = intersect(acc, set);
  }
  if (filters.communityF34 !== null) {
    const target = filters.communityF34;
    const set = new Set(
      searchIdx.filter((g) => g.c_macro_F34 === target).map((g) => g.s)
    );
    acc = intersect(acc, set);
  }

  // No active filter — return all genes from the search index
  const finalSet = acc ?? allGenes;
  return { genes: Array.from(finalSet).sort(), unresolved };
}

// ---------------------------------------------------------------------------
// Pre-computed query bundles (Q1/Q2/Q3 from pre-registration)
// ---------------------------------------------------------------------------

export interface GeneSummary {
  symbol: string;
  degree?: number;
  is_deg?: boolean;
  is_coloc?: boolean;
  is_drug?: boolean;
  c_macro_F01?: number;
  c_macro_F34?: number;
}

export async function executeQueryBundle(
  bundleId: string
): Promise<GeneSummary[] | null> {
  try {
    const res = await fetch(`${BASE_V2}/query_bundles/${bundleId}.json`);
    if (!res.ok) return null;
    const data = (await res.json()) as { genes: GeneSummary[] } | GeneSummary[];
    return Array.isArray(data) ? data : data.genes;
  } catch {
    return null;
  }
}

/** Hydrate a gene-symbol list into summary rows using the v2 search index. */
export async function hydrateGeneSummaries(
  symbols: string[]
): Promise<GeneSummary[]> {
  const idx = await loadSearchIndexV2().catch(() => [] as NetworkSearchEntryV2[]);
  const byName = new Map<string, NetworkSearchEntryV2>();
  for (const e of idx) if (!byName.has(e.s)) byName.set(e.s, e);
  return symbols.map((s) => {
    const e = byName.get(s);
    return {
      symbol: s,
      degree: e?.d,
      is_deg: e?.is_deg,
      is_coloc: e?.is_coloc,
      is_drug: e?.is_drug,
      c_macro_F01: e?.c_macro_F01,
      c_macro_F34: e?.c_macro_F34,
    };
  });
}

// ---------------------------------------------------------------------------
// URL param <-> filter state
// ---------------------------------------------------------------------------

export function filtersToParams(f: NetworkFilters): Record<string, string> {
  const out: Record<string, string> = {};
  if (f.fstage !== "all") out.f_stage = f.fstage;
  if (f.hasColoc) out.coloc = "1";
  if (f.druggable) out.druggable = "1";
  if (f.druggablePair) out.druggable_pair = "1";
  if (f.cohortSupport > 0) out.cohort_support = String(f.cohortSupport);
  if (f.communityF01 !== null) out.c_f01 = String(f.communityF01);
  if (f.communityF34 !== null) out.c_f34 = String(f.communityF34);
  return out;
}

export function filtersFromParams(params: URLSearchParams): NetworkFilters {
  const fstageRaw = params.get("f_stage");
  const allowed: FStageFilter[] = [
    "all",
    "F0-F1",
    "F2",
    "F3-F4",
    "F2_emerging",
    "F2_dissolving",
  ];
  const fstage = (allowed as string[]).includes(fstageRaw ?? "")
    ? (fstageRaw as FStageFilter)
    : "all";
  const cs = parseInt(params.get("cohort_support") ?? "0", 10);
  const cF01 = params.get("c_f01");
  const cF34 = params.get("c_f34");
  return {
    fstage,
    hasColoc: params.get("coloc") === "1",
    druggable: params.get("druggable") === "1",
    druggablePair: params.get("druggable_pair") === "1",
    cohortSupport: isNaN(cs) ? 0 : Math.max(0, Math.min(10, cs)),
    communityF01: cF01 !== null && cF01 !== "" ? parseInt(cF01, 10) : null,
    communityF34: cF34 !== null && cF34 !== "" ? parseInt(cF34, 10) : null,
  };
}
