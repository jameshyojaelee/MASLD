"use client";

import { useEffect, useState } from "react";
import { queryParquet } from "@/lib/duck";
import type { ModalityKey } from "@/lib/palette";

// ---------------------------------------------------------------------------
// Generic per-gene parquet query hook
// ---------------------------------------------------------------------------

export interface QueryResult<T> {
  data: T[] | null;
  error: string | null;
  loading: boolean;
}

/**
 * Register `file` (basename → view) and run `sql` against it, keyed by the SQL
 * string so it re-runs only when the gene (and thus the query) changes. Returns
 * `{ data, error, loading }`; `data` is null until the first result resolves.
 */
export function useParquetRows<T = Record<string, unknown>>(
  file: string,
  sql: string
): QueryResult<T> {
  const [data, setData] = useState<T[] | null>(null);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    let cancelled = false;
    setData(null);
    setError(null);
    queryParquet<T>(file, sql)
      .then((rows) => {
        if (!cancelled) setData(rows);
      })
      .catch((e: unknown) => {
        if (!cancelled) setError(e instanceof Error ? e.message : String(e));
      });
    return () => {
      cancelled = true;
    };
  }, [file, sql]);

  return { data, error, loading: data === null && error === null };
}

// ---------------------------------------------------------------------------
// atlas_core row (48-column UI subset) — fields used by the gene page
// ---------------------------------------------------------------------------

export interface AtlasCoreRow {
  human_symbol: string;
  ensembl_id: string | null;
  gene_biotype: string | null;
  mouse_ortholog: string | null;
  bulk_logFC: number | null;
  bulk_padj: number | null;
  treat_lfc: number | null;
  treat_fdr: number | null;
  is_deg: boolean | null;
  bulk_logFC_M: number | null;
  bulk_logFC_F: number | null;
  sex_class: string | null;
  sex_interaction_padj: number | null;
  coloc_best_susie_pp4: number | null;
  coloc_best_susie_gwas: string | null;
  coloc_abf_best_pp4: number | null;
  coloc_abf_best_gwas: string | null;
  n_coloc_sources: number | null;
  n_ancestry_gwas: number | null;
  coloc_cross_ancestry_replicated: boolean | null;
  coloc_susie_conf_tier: string | null;
  twas_z: number | null;
  twas_pval: number | null;
  convergence_rank: number | null;
  convergence_score: number | null;
  convergence_tier: string | null;
  concordance_state: string | null;
  essentiality_chronos: number | null;
  is_essential: boolean | null;
  spatial_is_svg: boolean | null;
  spatial_morans_i: number | null;
  spatial_consensus_direction: string | null;
  zonation_class: string | null;
  ferroptosis_class: string | null;
  dgidb_druggable: boolean | null;
  opentargets_drug: boolean | null;
  max_phase_masld: number | null;
  drug_dev_status: string | null;
  pharos_tdl: string | null;
  dominant_program_for_gene: string | null;
  dominant_program_logFC: number | null;
  layers_active: number | null;
  is_conserved: boolean | null;
}

// ---------------------------------------------------------------------------
// Numeric / label formatting
// ---------------------------------------------------------------------------

export const DASH = "—";

/** Coerce to a finite number or null. */
export function num(v: unknown): number | null {
  if (v == null) return null;
  const n = typeof v === "number" ? v : Number(v);
  return Number.isFinite(n) ? n : null;
}

export const clamp01 = (t: number): number => (t < 0 ? 0 : t > 1 ? 1 : t);

const mlog10 = (p: number | null): number =>
  p == null ? 0 : -Math.log10(Math.max(p, 1e-300));

/** Signed fold-change, 2 dp. */
export function fmtLogFC(v: number | null | undefined): string {
  const n = num(v);
  if (n == null) return DASH;
  return `${n >= 0 ? "+" : ""}${n.toFixed(2)}`;
}

/** p-value / FDR: fixed for moderate, scientific for small. */
export function fmtP(v: number | null | undefined): string {
  const n = num(v);
  if (n == null) return DASH;
  if (n === 0) return "0";
  if (n < 1e-3) return n.toExponential(1);
  return n.toFixed(3);
}

/** Posterior probability, 3 dp. */
export function fmtPP4(v: number | null | undefined): string {
  const n = num(v);
  return n == null ? DASH : n.toFixed(3);
}

/** ↑ / ↓ / — direction glyph from a signed value. */
export function directionArrow(v: number | null | undefined): string {
  const n = num(v);
  if (n == null || n === 0) return DASH;
  return n > 0 ? "↑" : "↓";
}

// ---------------------------------------------------------------------------
// s1–s8 evidence fingerprint derived from atlas_core columns
// ---------------------------------------------------------------------------

/**
 * Per-modality strengths in [0, 1] for the header EvidenceBar. Only the
 * modalities that atlas_core actually measures are populated (human bulk,
 * genetic, essentiality). Spatial context is categorical and dataset-qualified,
 * so it is deliberately excluded from this numeric evidence-strength object.
 */
export function evidenceStrengths(
  r: AtlasCoreRow
): Partial<Record<ModalityKey, number>> {
  const s: Partial<Record<ModalityKey, number>> = {};

  // S1 human bulk RNA-seq
  if (r.is_deg) s.s1_human = clamp01(Math.max(0.45, mlog10(r.treat_fdr) / 6));
  else if ((num(r.bulk_padj) ?? 1) < 0.05) s.s1_human = 0.3;
  else s.s1_human = 0;

  // S2 genetic (best colocalization PP.H4, already [0,1])
  const pp4 = Math.max(num(r.coloc_best_susie_pp4) ?? 0, num(r.coloc_abf_best_pp4) ?? 0);
  s.s2_genetic = clamp01(pp4);

  // S3 essentiality (DepMap CHRONOS; more negative = more essential)
  const chronos = num(r.essentiality_chronos);
  s.s3_essential = r.is_essential
    ? 0.85
    : chronos != null
      ? clamp01(-chronos) * 0.5
      : 0;

  return s;
}
