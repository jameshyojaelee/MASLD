"use client";

import { useEffect, useMemo, useState, useCallback } from "react";
import Link from "next/link";
import { Button } from "@/components/ui/button";
import { Badge } from "@/components/ui/badge";

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

interface ProteinRow {
  gene: string;
  logfc: number | null;
  padj: number | null;
  tstat: number | null;
  dataset: string;
}

interface StratumRow {
  stratum: string;
  n: number;
  rho: number | null;
  direction_pct: number | null;
}

interface ProteomicsSummary {
  n_proteins_tested: number;
  n_proteins_significant: number;
  n_atlas_tested: number;
  n_atlas_significant: number;
  overall_rho: number | null;
  conserved_core_rho: number | null;
  n_validated_targets: number;
  concordance_strata: StratumRow[];
  top_upregulated: ProteinRow[];
  top_downregulated: ProteinRow[];
  datasets: string[];
  provenance: Record<string, string>;
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

function fmtLFC(v: number | null): string {
  if (v == null) return "\u2014";
  return v >= 0 ? `+${v.toFixed(2)}` : v.toFixed(2);
}

function fmtPadj(v: number | null): string {
  if (v == null) return "\u2014";
  if (v === 0 || v < 1e-300) return "<1e-300";
  return v.toExponential(1);
}

function lfcColor(v: number | null): string {
  if (v == null) return "text-muted-foreground";
  if (v > 0) return "text-red-500 dark:text-red-400";
  if (v < 0) return "text-blue-500 dark:text-blue-400";
  return "text-muted-foreground";
}

type TableTab = "up" | "down";
type SortKey = "gene" | "logfc" | "padj";
type SortDir = "asc" | "desc";

// ---------------------------------------------------------------------------
// Stat card
// ---------------------------------------------------------------------------

function StatCard({
  label,
  value,
  accent,
  small,
}: {
  label: string;
  value: string;
  accent?: "green" | "blue" | "teal";
  small?: boolean;
}) {
  const accentClass =
    accent === "green"
      ? "text-green-600 dark:text-green-400"
      : accent === "blue"
        ? "text-blue-600 dark:text-blue-400"
        : accent === "teal"
          ? "text-teal-600 dark:text-teal-400"
          : "text-foreground";
  return (
    <div className="rounded-lg border border-border bg-card p-4">
      <p className="text-xs font-medium text-muted-foreground">{label}</p>
      <p
        className={`mt-1 font-semibold ${accentClass} ${small ? "text-sm" : "text-2xl"}`}
      >
        {value}
      </p>
    </div>
  );
}

// ---------------------------------------------------------------------------
// Page
// ---------------------------------------------------------------------------

export default function ProteomicsPage() {
  const [data, setData] = useState<ProteomicsSummary | null>(null);
  const [loading, setLoading] = useState(true);
  const [tab, setTab] = useState<TableTab>("up");
  const [sortKey, setSortKey] = useState<SortKey>("padj");
  const [sortDir, setSortDir] = useState<SortDir>("asc");

  useEffect(() => {
    fetch("/data/proteomics_summary.json")
      .then((r) => r.json())
      .then((d: ProteomicsSummary) => {
        setData(d);
        setLoading(false);
      })
      .catch(() => setLoading(false));
  }, []);

  const rows = useMemo(() => {
    if (!data) return [];
    const src = tab === "up" ? data.top_upregulated : data.top_downregulated;
    return [...src].sort((a, b) => {
      let av: string | number;
      let bv: string | number;
      if (sortKey === "gene") {
        av = a.gene;
        bv = b.gene;
      } else if (sortKey === "logfc") {
        av = a.logfc ?? 0;
        bv = b.logfc ?? 0;
      } else {
        av = a.padj ?? 1;
        bv = b.padj ?? 1;
      }
      let cmp: number;
      if (typeof av === "string" && typeof bv === "string") {
        cmp = av.localeCompare(bv);
      } else {
        cmp = (av as number) - (bv as number);
      }
      return sortDir === "asc" ? cmp : -cmp;
    });
  }, [data, tab, sortKey, sortDir]);

  const handleSort = useCallback(
    (key: SortKey) => {
      if (key === sortKey) {
        setSortDir((d) => (d === "asc" ? "desc" : "asc"));
      } else {
        setSortKey(key);
        setSortDir(key === "gene" ? "asc" : key === "padj" ? "asc" : "desc");
      }
    },
    [sortKey]
  );

  const sortIndicator = (key: SortKey) =>
    sortKey === key ? (sortDir === "asc" ? " \u25B2" : " \u25BC") : null;

  if (loading) {
    return (
      <div className="mx-auto max-w-6xl px-6 py-10">
        <h1 className="text-3xl font-bold tracking-tight">
          Plasma + Liver Proteomics
        </h1>
        <p className="mt-6 text-sm text-muted-foreground">
          Loading proteomics data...
        </p>
      </div>
    );
  }

  if (!data) {
    return (
      <div className="mx-auto max-w-6xl px-6 py-10">
        <h1 className="text-3xl font-bold tracking-tight">
          Plasma + Liver Proteomics
        </h1>
        <p className="mt-6 text-sm text-destructive">
          Failed to load proteomics data.
        </p>
      </div>
    );
  }

  return (
    <div className="mx-auto max-w-6xl px-6 py-8">
      {/* Header */}
      <div className="mb-8">
        <h1 className="text-3xl font-bold tracking-tight">
          Plasma + Liver Proteomics
        </h1>
        <p className="mt-2 text-muted-foreground">
          Protein-level evidence integrating Olink Explore 3072 plasma and
          DIA-MS liver proteomics, with mRNA-protein concordance across the
          MASLD atlas.
        </p>
      </div>

      {/* Summary stats */}
      <div className="mb-10 grid grid-cols-2 gap-4 sm:grid-cols-3 lg:grid-cols-5">
        <StatCard
          label="Proteins Tested"
          value={data.n_proteins_tested.toLocaleString()}
        />
        <StatCard
          label="Significant (atlas)"
          value={data.n_atlas_significant.toLocaleString()}
          accent="teal"
        />
        <StatCard
          label="Overall mRNA-Protein \u03C1"
          value={data.overall_rho != null ? data.overall_rho.toFixed(3) : "\u2014"}
          accent="blue"
        />
        <StatCard
          label="Conserved Core \u03C1"
          value={
            data.conserved_core_rho != null
              ? data.conserved_core_rho.toFixed(3)
              : "\u2014"
          }
          accent="green"
        />
        <StatCard
          label="Validated Drug Targets"
          value={data.n_validated_targets.toLocaleString()}
        />
      </div>

      {/* Concordance strata table */}
      <div className="mb-10">
        <h2 className="mb-4 text-xl font-semibold">mRNA-Protein Concordance</h2>
        <p className="mb-3 text-sm text-muted-foreground">
          Spearman correlation between dream logFC and protein logFC, stratified
          by gene class. Both-significant genes reach &rho; = 0.77 with 98.9%
          direction concordance, confirming protein-level support of the
          transcriptomic atlas.
        </p>
        <div className="overflow-x-auto rounded-lg border border-border">
          <table className="w-full text-sm">
            <thead>
              <tr className="border-b border-border bg-muted/50">
                <th className="px-3 py-2 text-left text-xs font-medium text-muted-foreground">
                  Stratum
                </th>
                <th className="px-3 py-2 text-right text-xs font-medium text-muted-foreground">
                  N
                </th>
                <th className="px-3 py-2 text-right text-xs font-medium text-muted-foreground">
                  Spearman &rho;
                </th>
                <th className="px-3 py-2 text-right text-xs font-medium text-muted-foreground">
                  Direction %
                </th>
              </tr>
            </thead>
            <tbody>
              {data.concordance_strata.map((r) => (
                <tr
                  key={r.stratum}
                  className="border-b border-border/50 transition-colors hover:bg-muted/30"
                >
                  <td className="px-3 py-1.5 text-xs font-medium">
                    {r.stratum}
                  </td>
                  <td className="px-3 py-1.5 text-right font-mono text-xs">
                    {r.n.toLocaleString()}
                  </td>
                  <td className="px-3 py-1.5 text-right font-mono text-xs">
                    {r.rho != null ? r.rho.toFixed(3) : "\u2014"}
                  </td>
                  <td className="px-3 py-1.5 text-right font-mono text-xs">
                    {r.direction_pct != null
                      ? r.direction_pct.toFixed(1) + "%"
                      : "\u2014"}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </div>

      {/* Top up/down */}
      <div className="mb-10">
        <div className="mb-4 flex items-center justify-between">
          <h2 className="text-xl font-semibold">Top Differential Proteins</h2>
          <div className="flex gap-2">
            <Button
              variant={tab === "up" ? "default" : "outline"}
              size="sm"
              onClick={() => setTab("up")}
            >
              Top Upregulated
            </Button>
            <Button
              variant={tab === "down" ? "default" : "outline"}
              size="sm"
              onClick={() => setTab("down")}
            >
              Top Downregulated
            </Button>
          </div>
        </div>

        <div className="overflow-x-auto rounded-lg border border-border">
          <table className="w-full text-sm">
            <thead>
              <tr className="border-b border-border bg-muted/50">
                <th
                  className="cursor-pointer select-none px-3 py-2 text-left text-xs font-medium text-muted-foreground hover:text-foreground"
                  onClick={() => handleSort("gene")}
                >
                  Gene{sortIndicator("gene")}
                </th>
                <th
                  className="cursor-pointer select-none px-3 py-2 text-right text-xs font-medium text-muted-foreground hover:text-foreground"
                  onClick={() => handleSort("logfc")}
                >
                  Protein logFC{sortIndicator("logfc")}
                </th>
                <th
                  className="cursor-pointer select-none px-3 py-2 text-right text-xs font-medium text-muted-foreground hover:text-foreground"
                  onClick={() => handleSort("padj")}
                >
                  padj{sortIndicator("padj")}
                </th>
                <th className="px-3 py-2 text-left text-xs font-medium text-muted-foreground">
                  Dataset
                </th>
              </tr>
            </thead>
            <tbody>
              {rows.map((r) => (
                <tr
                  key={r.gene}
                  className="border-b border-border/50 transition-colors hover:bg-muted/30"
                >
                  <td className="px-3 py-1.5">
                    <Link
                      href={`/gene/${encodeURIComponent(r.gene)}`}
                      className="font-mono font-semibold text-primary hover:underline"
                    >
                      {r.gene}
                    </Link>
                  </td>
                  <td
                    className={`px-3 py-1.5 text-right font-mono text-xs ${lfcColor(r.logfc)}`}
                  >
                    {fmtLFC(r.logfc)}
                  </td>
                  <td className="px-3 py-1.5 text-right font-mono text-xs">
                    {fmtPadj(r.padj)}
                  </td>
                  <td className="px-3 py-1.5">
                    <Badge variant="outline" className="text-[10px]">
                      {r.dataset}
                    </Badge>
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </div>

      {/* Provenance */}
      <div className="mb-10 rounded-lg border border-border bg-muted/20 p-5">
        <h3 className="mb-2 text-sm font-semibold">Data Provenance</h3>
        <ul className="space-y-1 text-xs text-muted-foreground">
          <li>
            <span className="font-mono">PXD052937</span> &mdash; DIA-MS liver
            proteomics, 72 samples
          </li>
          <li>
            <span className="font-mono">GSE276114</span> &mdash; liver fibrosis
            proteomics cohort
          </li>
          <li>Olink Explore 3072 plasma panel (Phase VIII plasma sweep)</li>
          <li>
            Datasets included:{" "}
            <span className="font-mono">{data.datasets.join(", ")}</span>
          </li>
        </ul>
      </div>
    </div>
  );
}
