"use client";

import { useEffect, useMemo, useState, useCallback } from "react";
import Link from "next/link";
import { Button } from "@/components/ui/button";
import { Badge } from "@/components/ui/badge";
import { Input } from "@/components/ui/input";
import { Textarea } from "@/components/ui/textarea";
import { EvidenceFingerprint } from "@/components/evidence-fingerprint";
import { EnrichmentPanel } from "@/components/enrichment-panel";
import { getGeneIndex } from "@/lib/search-index";
import type { GeneIndexEntry } from "@/lib/types";

type Preset = "all" | "deg" | "conserved" | "druggable";
type SortKey =
  | "symbol"
  | "dream_logfc"
  | "dream_padj"
  | "sex_class"
  | "layers_active";
type SortDir = "asc" | "desc";

const PAGE_SIZE = 50;

const PRESET_OPTIONS: { key: Preset; label: string }[] = [
  { key: "all", label: "All" },
  { key: "deg", label: "DEGs Only" },
  { key: "conserved", label: "Conserved Core" },
  { key: "druggable", label: "Druggable" },
];

const COLUMNS: { key: SortKey; label: string; className?: string }[] = [
  { key: "symbol", label: "Symbol" },
  { key: "dream_logfc", label: "logFC", className: "text-right" },
  { key: "dream_padj", label: "padj", className: "text-right" },
  { key: "sex_class", label: "Sex Class" },
  { key: "layers_active", label: "Layers", className: "text-right" },
];

function formatPadj(val: number | null | undefined): string {
  if (val == null) return "\u2014";
  if (val < 1e-300) return "<1e-300";
  return val.toExponential(1);
}

function formatLogFC(val: number | null | undefined): string {
  if (val == null) return "\u2014";
  return val >= 0 ? `+${val.toFixed(2)}` : val.toFixed(2);
}

function logfcColor(val: number | null | undefined): string {
  if (val == null) return "text-muted-foreground";
  if (val > 0) return "text-red-500 dark:text-red-400";
  if (val < 0) return "text-blue-500 dark:text-blue-400";
  return "text-muted-foreground";
}

function getSortValue(
  gene: GeneIndexEntry,
  key: SortKey
): string | number | null {
  switch (key) {
    case "symbol":
      return gene.symbol;
    case "dream_logfc":
      return gene.dream_logfc ?? null;
    case "dream_padj":
      return gene.dream_padj ?? null;
    case "sex_class":
      return gene.sex_class ?? "";
    case "layers_active":
      return gene.layers_active ?? 0;
  }
}

function compareFn(
  a: GeneIndexEntry,
  b: GeneIndexEntry,
  key: SortKey,
  dir: SortDir
): number {
  const av = getSortValue(a, key);
  const bv = getSortValue(b, key);

  // Nulls always sort last
  if (av == null && bv == null) return 0;
  if (av == null) return 1;
  if (bv == null) return -1;

  let cmp: number;
  if (typeof av === "string" && typeof bv === "string") {
    cmp = av.localeCompare(bv);
  } else {
    cmp = (av as number) - (bv as number);
  }
  return dir === "asc" ? cmp : -cmp;
}

/** Parse a pasted gene list (comma, newline, tab, or space separated) */
function parseGeneList(text: string): string[] {
  return text
    .split(/[\s,;\t\n]+/)
    .map((s) => s.trim().toUpperCase())
    .filter((s) => s.length > 0);
}

export default function ExplorePage() {
  const [allGenes, setAllGenes] = useState<GeneIndexEntry[]>([]);
  const [loading, setLoading] = useState(true);
  const [searchQuery, setSearchQuery] = useState("");
  const [preset, setPreset] = useState<Preset>("all");
  const [sortKey, setSortKey] = useState<SortKey>("layers_active");
  const [sortDir, setSortDir] = useState<SortDir>("desc");
  const [page, setPage] = useState(0);

  // Enrichment & gene list state
  const [showEnrichment, setShowEnrichment] = useState(false);
  const [showGeneListInput, setShowGeneListInput] = useState(false);
  const [customGeneText, setCustomGeneText] = useState("");
  const [customGenes, setCustomGenes] = useState<string[]>([]);

  useEffect(() => {
    getGeneIndex().then((data) => {
      setAllGenes(data);
      setLoading(false);
    });
  }, []);

  const filtered = useMemo(() => {
    let result = allGenes;

    // Custom gene list filter (takes precedence when active)
    if (customGenes.length > 0) {
      const gset = new Set(customGenes.map((g) => g.toUpperCase()));
      result = result.filter((g) => gset.has(g.symbol.toUpperCase()));
    }

    if (searchQuery.length >= 2) {
      const q = searchQuery.toUpperCase();
      result = result.filter((g) => g.symbol.toUpperCase().includes(q));
    }

    switch (preset) {
      case "deg":
        result = result.filter((g) => g.is_deg);
        break;
      case "conserved":
        result = result.filter((g) => g.is_conserved_core);
        break;
      case "druggable":
        result = result.filter((g) => g.dgidb_druggable);
        break;
    }

    return result;
  }, [allGenes, searchQuery, preset, customGenes]);

  const sorted = useMemo(() => {
    return [...filtered].sort((a, b) => compareFn(a, b, sortKey, sortDir));
  }, [filtered, sortKey, sortDir]);

  // Gene symbols for enrichment (from custom list or current filter)
  const enrichmentGenes = useMemo(() => {
    if (customGenes.length > 0) return customGenes;
    return filtered.map((g) => g.symbol);
  }, [customGenes, filtered]);

  const totalPages = Math.max(1, Math.ceil(sorted.length / PAGE_SIZE));
  const safePage = Math.min(page, totalPages - 1);
  const pageData = sorted.slice(
    safePage * PAGE_SIZE,
    (safePage + 1) * PAGE_SIZE
  );

  // Reset page when filters change
  useEffect(() => {
    setPage(0);
  }, [searchQuery, preset, sortKey, sortDir, customGenes]);

  const handleSort = useCallback(
    (key: SortKey) => {
      if (sortKey === key) {
        setSortDir((d) => (d === "asc" ? "desc" : "asc"));
      } else {
        setSortKey(key);
        setSortDir(key === "symbol" ? "asc" : "desc");
      }
    },
    [sortKey]
  );

  const sortIndicator = (key: SortKey) => {
    if (sortKey !== key) return null;
    return sortDir === "asc" ? " \u25B2" : " \u25BC";
  };

  return (
    <div className="w-full px-6 py-8">
      {/* Header */}
      <div className="mb-6">
        <h1 className="text-3xl font-bold tracking-tight">Gene Explorer</h1>
        <p className="mt-2 text-muted-foreground">
          Search, filter, and compare genes across 7 evidence layers.
        </p>
      </div>

      {/* Filter bar */}
      <div className="mb-4 flex flex-col gap-3 sm:flex-row sm:items-center sm:gap-4">
        <Input
          placeholder="Search by symbol..."
          value={searchQuery}
          onChange={(e) => setSearchQuery(e.target.value)}
          className="max-w-xs"
        />
        <div className="flex flex-wrap gap-1.5">
          {PRESET_OPTIONS.map((p) => (
            <Button
              key={p.key}
              variant={preset === p.key ? "default" : "outline"}
              size="sm"
              onClick={() => setPreset(p.key)}
            >
              {p.label}
            </Button>
          ))}

          {/* Gene list paste toggle */}
          <Button
            variant={showGeneListInput ? "secondary" : "outline"}
            size="sm"
            onClick={() => setShowGeneListInput((v) => !v)}
          >
            Paste Gene List
          </Button>

          {/* Run enrichment */}
          <Button
            variant={showEnrichment ? "secondary" : "outline"}
            size="sm"
            disabled={loading || enrichmentGenes.length === 0}
            onClick={() => setShowEnrichment((v) => !v)}
          >
            Run Enrichment
            {enrichmentGenes.length > 0 && !loading && (
              <Badge variant="secondary" className="ml-1.5 text-[10px]">
                {enrichmentGenes.length.toLocaleString()}
              </Badge>
            )}
          </Button>
        </div>
      </div>

      {/* Gene list paste area (collapsible) */}
      {showGeneListInput && (
        <div className="mb-4 rounded-lg border border-border bg-muted/20 p-4">
          <label className="mb-1.5 block text-sm font-medium">
            Paste gene symbols
            <span className="ml-1 font-normal text-muted-foreground">
              (comma, space, or newline separated)
            </span>
          </label>
          <Textarea
            placeholder={"TP53, BRCA1, PNPLA3\nor one gene per line..."}
            value={customGeneText}
            onChange={(e) => setCustomGeneText(e.target.value)}
            className="mb-2 max-h-32 font-mono text-xs"
            rows={3}
          />
          <div className="flex items-center gap-2">
            <Button
              size="sm"
              onClick={() => {
                const parsed = parseGeneList(customGeneText);
                setCustomGenes(parsed);
                setPage(0);
              }}
            >
              Apply ({parseGeneList(customGeneText).length} genes)
            </Button>
            <Button
              variant="ghost"
              size="sm"
              onClick={() => {
                setCustomGeneText("");
                setCustomGenes([]);
                setPage(0);
              }}
            >
              Clear
            </Button>
            {customGenes.length > 0 && (
              <span className="text-xs text-muted-foreground">
                {customGenes.length} genes applied &middot;{" "}
                {filtered.length} found in atlas
              </span>
            )}
          </div>
        </div>
      )}

      {/* Results count */}
      <div className="mb-3 text-sm text-muted-foreground">
        {loading ? (
          "Loading gene index..."
        ) : (
          <>
            Showing{" "}
            <span className="font-medium text-foreground">
              {sorted.length.toLocaleString()}
            </span>{" "}
            of{" "}
            <span className="font-medium text-foreground">
              {allGenes.length.toLocaleString()}
            </span>{" "}
            genes
          </>
        )}
      </div>

      {/* Table */}
      <div className="overflow-x-auto rounded-lg border border-border">
        <table className="w-full text-sm">
          <thead>
            <tr className="border-b border-border bg-muted/50">
              <th className="px-3 py-2 text-left font-medium text-muted-foreground">
                {/* Evidence column - no sort */}
                <span className="text-xs">Evidence</span>
              </th>
              {COLUMNS.map((col) => (
                <th
                  key={col.key}
                  className={`cursor-pointer select-none px-3 py-2 font-medium text-muted-foreground hover:text-foreground ${col.className ?? "text-left"}`}
                  onClick={() => handleSort(col.key)}
                >
                  <span className="text-xs">
                    {col.label}
                    {sortIndicator(col.key)}
                  </span>
                </th>
              ))}
              <th className="px-3 py-2 text-left font-medium text-muted-foreground">
                <span className="text-xs">Tags</span>
              </th>
            </tr>
          </thead>
          <tbody>
            {loading ? (
              <tr>
                <td
                  colSpan={7}
                  className="px-3 py-12 text-center text-muted-foreground"
                >
                  Loading...
                </td>
              </tr>
            ) : pageData.length === 0 ? (
              <tr>
                <td
                  colSpan={7}
                  className="px-3 py-12 text-center text-muted-foreground"
                >
                  No genes match the current filters.
                </td>
              </tr>
            ) : (
              pageData.map((gene) => (
                <tr
                  key={gene.symbol}
                  className="border-b border-border/50 transition-colors hover:bg-muted/30"
                >
                  {/* Evidence fingerprint */}
                  <td className="px-3 py-1.5">
                    <EvidenceFingerprint
                      evidence={gene.evidence}
                      size={24}
                      showTooltip={true}
                    />
                  </td>

                  {/* Symbol */}
                  <td className="px-3 py-1.5">
                    <Link
                      href={`/gene/${encodeURIComponent(gene.symbol)}`}
                      className="font-mono font-semibold text-primary hover:underline"
                    >
                      {gene.symbol}
                    </Link>
                  </td>

                  {/* logFC */}
                  <td
                    className={`px-3 py-1.5 text-right font-mono text-xs ${logfcColor(gene.dream_logfc)}`}
                  >
                    {formatLogFC(gene.dream_logfc)}
                  </td>

                  {/* padj */}
                  <td className="px-3 py-1.5 text-right font-mono text-xs text-muted-foreground">
                    {formatPadj(gene.dream_padj)}
                  </td>

                  {/* Sex Class */}
                  <td className="px-3 py-1.5 text-xs text-muted-foreground">
                    {gene.sex_class ?? "\u2014"}
                  </td>

                  {/* Layers */}
                  <td className="px-3 py-1.5 text-right font-mono text-xs">
                    <span className="text-foreground">
                      {gene.layers_active ?? 0}
                    </span>
                    <span className="text-muted-foreground">/7</span>
                  </td>

                  {/* Tags */}
                  <td className="px-3 py-1.5">
                    <div className="flex flex-wrap gap-1">
                      {gene.is_deg && (
                        <Badge variant="default" className="text-[10px]">
                          DEG
                        </Badge>
                      )}
                      {gene.is_conserved_core && (
                        <Badge variant="secondary" className="text-[10px]">
                          CC
                        </Badge>
                      )}
                      {gene.dgidb_druggable && (
                        <Badge variant="outline" className="text-[10px]">
                          Drug
                        </Badge>
                      )}
                    </div>
                  </td>
                </tr>
              ))
            )}
          </tbody>
        </table>
      </div>

      {/* Pagination */}
      {!loading && sorted.length > PAGE_SIZE && (
        <div className="mt-4 flex items-center justify-between">
          <p className="text-xs text-muted-foreground">
            Page {safePage + 1} of {totalPages}
          </p>
          <div className="flex gap-2">
            <Button
              variant="outline"
              size="sm"
              disabled={safePage === 0}
              onClick={() => setPage((p) => Math.max(0, p - 1))}
            >
              Previous
            </Button>
            <Button
              variant="outline"
              size="sm"
              disabled={safePage >= totalPages - 1}
              onClick={() => setPage((p) => Math.min(totalPages - 1, p + 1))}
            >
              Next
            </Button>
          </div>
        </div>
      )}

      {/* Enrichment panel (below table) */}
      {showEnrichment && enrichmentGenes.length > 0 && (
        <EnrichmentPanel
          genes={enrichmentGenes}
          universeSize={allGenes.length}
          onClose={() => setShowEnrichment(false)}
        />
      )}
    </div>
  );
}
