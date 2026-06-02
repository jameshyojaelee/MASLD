"use client";

import { useEffect, useMemo, useState, useCallback } from "react";
import Link from "next/link";
import { useRouter } from "next/navigation";
import { Button } from "@/components/ui/button";
import { Badge } from "@/components/ui/badge";

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

interface CrossSpeciesGene {
  symbol: string;
  human_logfc: number;
  mouse_logfc: number;
  category: string;
  n_concordant: number;
  n_diets_sig: number;
  translatability_score: number;
  diets_concordant: string;
}

interface PathwayConcordance {
  pathway: string;
  pathway_id: string;
  human_nes: number;
  mouse_nes: number;
  concordant: boolean;
}

interface CrossSpeciesSummary {
  total_orthologs: number;
  conserved_core: number;
  human_specific: number;
  mouse_specific: number;
  moderate_concordance: number;
  diet_selective: number;
  discordant: number;
  not_significant: number;
  diet_models: string[];
}

interface CrossSpeciesData {
  genes: CrossSpeciesGene[];
  summary: CrossSpeciesSummary;
  pathway_concordance: PathwayConcordance[];
}

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

const CATEGORY_CONFIG: Record<
  string,
  { label: string; color: string; dotClass: string }
> = {
  Conserved_Core: {
    label: "Conserved Core",
    color: "#22c55e",
    dotClass: "bg-green-500",
  },
  Moderate_Concordance: {
    label: "Moderate Concordance",
    color: "#86efac",
    dotClass: "bg-green-300",
  },
  Human_Enriched: {
    label: "Human Enriched",
    color: "#3b82f6",
    dotClass: "bg-blue-500",
  },
  Mouse_Specific: {
    label: "Mouse Specific",
    color: "#ec4899",
    dotClass: "bg-pink-500",
  },
  Diet_Selective: {
    label: "Diet Selective",
    color: "#f59e0b",
    dotClass: "bg-amber-500",
  },
  Species_Discordant: {
    label: "Discordant",
    color: "#ef4444",
    dotClass: "bg-red-500",
  },
  Not_Significant: {
    label: "Not Significant",
    color: "#9ca3af",
    dotClass: "bg-gray-400",
  },
  Unclassified: {
    label: "Unclassified",
    color: "#6b7280",
    dotClass: "bg-gray-500",
  },
};

const ALL_CATEGORIES = Object.keys(CATEGORY_CONFIG);

const PAGE_SIZE = 50;

type GeneSortKey =
  | "symbol"
  | "human_logfc"
  | "mouse_logfc"
  | "translatability_score"
  | "n_concordant";
type SortDir = "asc" | "desc";

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

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

function categoryBadgeVariant(
  cat: string
): "default" | "secondary" | "outline" | "destructive" {
  switch (cat) {
    case "Conserved_Core":
      return "default";
    case "Species_Discordant":
      return "destructive";
    case "Human_Enriched":
    case "Mouse_Specific":
      return "secondary";
    default:
      return "outline";
  }
}

// ---------------------------------------------------------------------------
// Scatterplot component
// ---------------------------------------------------------------------------

const PLOT_W = 500;
const PLOT_H = 500;
const MARGIN = { top: 20, right: 20, bottom: 45, left: 55 };
const INNER_W = PLOT_W - MARGIN.left - MARGIN.right;
const INNER_H = PLOT_H - MARGIN.top - MARGIN.bottom;

function Scatterplot({
  genes,
  activeCategories,
}: {
  genes: CrossSpeciesGene[];
  activeCategories: Set<string>;
}) {
  const router = useRouter();

  // Compute axis bounds: symmetric range covering data
  const { xMin, xMax, yMin, yMax } = useMemo(() => {
    const xVals = genes.map((g) => g.human_logfc);
    const yVals = genes.map((g) => g.mouse_logfc);
    const xAbs = Math.ceil(Math.max(Math.abs(Math.min(...xVals)), Math.abs(Math.max(...xVals))) * 10) / 10 + 0.3;
    const yAbs = Math.ceil(Math.max(Math.abs(Math.min(...yVals)), Math.abs(Math.max(...yVals))) * 10) / 10 + 0.3;
    return { xMin: -xAbs, xMax: xAbs, yMin: -yAbs, yMax: yAbs };
  }, [genes]);

  const scaleX = useCallback(
    (v: number) => ((v - xMin) / (xMax - xMin)) * INNER_W,
    [xMin, xMax]
  );
  const scaleY = useCallback(
    (v: number) => INNER_H - ((v - yMin) / (yMax - yMin)) * INNER_H,
    [yMin, yMax]
  );

  // Build ticks
  const xTicks = useMemo(() => {
    const ticks: number[] = [];
    const step = xMax > 2 ? 1 : 0.5;
    for (let v = Math.ceil(xMin / step) * step; v <= xMax; v += step) {
      ticks.push(Math.round(v * 10) / 10);
    }
    return ticks;
  }, [xMin, xMax]);

  const yTicks = useMemo(() => {
    const ticks: number[] = [];
    const step = yMax > 2 ? 1 : 0.5;
    for (let v = Math.ceil(yMin / step) * step; v <= yMax; v += step) {
      ticks.push(Math.round(v * 10) / 10);
    }
    return ticks;
  }, [yMin, yMax]);

  // Separate visible (active) and dimmed (inactive) genes
  const { visible, dimmed } = useMemo(() => {
    const vis: CrossSpeciesGene[] = [];
    const dim: CrossSpeciesGene[] = [];
    for (const g of genes) {
      if (activeCategories.has(g.category)) {
        vis.push(g);
      } else {
        dim.push(g);
      }
    }
    return { visible: vis, dimmed: dim };
  }, [genes, activeCategories]);

  return (
    <svg
      viewBox={`0 0 ${PLOT_W} ${PLOT_H}`}
      className="w-full max-w-[500px]"
      role="img"
      aria-label="Human vs Mouse logFC scatterplot"
    >
      <g transform={`translate(${MARGIN.left},${MARGIN.top})`}>
        {/* Grid lines */}
        {xTicks.map((t) => (
          <line
            key={`xg-${t}`}
            x1={scaleX(t)}
            x2={scaleX(t)}
            y1={0}
            y2={INNER_H}
            className="stroke-border"
            strokeWidth={t === 0 ? 1.5 : 0.5}
            strokeDasharray={t === 0 ? undefined : "2,2"}
          />
        ))}
        {yTicks.map((t) => (
          <line
            key={`yg-${t}`}
            x1={0}
            x2={INNER_W}
            y1={scaleY(t)}
            y2={scaleY(t)}
            className="stroke-border"
            strokeWidth={t === 0 ? 1.5 : 0.5}
            strokeDasharray={t === 0 ? undefined : "2,2"}
          />
        ))}

        {/* Quadrant labels */}
        <text
          x={INNER_W * 0.75}
          y={INNER_H * 0.08}
          textAnchor="middle"
          className="fill-muted-foreground text-[10px]"
        >
          Both Up
        </text>
        <text
          x={INNER_W * 0.25}
          y={INNER_H * 0.92}
          textAnchor="middle"
          className="fill-muted-foreground text-[10px]"
        >
          Both Down
        </text>
        <text
          x={INNER_W * 0.25}
          y={INNER_H * 0.08}
          textAnchor="middle"
          className="fill-muted-foreground/60 text-[10px]"
        >
          Discordant
        </text>
        <text
          x={INNER_W * 0.75}
          y={INNER_H * 0.92}
          textAnchor="middle"
          className="fill-muted-foreground/60 text-[10px]"
        >
          Discordant
        </text>

        {/* Dimmed dots (inactive categories) */}
        {dimmed.map((g) => (
          <circle
            key={`d-${g.symbol}`}
            cx={scaleX(g.human_logfc)}
            cy={scaleY(g.mouse_logfc)}
            r={2}
            fill="#9ca3af"
            opacity={0.12}
          />
        ))}

        {/* Active dots */}
        {visible.map((g) => {
          const cfg = CATEGORY_CONFIG[g.category];
          return (
            <circle
              key={g.symbol}
              cx={scaleX(g.human_logfc)}
              cy={scaleY(g.mouse_logfc)}
              r={g.category === "Conserved_Core" ? 3.5 : 2.5}
              fill={cfg?.color ?? "#6b7280"}
              opacity={0.75}
              className="cursor-pointer transition-opacity hover:opacity-100"
              onClick={() =>
                router.push(`/gene/${encodeURIComponent(g.symbol)}`)
              }
            >
              <title>
                {g.symbol} | Human: {formatLogFC(g.human_logfc)} | Mouse:{" "}
                {formatLogFC(g.mouse_logfc)} | {cfg?.label ?? g.category}
              </title>
            </circle>
          );
        })}

        {/* X axis ticks + labels */}
        {xTicks.map((t) => (
          <g key={`xt-${t}`} transform={`translate(${scaleX(t)},${INNER_H})`}>
            <line y2={4} className="stroke-muted-foreground" strokeWidth={0.5} />
            <text
              y={14}
              textAnchor="middle"
              className="fill-muted-foreground text-[10px]"
            >
              {t}
            </text>
          </g>
        ))}

        {/* Y axis ticks + labels */}
        {yTicks.map((t) => (
          <g key={`yt-${t}`} transform={`translate(0,${scaleY(t)})`}>
            <line x2={-4} className="stroke-muted-foreground" strokeWidth={0.5} />
            <text
              x={-8}
              textAnchor="end"
              dominantBaseline="middle"
              className="fill-muted-foreground text-[10px]"
            >
              {t}
            </text>
          </g>
        ))}

        {/* Axis labels */}
        <text
          x={INNER_W / 2}
          y={INNER_H + 36}
          textAnchor="middle"
          className="fill-foreground text-xs font-medium"
        >
          Human logFC
        </text>
        <text
          x={-INNER_H / 2}
          y={-42}
          textAnchor="middle"
          transform="rotate(-90)"
          className="fill-foreground text-xs font-medium"
        >
          Mouse logFC
        </text>
      </g>
    </svg>
  );
}

// ---------------------------------------------------------------------------
// NES bar helper
// ---------------------------------------------------------------------------

function NesBar({ value, maxAbs }: { value: number; maxAbs: number }) {
  const pct = Math.min(Math.abs(value) / maxAbs, 1) * 100;
  const isPos = value >= 0;
  return (
    <div className="flex items-center gap-1.5">
      <div className="relative h-3 w-20 overflow-hidden rounded-sm bg-muted">
        {isPos ? (
          <div
            className="absolute top-0 left-1/2 h-full rounded-sm bg-red-400/70"
            style={{ width: `${pct / 2}%` }}
          />
        ) : (
          <div
            className="absolute top-0 h-full rounded-sm bg-blue-400/70"
            style={{
              width: `${pct / 2}%`,
              right: "50%",
            }}
          />
        )}
      </div>
      <span
        className={`font-mono text-xs ${isPos ? "text-red-500 dark:text-red-400" : "text-blue-500 dark:text-blue-400"}`}
      >
        {value >= 0 ? "+" : ""}
        {value.toFixed(2)}
      </span>
    </div>
  );
}

// ---------------------------------------------------------------------------
// Main page
// ---------------------------------------------------------------------------

export default function SpeciesPage() {
  const [data, setData] = useState<CrossSpeciesData | null>(null);
  const [loading, setLoading] = useState(true);
  const [activeCategories, setActiveCategories] = useState<Set<string>>(
    () => new Set(ALL_CATEGORIES)
  );
  const [geneSortKey, setGeneSortKey] = useState<GeneSortKey>(
    "translatability_score"
  );
  const [geneSortDir, setGeneSortDir] = useState<SortDir>("desc");
  const [genePage, setGenePage] = useState(0);

  useEffect(() => {
    fetch("/data/cross_species.json")
      .then((r) => r.json())
      .then((d: CrossSpeciesData) => {
        setData(d);
        setLoading(false);
      })
      .catch(() => setLoading(false));
  }, []);

  // Category toggle
  const toggleCategory = useCallback((cat: string) => {
    setActiveCategories((prev) => {
      const next = new Set(prev);
      if (next.has(cat)) {
        next.delete(cat);
      } else {
        next.add(cat);
      }
      return next;
    });
    setGenePage(0);
  }, []);

  const selectOnlyCategory = useCallback((cat: string) => {
    setActiveCategories(new Set([cat]));
    setGenePage(0);
  }, []);

  const selectAllCategories = useCallback(() => {
    setActiveCategories(new Set(ALL_CATEGORIES));
    setGenePage(0);
  }, []);

  // Filtered + sorted genes
  const filteredGenes = useMemo(() => {
    if (!data) return [];
    return data.genes.filter((g) => activeCategories.has(g.category));
  }, [data, activeCategories]);

  const sortedGenes = useMemo(() => {
    return [...filteredGenes].sort((a, b) => {
      let av: string | number;
      let bv: string | number;
      switch (geneSortKey) {
        case "symbol":
          av = a.symbol;
          bv = b.symbol;
          break;
        case "human_logfc":
          av = a.human_logfc;
          bv = b.human_logfc;
          break;
        case "mouse_logfc":
          av = a.mouse_logfc;
          bv = b.mouse_logfc;
          break;
        case "translatability_score":
          av = a.translatability_score;
          bv = b.translatability_score;
          break;
        case "n_concordant":
          av = a.n_concordant;
          bv = b.n_concordant;
          break;
      }
      let cmp: number;
      if (typeof av === "string" && typeof bv === "string") {
        cmp = av.localeCompare(bv);
      } else {
        cmp = (av as number) - (bv as number);
      }
      return geneSortDir === "asc" ? cmp : -cmp;
    });
  }, [filteredGenes, geneSortKey, geneSortDir]);

  const totalPages = Math.max(1, Math.ceil(sortedGenes.length / PAGE_SIZE));
  const safePage = Math.min(genePage, totalPages - 1);
  const pageGenes = sortedGenes.slice(
    safePage * PAGE_SIZE,
    (safePage + 1) * PAGE_SIZE
  );

  // Pathway sorted by |human_nes|
  const sortedPathways = useMemo(() => {
    if (!data) return [];
    return [...data.pathway_concordance].sort(
      (a, b) => Math.abs(b.human_nes) - Math.abs(a.human_nes)
    );
  }, [data]);

  const maxPathwayNes = useMemo(() => {
    if (!sortedPathways.length) return 3;
    return Math.max(
      ...sortedPathways.map((p) =>
        Math.max(Math.abs(p.human_nes), Math.abs(p.mouse_nes))
      )
    );
  }, [sortedPathways]);

  const handleGeneSort = useCallback(
    (key: GeneSortKey) => {
      if (geneSortKey === key) {
        setGeneSortDir((d) => (d === "asc" ? "desc" : "asc"));
      } else {
        setGeneSortKey(key);
        setGeneSortDir(key === "symbol" ? "asc" : "desc");
      }
      setGenePage(0);
    },
    [geneSortKey]
  );

  const sortIndicator = (key: GeneSortKey) => {
    if (geneSortKey !== key) return null;
    return geneSortDir === "asc" ? " \u25B2" : " \u25BC";
  };

  // Loading state
  if (loading) {
    return (
      <div className="mx-auto max-w-6xl px-6 py-10">
        <h1 className="text-3xl font-bold tracking-tight">
          Cross-Species Mirror
        </h1>
        <p className="mt-6 text-sm text-muted-foreground">
          Loading concordance data...
        </p>
      </div>
    );
  }

  if (!data) {
    return (
      <div className="mx-auto max-w-6xl px-6 py-10">
        <h1 className="text-3xl font-bold tracking-tight">
          Cross-Species Mirror
        </h1>
        <p className="mt-6 text-sm text-destructive">
          Failed to load cross-species data.
        </p>
      </div>
    );
  }

  const { summary } = data;

  return (
    <div className="mx-auto max-w-6xl px-6 py-8">
      {/* Header */}
      <div className="mb-8">
        <h1 className="text-3xl font-bold tracking-tight">
          Cross-Species Mirror
        </h1>
        <p className="mt-2 text-muted-foreground">
          Human-mouse concordance validation across {summary.total_orthologs.toLocaleString()} orthologs,{" "}
          {summary.diet_models.length} diet models, and{" "}
          {summary.conserved_core.toLocaleString()} Conserved Core genes.
        </p>
      </div>

      {/* Section 1: Summary Stats */}
      <div className="mb-10 grid grid-cols-2 gap-4 sm:grid-cols-3 lg:grid-cols-5">
        <StatCard
          label="Total Orthologs"
          value={summary.total_orthologs.toLocaleString()}
        />
        <StatCard
          label="Conserved Core"
          value={summary.conserved_core.toLocaleString()}
          accent="green"
        />
        <StatCard
          label="Human Enriched"
          value={summary.human_specific.toLocaleString()}
          accent="blue"
        />
        <StatCard
          label="Discordant"
          value={summary.discordant.toLocaleString()}
          accent="red"
        />
        <StatCard
          label="Diet Models"
          value={summary.diet_models.join(", ")}
          small
        />
      </div>

      {/* Section 2 + 3: Scatterplot with category filters */}
      <div className="mb-10">
        <h2 className="mb-4 text-xl font-semibold">
          Concordance Scatterplot
        </h2>

        {/* Category filters */}
        <div className="mb-4 flex flex-wrap items-center gap-2">
          <Button
            variant={
              activeCategories.size === ALL_CATEGORIES.length
                ? "default"
                : "outline"
            }
            size="sm"
            onClick={selectAllCategories}
          >
            All
          </Button>
          {ALL_CATEGORIES.map((cat) => {
            const cfg = CATEGORY_CONFIG[cat];
            const isActive = activeCategories.has(cat);
            const count = data.genes.filter(
              (g) => g.category === cat
            ).length;
            return (
              <Button
                key={cat}
                variant={isActive ? "default" : "outline"}
                size="sm"
                onClick={() => toggleCategory(cat)}
                onDoubleClick={() => selectOnlyCategory(cat)}
              >
                <span
                  className={`mr-1.5 inline-block h-2.5 w-2.5 rounded-full ${cfg.dotClass}`}
                  style={{ opacity: isActive ? 1 : 0.3 }}
                />
                {cfg.label}
                <span className="ml-1 text-xs text-muted-foreground">
                  ({count})
                </span>
              </Button>
            );
          })}
        </div>
        <p className="mb-3 text-xs text-muted-foreground">
          Click a gene point to view its profile. Double-click a filter to
          isolate that category.
        </p>

        {/* Scatterplot */}
        <div className="flex justify-center">
          <Scatterplot
            genes={data.genes}
            activeCategories={activeCategories}
          />
        </div>
      </div>

      {/* Section 4: Pathway Concordance */}
      <div className="mb-10">
        <h2 className="mb-4 text-xl font-semibold">Pathway Concordance</h2>
        <p className="mb-3 text-sm text-muted-foreground">
          Hallmark pathway enrichment (NES) in human vs mouse. Sorted by
          |Human NES|.
        </p>
        <div className="overflow-x-auto rounded-lg border border-border">
          <table className="w-full text-sm">
            <thead>
              <tr className="border-b border-border bg-muted/50">
                <th className="px-3 py-2 text-left text-xs font-medium text-muted-foreground">
                  Pathway
                </th>
                <th className="px-3 py-2 text-left text-xs font-medium text-muted-foreground">
                  Human NES
                </th>
                <th className="px-3 py-2 text-left text-xs font-medium text-muted-foreground">
                  Mouse NES
                </th>
                <th className="px-3 py-2 text-center text-xs font-medium text-muted-foreground">
                  Concordant
                </th>
              </tr>
            </thead>
            <tbody>
              {sortedPathways.map((pw) => (
                <tr
                  key={pw.pathway_id}
                  className="border-b border-border/50 transition-colors hover:bg-muted/30"
                >
                  <td className="px-3 py-1.5 text-xs font-medium">
                    {pw.pathway}
                  </td>
                  <td className="px-3 py-1.5">
                    <NesBar value={pw.human_nes} maxAbs={maxPathwayNes} />
                  </td>
                  <td className="px-3 py-1.5">
                    <NesBar value={pw.mouse_nes} maxAbs={maxPathwayNes} />
                  </td>
                  <td className="px-3 py-1.5 text-center">
                    {pw.concordant ? (
                      <span
                        className="text-green-600 dark:text-green-400"
                        title="Concordant"
                      >
                        &#10003;
                      </span>
                    ) : (
                      <span
                        className="text-red-500 dark:text-red-400"
                        title="Discordant"
                      >
                        &#10007;
                      </span>
                    )}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </div>

      {/* Section 5: Gene List */}
      <div className="mb-10">
        <h2 className="mb-4 text-xl font-semibold">Gene List</h2>
        <div className="mb-3 text-sm text-muted-foreground">
          Showing{" "}
          <span className="font-medium text-foreground">
            {sortedGenes.length.toLocaleString()}
          </span>{" "}
          of{" "}
          <span className="font-medium text-foreground">
            {data.genes.length.toLocaleString()}
          </span>{" "}
          genes matching active filters.
        </div>

        <div className="overflow-x-auto rounded-lg border border-border">
          <table className="w-full text-sm">
            <thead>
              <tr className="border-b border-border bg-muted/50">
                {(
                  [
                    { key: "symbol" as GeneSortKey, label: "Symbol", align: "text-left" },
                    {
                      key: "human_logfc" as GeneSortKey,
                      label: "Human logFC",
                      align: "text-right",
                    },
                    {
                      key: "mouse_logfc" as GeneSortKey,
                      label: "Mouse logFC",
                      align: "text-right",
                    },
                    {
                      key: "translatability_score" as GeneSortKey,
                      label: "Translatability",
                      align: "text-right",
                    },
                    {
                      key: "n_concordant" as GeneSortKey,
                      label: "Diets Concordant",
                      align: "text-right",
                    },
                  ] as const
                ).map((col) => (
                  <th
                    key={col.key}
                    className={`cursor-pointer select-none px-3 py-2 text-xs font-medium text-muted-foreground hover:text-foreground ${col.align}`}
                    onClick={() => handleGeneSort(col.key)}
                  >
                    {col.label}
                    {sortIndicator(col.key)}
                  </th>
                ))}
                <th className="px-3 py-2 text-left text-xs font-medium text-muted-foreground">
                  Category
                </th>
              </tr>
            </thead>
            <tbody>
              {pageGenes.length === 0 ? (
                <tr>
                  <td
                    colSpan={6}
                    className="px-3 py-12 text-center text-muted-foreground"
                  >
                    No genes match the current filters.
                  </td>
                </tr>
              ) : (
                pageGenes.map((gene) => {
                  const cfg = CATEGORY_CONFIG[gene.category];
                  return (
                    <tr
                      key={gene.symbol}
                      className="border-b border-border/50 transition-colors hover:bg-muted/30"
                    >
                      <td className="px-3 py-1.5">
                        <Link
                          href={`/gene/${encodeURIComponent(gene.symbol)}`}
                          className="font-mono font-semibold text-primary hover:underline"
                        >
                          {gene.symbol}
                        </Link>
                      </td>
                      <td
                        className={`px-3 py-1.5 text-right font-mono text-xs ${logfcColor(gene.human_logfc)}`}
                      >
                        {formatLogFC(gene.human_logfc)}
                      </td>
                      <td
                        className={`px-3 py-1.5 text-right font-mono text-xs ${logfcColor(gene.mouse_logfc)}`}
                      >
                        {formatLogFC(gene.mouse_logfc)}
                      </td>
                      <td className="px-3 py-1.5 text-right font-mono text-xs">
                        {gene.translatability_score.toFixed(3)}
                      </td>
                      <td className="px-3 py-1.5 text-right font-mono text-xs">
                        {gene.n_concordant}/{summary.diet_models.length}
                      </td>
                      <td className="px-3 py-1.5">
                        <Badge
                          variant={categoryBadgeVariant(gene.category)}
                          className="text-[10px]"
                        >
                          {cfg?.label ?? gene.category}
                        </Badge>
                      </td>
                    </tr>
                  );
                })
              )}
            </tbody>
          </table>
        </div>

        {/* Pagination */}
        {sortedGenes.length > PAGE_SIZE && (
          <div className="mt-4 flex items-center justify-between">
            <p className="text-xs text-muted-foreground">
              Page {safePage + 1} of {totalPages}
            </p>
            <div className="flex gap-2">
              <Button
                variant="outline"
                size="sm"
                disabled={safePage === 0}
                onClick={() => setGenePage((p) => Math.max(0, p - 1))}
              >
                Previous
              </Button>
              <Button
                variant="outline"
                size="sm"
                disabled={safePage >= totalPages - 1}
                onClick={() =>
                  setGenePage((p) => Math.min(totalPages - 1, p + 1))
                }
              >
                Next
              </Button>
            </div>
          </div>
        )}
      </div>
    </div>
  );
}

// ---------------------------------------------------------------------------
// Stat card sub-component
// ---------------------------------------------------------------------------

function StatCard({
  label,
  value,
  accent,
  small,
}: {
  label: string;
  value: string;
  accent?: "green" | "blue" | "red";
  small?: boolean;
}) {
  const accentClass =
    accent === "green"
      ? "text-green-600 dark:text-green-400"
      : accent === "blue"
        ? "text-blue-600 dark:text-blue-400"
        : accent === "red"
          ? "text-red-600 dark:text-red-400"
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
