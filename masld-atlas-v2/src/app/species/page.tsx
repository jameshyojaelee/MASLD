"use client";

import { useEffect, useMemo, useState, useCallback, type CSSProperties } from "react";
import { HashLink as Link } from "@/components/hash-link";
import { useHashNavigate } from "@/lib/hash-router";
import { Button } from "@/components/ui/button";
import { Badge } from "@/components/ui/badge";
import { PageContainer } from "@/components/page-container";
import { PageHeader } from "@/components/page-header";
import { SkeletonBlock } from "@/components/states";
import { Scatter, type ScatterPoint } from "@/components/charts";
import { dataUrl } from "@/lib/data-base";
import { speciesCategoryColor, categoricalColor, CONTROL } from "@/lib/palette";

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
// Category config. Colors flow through the palette authority:
// speciesCategoryColor maps the canonical categories (Not_Significant /
// Unclassified -> CONTROL); Moderate_Concordance is not a palette category, so
// it takes a distinct qualitative slot (olive) — never control gray.
// ---------------------------------------------------------------------------

const CATEGORY_LABELS: Record<string, string> = {
  Conserved_Core: "Conserved Core",
  Moderate_Concordance: "Moderate Concordance",
  Human_Enriched: "Human Enriched",
  Mouse_Specific: "Mouse Specific",
  Diet_Selective: "Diet Selective",
  Species_Discordant: "Discordant",
  Not_Significant: "Not Significant",
  Unclassified: "Unclassified",
};

const ALL_CATEGORIES = Object.keys(CATEGORY_LABELS);

function categoryColor(cat: string): string {
  if (cat === "Moderate_Concordance") return categoricalColor(9); // olive
  return speciesCategoryColor(cat);
}

function categoryLabel(cat: string): string {
  return CATEGORY_LABELS[cat] ?? cat;
}

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
  if (val == null) return "—";
  return val >= 0 ? `+${val.toFixed(2)}` : val.toFixed(2);
}

/** Sign-encoded logFC color for table cells (direction = data mark). */
function lfcStyle(val: number | null | undefined): CSSProperties {
  if (val == null || val === 0) return { color: "var(--color-muted-foreground)" };
  return { color: val > 0 ? "var(--color-effect-up)" : "var(--color-effect-down)" };
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
// NES bar helper (diverging: down = blue, up = red, via effect tokens)
// ---------------------------------------------------------------------------

function NesBar({ value, maxAbs }: { value: number; maxAbs: number }) {
  const pct = Math.min(Math.abs(value) / maxAbs, 1) * 100;
  const isPos = value >= 0;
  const barColor = isPos ? "var(--color-effect-up)" : "var(--color-effect-down)";
  return (
    <div className="flex items-center gap-1.5">
      <div className="relative h-3 w-20 overflow-hidden rounded-sm bg-muted">
        <div
          className="absolute top-0 h-full rounded-sm"
          style={{
            width: `${pct / 2}%`,
            backgroundColor: barColor,
            opacity: 0.7,
            ...(isPos ? { left: "50%" } : { right: "50%" }),
          }}
        />
      </div>
      <span className="font-mono text-xs" style={{ color: barColor }}>
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
  const navigate = useHashNavigate();
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
    fetch(dataUrl("cross_species.json"))
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

  // Scatter points (color via palette; inactive categories dim by opacity).
  const scatterPoints = useMemo<ScatterPoint[]>(() => {
    if (!data) return [];
    return data.genes.map((g) => ({
      x: g.human_logfc,
      y: g.mouse_logfc,
      label: g.symbol,
      color: categoryColor(g.category),
      category: g.category,
      dimmed: !activeCategories.has(g.category),
      emphasized: g.category === "Conserved_Core",
    }));
  }, [data, activeCategories]);

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
    return geneSortDir === "asc" ? " ▲" : " ▼";
  };

  const header = (
    <PageHeader
      eyebrow="Cross-species validation"
      title="Cross-Species Mirror"
      description={
        data
          ? `Human-mouse concordance validation across ${data.summary.total_orthologs.toLocaleString()} orthologs, ${data.summary.diet_models.length} diet models, and ${data.summary.conserved_core.toLocaleString()} Conserved Core genes.`
          : "Human-mouse concordance validation across orthologs and diet models."
      }
    />
  );

  // Loading state
  if (loading) {
    return (
      <PageContainer>
        {header}
        <div className="grid grid-cols-2 gap-4 sm:grid-cols-3 lg:grid-cols-5">
          {Array.from({ length: 5 }).map((_, i) => (
            <SkeletonBlock key={i} className="h-20 rounded-lg" />
          ))}
        </div>
        <SkeletonBlock className="mt-8 h-[420px] w-full max-w-2xl rounded-lg" />
      </PageContainer>
    );
  }

  if (!data) {
    return (
      <PageContainer>
        {header}
        <p className="mt-6 text-sm text-destructive">
          Failed to load cross-species data.
        </p>
      </PageContainer>
    );
  }

  const { summary } = data;

  return (
    <PageContainer>
      {header}

      {/* Section 1: Summary Stats */}
      <div className="mb-10 grid grid-cols-2 gap-4 sm:grid-cols-3 lg:grid-cols-5">
        <StatCard
          label="Total Orthologs"
          value={summary.total_orthologs.toLocaleString()}
        />
        <StatCard
          label="Conserved Core"
          value={summary.conserved_core.toLocaleString()}
          accentColor={categoryColor("Conserved_Core")}
        />
        <StatCard
          label="Human Enriched"
          value={summary.human_specific.toLocaleString()}
          accentColor={categoryColor("Human_Enriched")}
        />
        <StatCard
          label="Discordant"
          value={summary.discordant.toLocaleString()}
          accentColor={categoryColor("Species_Discordant")}
        />
        <StatCard
          label="Diet Models"
          value={summary.diet_models.join(", ")}
          small
        />
      </div>

      {/* Section 2 + 3: Scatterplot with category filters */}
      <div className="mb-10">
        <h2 className="mb-4 text-xl font-semibold">Concordance Scatterplot</h2>

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
            const isActive = activeCategories.has(cat);
            const count = data.genes.filter((g) => g.category === cat).length;
            return (
              <Button
                key={cat}
                variant={isActive ? "default" : "outline"}
                size="sm"
                onClick={() => toggleCategory(cat)}
                onDoubleClick={() => selectOnlyCategory(cat)}
              >
                <span
                  className="mr-1.5 inline-block h-2.5 w-2.5 rounded-full"
                  style={{
                    backgroundColor: categoryColor(cat),
                    opacity: isActive ? 1 : 0.3,
                  }}
                />
                {categoryLabel(cat)}
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
        <div className="mx-auto w-full max-w-2xl">
          <Scatter
            data={scatterPoints}
            xLabel="Human logFC"
            yLabel="Mouse logFC"
            height={460}
            ariaLabel="Human vs mouse logFC concordance scatterplot"
            quadrantLabels={[
              { x: "right", y: "top", text: "Both Up" },
              { x: "left", y: "bottom", text: "Both Down" },
              { x: "left", y: "top", text: "Discordant", muted: true },
              { x: "right", y: "bottom", text: "Discordant", muted: true },
            ]}
            onPointClick={(sym) =>
              navigate(`#/gene?symbol=${encodeURIComponent(sym)}`)
            }
            tooltipLines={(p) => (
              <>
                <div className="font-medium italic">{p.label}</div>
                <div>
                  Human {formatLogFC(p.x)} · Mouse {formatLogFC(p.y)}
                </div>
                <div className="text-muted-foreground">
                  {categoryLabel(p.category ?? "")}
                </div>
              </>
            )}
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
                        style={{ color: categoricalColor(2) }}
                        title="Concordant"
                      >
                        &#10003;
                      </span>
                    ) : (
                      <span style={{ color: CONTROL }} title="Discordant">
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
                pageGenes.map((gene) => (
                  <tr
                    key={gene.symbol}
                    className="border-b border-border/50 transition-colors hover:bg-muted/30"
                  >
                    <td className="px-3 py-1.5">
                      <Link
                        href={`#/gene?symbol=${encodeURIComponent(gene.symbol)}`}
                        className="font-mono font-semibold text-primary hover:underline"
                      >
                        {gene.symbol}
                      </Link>
                    </td>
                    <td
                      className="px-3 py-1.5 text-right font-mono text-xs"
                      style={lfcStyle(gene.human_logfc)}
                    >
                      {formatLogFC(gene.human_logfc)}
                    </td>
                    <td
                      className="px-3 py-1.5 text-right font-mono text-xs"
                      style={lfcStyle(gene.mouse_logfc)}
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
                        {categoryLabel(gene.category)}
                      </Badge>
                    </td>
                  </tr>
                ))
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
    </PageContainer>
  );
}

// ---------------------------------------------------------------------------
// Stat card sub-component
// ---------------------------------------------------------------------------

function StatCard({
  label,
  value,
  accentColor,
  small,
}: {
  label: string;
  value: string;
  accentColor?: string;
  small?: boolean;
}) {
  return (
    <div className="rounded-lg border border-border bg-card p-4 hover-lift">
      <p className="text-xs font-medium text-muted-foreground">{label}</p>
      <p
        className={`mt-1 font-semibold ${small ? "text-sm" : "text-2xl"} ${accentColor ? "" : "text-foreground"}`}
        style={accentColor ? { color: accentColor } : undefined}
      >
        {value}
      </p>
    </div>
  );
}
