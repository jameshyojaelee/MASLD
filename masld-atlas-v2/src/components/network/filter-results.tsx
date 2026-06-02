"use client";

import { useCallback, useMemo } from "react";
import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import type { GeneSummary } from "@/lib/network-filters";

interface FilterResultsProps {
  loading: boolean;
  results: GeneSummary[];
  unresolved: string[];
  selectedGene?: string | null;
  onSelect: (symbol: string) => void;
  onClear: () => void;
}

const MAX_DISPLAY = 300;

function toCsv(rows: GeneSummary[]): string {
  const header = [
    "symbol",
    "degree",
    "is_deg",
    "is_coloc",
    "is_drug",
    "c_macro_F01",
    "c_macro_F34",
  ].join(",");
  const lines = rows.map((r) =>
    [
      r.symbol,
      r.degree ?? "",
      r.is_deg ?? "",
      r.is_coloc ?? "",
      r.is_drug ?? "",
      r.c_macro_F01 ?? "",
      r.c_macro_F34 ?? "",
    ].join(",")
  );
  return [header, ...lines].join("\n");
}

export function FilterResults({
  loading,
  results,
  unresolved,
  selectedGene,
  onSelect,
  onClear,
}: FilterResultsProps) {
  const display = useMemo(() => results.slice(0, MAX_DISPLAY), [results]);

  const handleDownload = useCallback(() => {
    const csv = toCsv(results);
    const blob = new Blob([csv], { type: "text/csv;charset=utf-8" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `masld_network_filter_${Date.now()}.csv`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  }, [results]);

  return (
    <section className="mt-4 border-t border-border pt-3">
      <div className="mb-2 flex items-center justify-between">
        <h3 className="text-xs font-medium uppercase text-muted-foreground">
          Results
        </h3>
        <div className="flex gap-1">
          <Button
            size="sm"
            variant="outline"
            onClick={handleDownload}
            disabled={results.length === 0}
            title="Download matching genes as CSV"
          >
            CSV
          </Button>
          <Button size="sm" variant="ghost" onClick={onClear} title="Clear filters">
            Clear
          </Button>
        </div>
      </div>

      <div className="mb-2 text-xs">
        {loading ? (
          <span className="text-muted-foreground animate-pulse">
            Filtering...
          </span>
        ) : (
          <span>
            <span className="font-semibold">{results.length.toLocaleString()}</span>{" "}
            <span className="text-muted-foreground">
              gene{results.length === 1 ? "" : "s"} matching filters
            </span>
          </span>
        )}
      </div>

      {unresolved.length > 0 && (
        <div className="mb-2 rounded border border-amber-500/40 bg-amber-500/10 px-2 py-1 text-[11px] text-amber-700 dark:text-amber-300">
          Awaiting filter indexes for: {unresolved.join(", ")}. Retrying...
        </div>
      )}

      <ul className="max-h-64 overflow-y-auto rounded border border-border bg-background/60 text-xs">
        {display.length === 0 && !loading && (
          <li className="px-2 py-3 text-center text-muted-foreground">
            No genes match current filters.
          </li>
        )}
        {display.map((g) => {
          const active = selectedGene === g.symbol;
          return (
            <li key={g.symbol}>
              <button
                type="button"
                onClick={() => onSelect(g.symbol)}
                className={
                  "flex w-full items-center justify-between gap-2 px-2 py-1 text-left transition-colors hover:bg-muted " +
                  (active ? "bg-muted font-semibold" : "")
                }
              >
                <span className="font-mono">{g.symbol}</span>
                <span className="flex items-center gap-1">
                  {g.is_deg && (
                    <Badge variant="secondary" className="text-[9px] px-1 py-0">
                      DEG
                    </Badge>
                  )}
                  {g.is_coloc && (
                    <Badge variant="secondary" className="text-[9px] px-1 py-0">
                      COLOC
                    </Badge>
                  )}
                  {g.is_drug && (
                    <Badge variant="secondary" className="text-[9px] px-1 py-0">
                      Drug
                    </Badge>
                  )}
                  {typeof g.degree === "number" && (
                    <span className="text-muted-foreground">{g.degree}</span>
                  )}
                </span>
              </button>
            </li>
          );
        })}
      </ul>

      {results.length > MAX_DISPLAY && (
        <div className="mt-1 text-[11px] text-muted-foreground">
          Showing first {MAX_DISPLAY.toLocaleString()} of{" "}
          {results.length.toLocaleString()}. Use CSV export for the full list.
        </div>
      )}
    </section>
  );
}
