"use client";

/**
 * Landing mini-viz: colocalized effector genes per GWAS ancestry.
 *
 * Self-contained and fail-soft — fetches `coloc_ancestry_summary.json` and
 * reads `by_ancestry[*].n_genes_h4_08` (genes passing PP.H4 >= 0.8). Horizontal
 * SVG bars ordered EUR -> SAS, each filled with its colorblind-safe Okabe-Ito
 * ancestry color; the ancestry code and gene count are drawn in theme-token
 * text (never colored) with tabular monospaced numerals.
 */

import { useEffect, useState } from "react";
import { dataUrl } from "@/lib/data-base";
import { ancestryColor, ANCESTRY_ORDER } from "@/lib/palette";
import { useMeasure } from "./use-measure";

interface AncestryRow {
  ancestry: string;
  n_genes_h4_08: number;
}

const LABEL_W = 26;
const VALUE_W = 30;
const GAP = 4;

export function AncestryBars() {
  const [ref, { width, height }] = useMeasure<HTMLDivElement>();
  const [rows, setRows] = useState<AncestryRow[] | null>(null);
  const [failed, setFailed] = useState(false);

  useEffect(() => {
    let alive = true;
    fetch(dataUrl("coloc_ancestry_summary.json"))
      .then((r) => (r.ok ? r.json() : Promise.reject(r.status)))
      .then((d: { by_ancestry?: AncestryRow[] }) => {
        if (!alive) return;
        const by = d.by_ancestry ?? [];
        const order = new Map(ANCESTRY_ORDER.map((a, i) => [a, i]));
        const sorted = [...by].sort(
          (a, b) =>
            (order.get(a.ancestry as never) ?? 99) -
            (order.get(b.ancestry as never) ?? 99)
        );
        setRows(sorted);
      })
      .catch(() => alive && setFailed(true));
    return () => {
      alive = false;
    };
  }, []);

  if (failed) return null;

  const ready = rows && rows.length > 0 && width > 0 && height > 0;

  return (
    <div ref={ref} className="h-full w-full">
      {!ready ? (
        <div className="h-full w-full animate-pulse rounded-md bg-muted/40 motion-reduce:animate-none" />
      ) : (
        (() => {
          const n = rows!.length;
          const rowH = height / n;
          const barH = Math.max(4, Math.min(14, rowH * 0.52));
          const track = Math.max(1, width - LABEL_W - VALUE_W - 2 * GAP);
          const max = Math.max(1, ...rows!.map((r) => r.n_genes_h4_08));
          const fs = Math.max(8, Math.min(10, rowH * 0.4));
          return (
            <svg
              width={width}
              height={height}
              className="block"
              role="img"
              aria-label="Colocalized effector genes by ancestry"
            >
              {rows!.map((r, i) => {
                const cy = i * rowH + rowH / 2;
                const barW = (r.n_genes_h4_08 / max) * track;
                const barX = LABEL_W + GAP;
                return (
                  <g key={r.ancestry}>
                    <text
                      x={0}
                      y={cy}
                      dominantBaseline="central"
                      fontSize={fs}
                      fill="var(--color-muted-foreground)"
                    >
                      {r.ancestry}
                    </text>
                    <rect
                      x={barX}
                      y={cy - barH / 2}
                      width={Math.max(1, barW)}
                      height={barH}
                      rx={2}
                      fill={ancestryColor(r.ancestry)}
                    />
                    <text
                      x={barX + barW + GAP}
                      y={cy}
                      dominantBaseline="central"
                      fontSize={fs}
                      fill="var(--color-muted-foreground)"
                      style={{
                        fontFamily: "var(--font-mono)",
                        fontVariantNumeric: "tabular-nums",
                      }}
                    >
                      {r.n_genes_h4_08}
                    </text>
                  </g>
                );
              })}
            </svg>
          );
        })()
      )}
    </div>
  );
}
