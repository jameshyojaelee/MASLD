"use client";

/**
 * Regulatory locus-zoom.
 *
 * Zooms into the chromosome + window of a selected fine-mapped variant: y = PIP,
 * x = genomic position (bp), color = colocalization PP.H4 (sequential ramp).
 * The selected variant is the lead (foreground ring); motif-disrupting variants
 * carry a regulatory-hue ring. A compact gene track beneath the axis places each
 * distinct nearest-gene at the span of its variants — honest to the client data
 * (nearest-gene labels; no exon models are shipped).
 *
 * Brush-linked: the caller drives `selected` from the Manhattan; changing it
 * re-centres this view. SCHEMA-FORWARD like {@link Manhattan} — accepts
 * `neglog10p` / `r2_to_lead` for a future genome-wide upgrade.
 */

import { useMemo, useState } from "react";
import { scaleLinear } from "@visx/scale";
import { Text } from "@visx/text";
import { sequentialColor, CONTROL, MODALITY_HEX } from "@/lib/palette";
import type { ManhattanVariant } from "./manhattan";
import { ChartFrame } from "./chart-frame";
import { useChartTooltip, ChartTooltip } from "./chart-tooltip";
import { BottomAxis, LeftAxis, GridRows, AXIS_TEXT, MarkGroup } from "./chart-primitives";

export interface LocusZoomProps {
  variants: ManhattanVariant[];
  selected: ManhattanVariant | null;
  /** Half-window in bp around the selected variant (default 500 kb). */
  windowBp?: number;
  yField?: "pip" | "neglog10p";
  colorField?: "coloc_pp4" | "r2_to_lead";
  onVariantClick?: (v: ManhattanVariant) => void;
  title?: React.ReactNode;
  caption?: React.ReactNode;
  height?: number;
  ariaLabel?: string;
}

const MOTIF_HEX = MODALITY_HEX.s4_epigenomic;

function chrShort(chr: string): string {
  return chr.replace(/^chr/i, "");
}

export function LocusZoom({
  variants,
  selected,
  windowBp = 500_000,
  yField,
  colorField = "coloc_pp4",
  onVariantClick,
  title,
  caption,
  height = 300,
  ariaLabel,
}: LocusZoomProps) {
  const { wrapperRef, tooltip, show, hide } = useChartTooltip<ManhattanVariant>();
  const [hoverId, setHoverId] = useState<string | null>(null);

  const yMode: "pip" | "neglog10p" = useMemo(() => {
    if (yField) return yField;
    return variants.some((v) => v.pip != null) ? "pip" : "neglog10p";
  }, [yField, variants]);
  const yOf = (v: ManhattanVariant) =>
    (yMode === "pip" ? v.pip : v.neglog10p) ?? 0;
  const colorValOf = (v: ManhattanVariant) =>
    (colorField === "coloc_pp4" ? v.coloc_pp4 : v.r2_to_lead) ?? null;

  // Variants inside the window on the selected chromosome. Widen to the whole
  // chromosome if the fixed window captures fewer than two variants.
  const { local, lo, hi } = useMemo(() => {
    if (!selected) return { local: [] as ManhattanVariant[], lo: 0, hi: 1 };
    const sameChr = variants.filter((v) => v.chr === selected.chr);
    let win = sameChr.filter((v) => Math.abs(v.pos - selected.pos) <= windowBp);
    if (win.length < 2) win = sameChr;
    const ps = win.map((v) => v.pos);
    const rawLo = Math.min(...ps);
    const rawHi = Math.max(...ps);
    const pad = Math.max(5_000, (rawHi - rawLo) * 0.08);
    return { local: win, lo: rawLo - pad, hi: rawHi + pad };
  }, [variants, selected, windowBp]);

  const yMax = useMemo(
    () => Math.max(yMode === "pip" ? 1 : 1, ...local.map((v) => yOf(v))),
    // eslint-disable-next-line react-hooks/exhaustive-deps
    [local, yMode]
  );

  // Gene track: distinct nearest-genes → [minPos, maxPos] of their variants.
  const geneSpans = useMemo(() => {
    const m = new Map<string, { lo: number; hi: number }>();
    for (const v of local) {
      const g = v.nearest_gene;
      if (!g) continue;
      const cur = m.get(g) ?? { lo: v.pos, hi: v.pos };
      cur.lo = Math.min(cur.lo, v.pos);
      cur.hi = Math.max(cur.hi, v.pos);
      m.set(g, cur);
    }
    return [...m.entries()].sort((a, b) => a[1].lo - b[1].lo);
  }, [local]);

  if (!selected || local.length === 0) {
    return (
      <div ref={wrapperRef} className="relative w-full">
        <div className="flex h-[220px] items-center justify-center rounded-lg border border-dashed border-border bg-card/40 text-sm text-muted-foreground">
          Select a variant in the Manhattan to zoom its locus.
        </div>
      </div>
    );
  }

  const mb = (bp: number) => (bp / 1_000_000).toFixed(2);

  return (
    <div ref={wrapperRef} className="relative w-full">
      <ChartFrame
        title={title}
        caption={caption}
        height={height}
        ariaLabel={ariaLabel ?? "Regulatory locus-zoom"}
        margin={{ top: 12, right: 16, bottom: 74, left: 48 }}
      >
        {({ innerWidth, innerHeight }) => {
          const xScale = scaleLinear({ domain: [lo, hi], range: [0, innerWidth] });
          const yScale = scaleLinear({
            domain: [0, yMax * (yMode === "pip" ? 1 : 1.05)],
            range: [innerHeight, 0],
            nice: yMode !== "pip",
          });
          const yTicks = yScale.ticks(5);
          const geneLaneY = innerHeight + 14;
          const axisTop = innerHeight + 34;

          return (
            <>
              <GridRows scale={yScale} ticks={yTicks} width={innerWidth} />

              <MarkGroup>
                {local.map((v) => {
                  const px = xScale(v.pos);
                  const py = yScale(yOf(v));
                  const isSel = v.id === selected.id;
                  const isHover = v.id === hoverId;
                  const cval = colorValOf(v);
                  const fill = cval != null ? sequentialColor(cval) : CONTROL;
                  const hasMotif = (v.motifCount ?? 0) > 0;
                  return (
                    <circle
                      key={v.id}
                      cx={px}
                      cy={py}
                      r={isSel ? 6 : 4}
                      fill={fill}
                      fillOpacity={cval != null ? 1 : 0.55}
                      stroke={
                        isSel
                          ? "var(--color-foreground)"
                          : hasMotif
                            ? MOTIF_HEX
                            : isHover
                              ? "var(--color-foreground)"
                              : "none"
                      }
                      strokeWidth={isSel ? 2 : hasMotif ? 1.4 : isHover ? 1 : 0}
                      className="cursor-pointer"
                      onMouseMove={(e) => {
                        show(e, v);
                        setHoverId(v.id);
                      }}
                      onMouseLeave={() => {
                        hide();
                        setHoverId(null);
                      }}
                      onClick={() => onVariantClick?.(v)}
                    />
                  );
                })}
              </MarkGroup>

              {/* Gene track (neutral annotation lane) */}
              <line
                x1={0}
                x2={innerWidth}
                y1={geneLaneY}
                y2={geneLaneY}
                stroke="var(--color-border)"
                strokeWidth={1}
              />
              {geneSpans.map(([gene, span]) => {
                const x1 = xScale(span.lo);
                const x2 = xScale(span.hi);
                const cxg = (x1 + x2) / 2;
                return (
                  <g key={`gene-${gene}`}>
                    <line
                      x1={x1}
                      x2={x2}
                      y1={geneLaneY}
                      y2={geneLaneY}
                      stroke="var(--color-muted-foreground)"
                      strokeWidth={3}
                      strokeLinecap="round"
                    />
                    <circle cx={cxg} cy={geneLaneY} r={2.2} fill="var(--color-muted-foreground)" />
                    <Text
                      x={cxg}
                      y={geneLaneY + 6}
                      textAnchor="middle"
                      verticalAnchor="start"
                      fontSize={9}
                      fontStyle="italic"
                      fill={AXIS_TEXT}
                    >
                      {gene}
                    </Text>
                  </g>
                );
              })}

              <BottomAxis
                scale={xScale}
                top={axisTop}
                numTicks={5}
                tickFormat={((v: number) => mb(v)) as never}
                label={`chr${chrShort(selected.chr)} position (Mb)`}
              />
              <LeftAxis
                scale={yScale}
                numTicks={5}
                label={yMode === "pip" ? "PIP" : "−log₁₀ p"}
              />
            </>
          );
        }}
      </ChartFrame>

      {tooltip && (
        <ChartTooltip left={tooltip.left} top={tooltip.top}>
          <div className="font-medium">
            {tooltip.data.nearest_gene ? (
              <span className="italic">{tooltip.data.nearest_gene}</span>
            ) : (
              tooltip.data.id
            )}
          </div>
          <div className="font-mono text-[10px]">
            {chrShort(tooltip.data.chr)}:{tooltip.data.pos.toLocaleString()}
          </div>
          {tooltip.data.pip != null && <div>PIP {tooltip.data.pip.toFixed(3)}</div>}
          {tooltip.data.coloc_pp4 != null && (
            <div>PP.H4 {tooltip.data.coloc_pp4.toFixed(2)}</div>
          )}
          {(tooltip.data.motifCount ?? 0) > 0 && (
            <div>{tooltip.data.motifCount} motif(s) disrupted</div>
          )}
        </ChartTooltip>
      )}
    </div>
  );
}
