"use client";

/**
 * Fine-mapping PIP Manhattan.
 *
 * Honest to the shipped client data: y = posterior inclusion probability (PIP)
 * of each fine-mapped variant, x = genomic position laid out per chromosome,
 * color = colocalization PP.H4 (sequential ramp). Lead variants (top PIP per
 * chromosome) render at a larger radius; motif-disrupting variants get a stroke
 * ring AND a tick in the marginal rug strip below the axis (a rug, never a
 * lollipop). Clicking a point selects it (brush-links the locus-zoom).
 *
 * SCHEMA-FORWARD: the variant type also accepts `neglog10p` and `r2_to_lead`,
 * so when a genome-wide per-SNP export ships, switching `yField`/`colorField`
 * upgrades this to a classic GWAS Manhattan with no structural rewrite.
 */

import { useMemo, useState } from "react";
import { scaleLinear } from "@visx/scale";
import { Text } from "@visx/text";
import { sequentialColor, CONTROL, MODALITY_HEX } from "@/lib/palette";
import { ChartFrame } from "./chart-frame";
import { useChartTooltip, ChartTooltip } from "./chart-tooltip";
import { LeftAxis, GridRows, AXIS_LINE, AXIS_TEXT, MarkGroup } from "./chart-primitives";

export interface ManhattanVariant {
  id: string;
  chr: string;
  pos: number;
  /** Fine-mapping posterior inclusion probability (current y). */
  pip?: number;
  /** Genome-wide −log10 p (future y; schema-forward). */
  neglog10p?: number;
  /** Colocalization PP.H4 (current color). */
  coloc_pp4?: number | null;
  /** LD r² to the lead SNP (future color; schema-forward). */
  r2_to_lead?: number | null;
  nearest_gene?: string;
  /** Count of disrupted TF motifs (ring + rug tick when > 0). */
  motifCount?: number;
}

export interface ManhattanProps {
  variants: ManhattanVariant[];
  /** Which field drives the y axis. Auto-detects PIP vs −log10 p if unset. */
  yField?: "pip" | "neglog10p";
  /** Which field drives the point color. Default colocalization PP.H4. */
  colorField?: "coloc_pp4" | "r2_to_lead";
  selectedId?: string | null;
  onVariantClick?: (v: ManhattanVariant) => void;
  title?: React.ReactNode;
  caption?: React.ReactNode;
  height?: number;
  ariaLabel?: string;
}

const MOTIF_HEX = MODALITY_HEX.s4_epigenomic; // regulatory / epigenomic hue

/** Numeric sort key for a chromosome label ("chr12" / "12" / "X" / "Y"). */
function chrRank(chr: string): number {
  const c = chr.replace(/^chr/i, "").toUpperCase();
  if (c === "X") return 23;
  if (c === "Y") return 24;
  if (c === "M" || c === "MT") return 25;
  const n = parseInt(c, 10);
  return Number.isNaN(n) ? 99 : n;
}

function chrShort(chr: string): string {
  return chr.replace(/^chr/i, "");
}

export function Manhattan({
  variants,
  yField,
  colorField = "coloc_pp4",
  selectedId,
  onVariantClick,
  title,
  caption,
  height = 320,
  ariaLabel,
}: ManhattanProps) {
  const { wrapperRef, tooltip, show, hide } = useChartTooltip<ManhattanVariant>();
  const [hoverId, setHoverId] = useState<string | null>(null);

  // Auto-pick the y field: PIP if any variant has one, else −log10 p.
  const yMode: "pip" | "neglog10p" = useMemo(() => {
    if (yField) return yField;
    return variants.some((v) => v.pip != null) ? "pip" : "neglog10p";
  }, [yField, variants]);

  const yOf = (v: ManhattanVariant) =>
    (yMode === "pip" ? v.pip : v.neglog10p) ?? 0;
  const colorValOf = (v: ManhattanVariant) =>
    (colorField === "coloc_pp4" ? v.coloc_pp4 : v.r2_to_lead) ?? null;

  // Chromosome band layout (equal-width bands, position-linear within a band).
  const { chrs, bandOf, xInBand, yMax, leadIds } = useMemo(() => {
    const byChr = new Map<string, ManhattanVariant[]>();
    for (const v of variants) {
      const arr = byChr.get(v.chr) ?? [];
      arr.push(v);
      byChr.set(v.chr, arr);
    }
    const chrList = [...byChr.keys()].sort((a, b) => chrRank(a) - chrRank(b));
    const band = new Map<string, number>();
    chrList.forEach((c, i) => band.set(c, i));

    // per-chr [min,max] pos for within-band positioning
    const range = new Map<string, [number, number]>();
    for (const c of chrList) {
      const ps = byChr.get(c)!.map((v) => v.pos);
      range.set(c, [Math.min(...ps), Math.max(...ps)]);
    }
    const xFrac = (v: ManhattanVariant) => {
      const [lo, hi] = range.get(v.chr)!;
      return hi > lo ? (v.pos - lo) / (hi - lo) : 0.5;
    };

    // lead = top-y variant per chromosome
    const leads = new Set<string>();
    for (const c of chrList) {
      const best = byChr
        .get(c)!
        .reduce((m, v) => (yOf(v) > yOf(m) ? v : m));
      leads.add(best.id);
    }

    const ymax = Math.max(
      yMode === "pip" ? 1 : 1,
      ...variants.map((v) => yOf(v))
    );

    return { chrs: chrList, bandOf: band, xInBand: xFrac, yMax: ymax, leadIds: leads };
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [variants, yMode]);

  const nBands = Math.max(1, chrs.length);

  return (
    <div ref={wrapperRef} className="relative w-full">
      <ChartFrame
        title={title}
        caption={caption}
        height={height}
        ariaLabel={ariaLabel ?? "Fine-mapping PIP Manhattan"}
        margin={{ top: 12, right: 14, bottom: 52, left: 48 }}
      >
        {({ innerWidth, innerHeight }) => {
          const bandW = innerWidth / nBands;
          const yScale = scaleLinear({
            domain: [0, yMax * (yMode === "pip" ? 1 : 1.05)],
            range: [innerHeight, 0],
            nice: yMode !== "pip",
          });
          const yTicks = yScale.ticks(5);
          const cx = (v: ManhattanVariant) =>
            (bandOf.get(v.chr) ?? 0) * bandW + (0.12 + 0.76 * xInBand(v)) * bandW;
          const rugY0 = innerHeight + 6;
          const rugY1 = innerHeight + 13;

          return (
            <>
              {/* Alternating chromosome band backgrounds (chrome, not data) */}
              {chrs.map((c, i) =>
                i % 2 === 1 ? (
                  <rect
                    key={`band-${c}`}
                    x={i * bandW}
                    y={0}
                    width={bandW}
                    height={innerHeight}
                    fill="var(--color-muted)"
                    opacity={0.35}
                  />
                ) : null
              )}

              <GridRows scale={yScale} ticks={yTicks} width={innerWidth} />

              {/* Points */}
              <MarkGroup>
                {variants.map((v) => {
                  const px = cx(v);
                  const py = yScale(yOf(v));
                  const isLead = leadIds.has(v.id);
                  const isSel = v.id === selectedId;
                  const isHover = v.id === hoverId;
                  const cval = colorValOf(v);
                  const fill = cval != null ? sequentialColor(cval) : CONTROL;
                  const r = isSel ? 5.5 : isLead ? 4.5 : 2.8;
                  const hasMotif = (v.motifCount ?? 0) > 0;
                  return (
                    <circle
                      key={v.id}
                      cx={px}
                      cy={py}
                      r={r}
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

              {/* Motif-disruption rug strip (ticks, never lollipops) */}
              {variants
                .filter((v) => (v.motifCount ?? 0) > 0)
                .map((v) => (
                  <line
                    key={`rug-${v.id}`}
                    x1={cx(v)}
                    x2={cx(v)}
                    y1={rugY0}
                    y2={rugY1}
                    stroke={MOTIF_HEX}
                    strokeWidth={1.1}
                    strokeOpacity={0.85}
                  />
                ))}

              {/* Chromosome labels */}
              {chrs.map((c, i) => (
                <Text
                  key={`lbl-${c}`}
                  x={i * bandW + bandW / 2}
                  y={innerHeight + 26}
                  textAnchor="middle"
                  verticalAnchor="start"
                  fontSize={9}
                  fill={AXIS_TEXT}
                >
                  {chrShort(c)}
                </Text>
              ))}

              {/* Baseline */}
              <line
                x1={0}
                x2={innerWidth}
                y1={innerHeight}
                y2={innerHeight}
                stroke={AXIS_LINE}
                strokeWidth={1}
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
