"use client";

/**
 * Forest plot: per-row point estimate with optional confidence interval,
 * against a vertical zero (null-effect) line. Supports faceting into labelled
 * groups and coloring points by ancestry (colorblind-safe Okabe-Ito).
 *
 * Significance is encoded redundantly: significant rows are filled at full
 * opacity, non-significant rows are hollow (white fill, colored ring) at reduced
 * opacity — never color alone.
 */

import { useMemo } from "react";
import { scaleLinear, scaleBand } from "@visx/scale";
import { Line } from "@visx/shape";
import { Text } from "@visx/text";
import {
  ancestryColor,
  divergingColor,
  CONTROL,
  NONSIG_OPACITY,
} from "@/lib/palette";
import { ChartFrame } from "./chart-frame";
import { useChartTooltip, ChartTooltip } from "./chart-tooltip";
import { BottomAxis, GridCols, AXIS_LINE, AXIS_TEXT, MarkGroup } from "./chart-primitives";

export interface ForestRow {
  label: string;
  estimate: number;
  ciLow?: number;
  ciHigh?: number;
  /** Facet group; rows are grouped under a labelled subheader. */
  group?: string;
  /** Ancestry key for coloring (EUR/AFR/EAS/AMR/SAS). */
  ancestry?: string;
  /** Significant? Controls fill vs hollow + opacity. Defaults to true. */
  sig?: boolean;
}

export interface ForestPlotProps {
  data: ForestRow[];
  /** X axis label. */
  xLabel?: string;
  /** Value at the reference line (default 0). */
  nullValue?: number;
  valueFormat?: (v: number) => string;
  title?: React.ReactNode;
  caption?: React.ReactNode;
  height?: number;
  ariaLabel?: string;
}

type Slot =
  | { kind: "header"; label: string }
  | { kind: "row"; row: ForestRow };

export function ForestPlot({
  data,
  xLabel = "effect estimate",
  nullValue = 0,
  valueFormat = (v) => v.toFixed(2),
  title,
  caption,
  height,
  ariaLabel,
}: ForestPlotProps) {
  const { wrapperRef, tooltip, show, hide } = useChartTooltip<ForestRow>();

  const faceted = useMemo(() => data.some((d) => d.group != null), [data]);

  const slots: Slot[] = useMemo(() => {
    if (!faceted) return data.map((row) => ({ kind: "row", row }) as Slot);
    const groups: string[] = [];
    const byGroup = new Map<string, ForestRow[]>();
    for (const d of data) {
      const g = d.group ?? "";
      if (!byGroup.has(g)) {
        byGroup.set(g, []);
        groups.push(g);
      }
      byGroup.get(g)!.push(d);
    }
    const out: Slot[] = [];
    for (const g of groups) {
      out.push({ kind: "header", label: g });
      for (const row of byGroup.get(g)!) out.push({ kind: "row", row });
    }
    return out;
  }, [data, faceted]);

  const xDomain = useMemo(() => {
    const vals: number[] = [nullValue];
    for (const d of data) {
      vals.push(d.estimate);
      if (d.ciLow != null) vals.push(d.ciLow);
      if (d.ciHigh != null) vals.push(d.ciHigh);
    }
    const lo = Math.min(...vals);
    const hi = Math.max(...vals);
    const pad = (hi - lo || 1) * 0.08;
    return [lo - pad, hi + pad] as [number, number];
  }, [data, nullValue]);

  const resolvedHeight =
    height ?? Math.max(140, slots.length * 24 + 52);
  const rowLabelW = 128;

  return (
    <div ref={wrapperRef} className="relative w-full">
      <ChartFrame
        title={title}
        caption={caption}
        height={resolvedHeight}
        ariaLabel={ariaLabel ?? "Forest plot"}
        margin={{ top: 10, right: 20, bottom: 40, left: rowLabelW }}
      >
        {({ innerWidth, innerHeight }) => {
          const xScale = scaleLinear({ domain: xDomain, range: [0, innerWidth] });
          const yScale = scaleBand({
            domain: slots.map((_, i) => i),
            range: [0, innerHeight],
            padding: 0.28,
          });
          const rowH = yScale.bandwidth();
          const xTicks = xScale.ticks(6);
          const zeroX = xScale(nullValue);

          return (
            <>
              <GridCols scale={xScale} ticks={xTicks} height={innerHeight} />

              {/* Null-effect reference line */}
              <Line
                from={{ x: zeroX, y: 0 }}
                to={{ x: zeroX, y: innerHeight }}
                stroke={AXIS_LINE}
                strokeWidth={1.25}
              />

              <MarkGroup>
              {slots.map((slot, i) => {
                const cy = (yScale(i) ?? 0) + rowH / 2;
                if (slot.kind === "header") {
                  return (
                    <Text
                      key={`h-${i}`}
                      x={-rowLabelW + 6}
                      y={cy}
                      verticalAnchor="middle"
                      textAnchor="start"
                      fontSize={10}
                      fontWeight={600}
                      fill="var(--color-foreground)"
                    >
                      {slot.label}
                    </Text>
                  );
                }
                const { row } = slot;
                const sig = row.sig ?? true;
                const baseColor = row.ancestry
                  ? ancestryColor(row.ancestry)
                  : sig
                  ? divergingColor(row.estimate - nullValue)
                  : CONTROL;
                const cx = xScale(row.estimate);
                const x1 = row.ciLow != null ? xScale(row.ciLow) : cx;
                const x2 = row.ciHigh != null ? xScale(row.ciHigh) : cx;

                return (
                  <g
                    key={`r-${i}`}
                    onMouseMove={(e) => show(e, row)}
                    onMouseLeave={hide}
                  >
                    {/* Row label */}
                    <Text
                      x={faceted ? -rowLabelW + 16 : -8}
                      y={cy}
                      verticalAnchor="middle"
                      textAnchor={faceted ? "start" : "end"}
                      fontSize={10}
                      fill={AXIS_TEXT}
                      width={rowLabelW - (faceted ? 20 : 12)}
                    >
                      {row.label}
                    </Text>
                    {/* CI whisker */}
                    {(row.ciLow != null || row.ciHigh != null) && (
                      <Line
                        from={{ x: x1, y: cy }}
                        to={{ x: x2, y: cy }}
                        stroke={baseColor}
                        strokeWidth={1.5}
                        strokeOpacity={sig ? 1 : NONSIG_OPACITY}
                      />
                    )}
                    {/* Point */}
                    <circle
                      cx={cx}
                      cy={cy}
                      r={4}
                      fill={sig ? baseColor : "var(--color-background)"}
                      fillOpacity={sig ? 1 : NONSIG_OPACITY}
                      stroke={baseColor}
                      strokeWidth={sig ? 0 : 1.5}
                    />
                  </g>
                );
              })}
              </MarkGroup>

              <BottomAxis
                scale={xScale}
                top={innerHeight}
                numTicks={6}
                label={xLabel}
              />
            </>
          );
        }}
      </ChartFrame>

      {tooltip && (
        <ChartTooltip left={tooltip.left} top={tooltip.top}>
          <div className="font-medium">{tooltip.data.label}</div>
          <div>
            {valueFormat(tooltip.data.estimate)}
            {tooltip.data.ciLow != null && tooltip.data.ciHigh != null && (
              <>
                {" "}
                [{valueFormat(tooltip.data.ciLow)},{" "}
                {valueFormat(tooltip.data.ciHigh)}]
              </>
            )}
          </div>
          {tooltip.data.ancestry && <div>{tooltip.data.ancestry}</div>}
        </ChartTooltip>
      )}
    </div>
  );
}
