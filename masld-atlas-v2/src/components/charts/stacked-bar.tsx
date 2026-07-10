"use client";

/**
 * Stacked bar chart for composition (e.g. cell-type fractions per disease
 * stage). Long-format input; optional normalisation to 100%. Any category named
 * Healthy/Control is forced to control gray; remaining categories draw from a
 * colorblind-safe categorical palette (override via `colorFor`).
 */

import { useMemo, useState } from "react";
import { scaleLinear, scaleBand } from "@visx/scale";
import { CONTROL, categoricalColor } from "@/lib/palette";
import { ChartFrame } from "./chart-frame";
import { useChartTooltip, ChartTooltip } from "./chart-tooltip";
import { BottomAxis, LeftAxis, GridRows, ChartLegend, MarkGroup } from "./chart-primitives";

export interface StackedBarDatum {
  group: string;
  category: string;
  value: number;
}

export interface StackedBarProps {
  data: StackedBarDatum[];
  groups?: string[];
  categories?: string[];
  /** Normalise each group's segments to sum to 1 (100% composition). */
  normalize?: boolean;
  /** Category → color override. Defaults handle Healthy/Control = gray. */
  colorFor?: (category: string) => string;
  xLabel?: string;
  yLabel?: string;
  showLegend?: boolean;
  title?: React.ReactNode;
  caption?: React.ReactNode;
  height?: number;
  ariaLabel?: string;
}

const isControlCategory = (c: string) =>
  /^(healthy|control|normal)$/i.test(c.trim());

function uniqueInOrder(values: string[]): string[] {
  const seen = new Set<string>();
  const out: string[] = [];
  for (const v of values) {
    if (!seen.has(v)) {
      seen.add(v);
      out.push(v);
    }
  }
  return out;
}

interface Segment {
  group: string;
  category: string;
  value: number;
  y0: number;
  y1: number;
}

export function StackedBar({
  data,
  groups: groupsProp,
  categories: catsProp,
  normalize = false,
  colorFor: colorForProp,
  xLabel,
  yLabel = "count",
  showLegend = true,
  title,
  caption,
  height = 320,
  ariaLabel,
}: StackedBarProps) {
  const { wrapperRef, tooltip, show, hide } = useChartTooltip<Segment>();
  // Legend click-to-isolate: hidden categories are dropped from the stacks.
  const [hidden, setHidden] = useState<Set<string>>(new Set());

  const groups = useMemo(
    () => groupsProp ?? uniqueInOrder(data.map((d) => d.group)),
    [groupsProp, data]
  );
  const categories = useMemo(
    () => catsProp ?? uniqueInOrder(data.map((d) => d.category)),
    [catsProp, data]
  );

  const colorFor = useMemo(() => {
    if (colorForProp) return colorForProp;
    const nonControl = categories.filter((c) => !isControlCategory(c));
    const idx = new Map(nonControl.map((c, i) => [c, i]));
    return (c: string) =>
      isControlCategory(c) ? CONTROL : categoricalColor(idx.get(c) ?? 0);
  }, [colorForProp, categories]);

  // Categories actually stacked = all minus the legend-isolated ones.
  const shownCategories = useMemo(
    () => categories.filter((c) => !hidden.has(c)),
    [categories, hidden]
  );

  const { segments, maxTotal } = useMemo(() => {
    const lookup = new Map<string, number>();
    for (const d of data) lookup.set(`${d.group} ${d.category}`, d.value);
    const segs: Segment[] = [];
    let maxT = 0;
    for (const g of groups) {
      const raw = shownCategories.map((c) => lookup.get(`${g} ${c}`) ?? 0);
      const total = raw.reduce((a, b) => a + b, 0);
      const denom = normalize ? total || 1 : 1;
      let cum = 0;
      shownCategories.forEach((c, i) => {
        const v = raw[i] / denom;
        if (v > 0) {
          segs.push({ group: g, category: c, value: v, y0: cum, y1: cum + v });
        }
        cum += v;
      });
      maxT = Math.max(maxT, cum);
    }
    return { segments: segs, maxTotal: normalize ? 1 : maxT };
  }, [data, groups, shownCategories, normalize]);

  const toggleCategory = (label: string) =>
    setHidden((prev) => {
      const next = new Set(prev);
      // Never let the user hide every series.
      if (next.has(label)) next.delete(label);
      else if (categories.length - next.size > 1) next.add(label);
      return next;
    });

  const yFormat = (v: number) =>
    normalize ? `${Math.round(v * 100)}%` : `${v}`;

  return (
    <div ref={wrapperRef} className="relative w-full">
      <ChartFrame
        title={title}
        caption={caption}
        height={height}
        ariaLabel={ariaLabel ?? "Stacked bar chart"}
        margin={{ top: 10, right: 14, bottom: 40, left: 52 }}
      >
        {({ innerWidth, innerHeight }) => {
          const xScale = scaleBand({
            domain: groups,
            range: [0, innerWidth],
            padding: 0.28,
          });
          const yScale = scaleLinear({
            domain: [0, maxTotal || 1],
            range: [innerHeight, 0],
            nice: !normalize,
          });
          const bw = xScale.bandwidth();
          const yTicks = yScale.ticks(5);

          return (
            <>
              <GridRows scale={yScale} ticks={yTicks} width={innerWidth} />

              <MarkGroup>
                {segments.map((s, i) => {
                  const x = xScale(s.group) ?? 0;
                  const yTop = yScale(s.y1);
                  const yBot = yScale(s.y0);
                  return (
                    <rect
                      key={`${s.group}-${s.category}-${i}`}
                      x={x}
                      y={yTop}
                      width={bw}
                      height={Math.max(0, yBot - yTop)}
                      fill={colorFor(s.category)}
                      stroke="var(--color-background)"
                      strokeWidth={0.75}
                      onMouseMove={(e) => show(e, s)}
                      onMouseLeave={hide}
                    />
                  );
                })}
              </MarkGroup>

              <BottomAxis
                scale={xScale}
                top={innerHeight}
                hideAxisLine={false}
                label={xLabel}
              />
              <LeftAxis
                scale={yScale}
                numTicks={5}
                tickFormat={((v: number) => yFormat(v)) as never}
                label={yLabel}
              />
            </>
          );
        }}
      </ChartFrame>

      {showLegend && (
        <ChartLegend
          items={categories.map((c) => ({ label: c, color: colorFor(c) }))}
          mutedLabels={hidden}
          onItemClick={toggleCategory}
        />
      )}

      {tooltip && (
        <ChartTooltip left={tooltip.left} top={tooltip.top}>
          <div className="font-medium">{tooltip.data.group}</div>
          <div>{tooltip.data.category}</div>
          <div>{yFormat(tooltip.data.value)}</div>
        </ChartTooltip>
      )}
    </div>
  );
}
