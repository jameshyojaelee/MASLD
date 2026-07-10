/**
 * Shared chart component library (visx-based).
 *
 * Every portal chart imports colors from `@/lib/palette` (control `#9E9E9E`,
 * diverging effect scale, ancestry, modality, dev-stage, disease-stage) and
 * builds on `ChartFrame` for responsive sizing + consistent chrome. Design
 * contract enforced across all of them:
 *   - control / healthy / non-significant is ALWAYS `#9E9E9E`;
 *   - significance is encoded redundantly (fill/hollow + opacity), never color
 *     alone;
 *   - text uses theme foreground/muted tokens, never a chart hue;
 *   - responsive + light/dark aware; numbers use tabular-nums.
 *
 * Constraints: bars / heatmaps / scatters / strips only — no lollipop, no 3D.
 */

// Frame + shared primitives
export { ChartFrame } from "./chart-frame";
export type {
  ChartFrameProps,
  ChartFrameRenderArgs,
  ChartMargin,
} from "./chart-frame";

export { useChartTooltip, ChartTooltip } from "./chart-tooltip";
export type { TooltipState, ChartTooltipProps } from "./chart-tooltip";

export {
  BottomAxis,
  LeftAxis,
  GridRows,
  GridCols,
  ChartLegend,
  ChartCanvas,
  AXIS_LINE,
  AXIS_TEXT,
  GRID_LINE,
} from "./chart-primitives";
export type { LegendItem } from "./chart-primitives";

// Concrete charts
export { Volcano } from "./volcano";
export type { VolcanoProps, VolcanoPoint } from "./volcano";

export { Heatmap } from "./heatmap";
export type { HeatmapProps, HeatmapCell } from "./heatmap";

export { DotPlot } from "./dotplot";
export type { DotPlotProps, DotPlotPoint } from "./dotplot";

export { ForestPlot } from "./forest-plot";
export type { ForestPlotProps, ForestRow } from "./forest-plot";

export { DivergingBar } from "./diverging-bar";
export type { DivergingBarProps, DivergingBarRow } from "./diverging-bar";

export { Scatter } from "./scatter";
export type { ScatterProps, ScatterPoint, QuadrantLabel } from "./scatter";

export { Bar } from "./bar";
export type { BarProps, BarDatum } from "./bar";

export { StackedBar } from "./stacked-bar";
export type { StackedBarProps, StackedBarDatum } from "./stacked-bar";

export { EvidenceBar } from "./evidence-bar";
export type { EvidenceBarProps } from "./evidence-bar";
