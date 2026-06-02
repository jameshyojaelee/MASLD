/**
 * Type stub for react-force-graph-2d.
 *
 * Install the actual package to get full types:
 *   npm install react-force-graph
 *
 * This declaration provides enough surface area for the network explorer.
 */
declare module "react-force-graph-2d" {
  import { Component } from "react";

  /* eslint-disable @typescript-eslint/no-explicit-any */

  interface GraphData {
    nodes: any[];
    links: any[];
  }

  interface ForceGraph2DProps {
    graphData?: GraphData;
    width?: number;
    height?: number;
    backgroundColor?: string;
    nodeCanvasObject?: (node: any, ctx: CanvasRenderingContext2D, globalScale: number) => void;
    nodePointerAreaPaint?: (node: any, color: string, ctx: CanvasRenderingContext2D) => void;
    linkCanvasObject?: (link: any, ctx: CanvasRenderingContext2D, globalScale: number) => void;
    onNodeClick?: (node: any, event: MouseEvent) => void;
    onNodeHover?: (node: any | null, prevNode: any | null) => void;
    cooldownTicks?: number;
    enableNodeDrag?: boolean;
    nodeVal?: number | ((node: any) => number);
    d3VelocityDecay?: number;
    [key: string]: unknown;
  }

  export default class ForceGraph2D extends Component<ForceGraph2DProps> {
    d3Force(name: string, force?: unknown): any;
    zoomToFit(ms?: number, padding?: number): void;
  }
}
