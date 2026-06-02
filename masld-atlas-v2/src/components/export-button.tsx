"use client";

import { useCallback, useState, useEffect, useRef } from "react";

interface ExportButtonProps {
  /** Ref to the SVG element to export (for SVG/PNG) */
  svgRef?: React.RefObject<SVGSVGElement | null>;
  /** CSV data to export (array of objects) */
  csvData?: Record<string, unknown>[];
  /** Filename prefix */
  filename?: string;
}

function downloadBlob(blob: Blob, filename: string) {
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  a.click();
  URL.revokeObjectURL(url);
}

export function ExportButton({
  svgRef,
  csvData,
  filename = "masld-atlas",
}: ExportButtonProps) {
  const [open, setOpen] = useState(false);
  const menuRef = useRef<HTMLDivElement>(null);

  // Close dropdown on outside click
  useEffect(() => {
    if (!open) return;
    function handleClick(e: MouseEvent) {
      if (menuRef.current && !menuRef.current.contains(e.target as Node)) {
        setOpen(false);
      }
    }
    document.addEventListener("mousedown", handleClick);
    return () => document.removeEventListener("mousedown", handleClick);
  }, [open]);

  const exportSVG = useCallback(() => {
    if (!svgRef?.current) return;
    const svgData = new XMLSerializer().serializeToString(svgRef.current);
    const blob = new Blob([svgData], { type: "image/svg+xml" });
    downloadBlob(blob, `${filename}.svg`);
    setOpen(false);
  }, [svgRef, filename]);

  const exportPNG = useCallback(() => {
    if (!svgRef?.current) return;
    const svgEl = svgRef.current;
    const svgData = new XMLSerializer().serializeToString(svgEl);
    const canvas = document.createElement("canvas");
    const ctx = canvas.getContext("2d");
    if (!ctx) return;

    const img = new Image();
    const svgBlob = new Blob([svgData], {
      type: "image/svg+xml;charset=utf-8",
    });
    const url = URL.createObjectURL(svgBlob);

    img.onload = () => {
      // 2x resolution for retina displays
      canvas.width = img.width * 2;
      canvas.height = img.height * 2;
      ctx.scale(2, 2);
      ctx.drawImage(img, 0, 0);
      URL.revokeObjectURL(url);
      canvas.toBlob((blob) => {
        if (blob) downloadBlob(blob, `${filename}.png`);
      });
    };
    img.src = url;
    setOpen(false);
  }, [svgRef, filename]);

  const exportCSV = useCallback(() => {
    if (!csvData || csvData.length === 0) return;
    const headers = Object.keys(csvData[0]);
    const rows = csvData.map((row) =>
      headers.map((h) => JSON.stringify(row[h] ?? "")).join(","),
    );
    const csv = [headers.join(","), ...rows].join("\n");
    downloadBlob(new Blob([csv], { type: "text/csv" }), `${filename}.csv`);
    setOpen(false);
  }, [csvData, filename]);

  const hasSvg = !!svgRef;
  const hasCsv = !!csvData;

  if (!hasSvg && !hasCsv) return null;

  return (
    <div className="relative inline-block" ref={menuRef}>
      <button
        onClick={() => setOpen(!open)}
        className="inline-flex items-center gap-1 rounded-md border px-2 py-1 text-xs text-muted-foreground hover:text-foreground hover:bg-accent transition-colors"
        aria-label="Export"
      >
        <svg
          className="h-3 w-3"
          fill="none"
          viewBox="0 0 24 24"
          stroke="currentColor"
          strokeWidth={2}
        >
          <path
            strokeLinecap="round"
            strokeLinejoin="round"
            d="M4 16v1a3 3 0 003 3h10a3 3 0 003-3v-1m-4-4l-4 4m0 0l-4-4m4 4V4"
          />
        </svg>
        Export
      </button>
      {open && (
        <div className="absolute right-0 top-full z-50 mt-1 min-w-[100px] rounded-md border bg-popover p-1 shadow-md">
          {hasSvg && (
            <button
              onClick={exportSVG}
              className="block w-full rounded px-3 py-1.5 text-left text-xs hover:bg-accent"
            >
              SVG
            </button>
          )}
          {hasSvg && (
            <button
              onClick={exportPNG}
              className="block w-full rounded px-3 py-1.5 text-left text-xs hover:bg-accent"
            >
              PNG @2x
            </button>
          )}
          {hasCsv && (
            <button
              onClick={exportCSV}
              className="block w-full rounded px-3 py-1.5 text-left text-xs hover:bg-accent"
            >
              CSV
            </button>
          )}
        </div>
      )}
    </div>
  );
}
