"use client";

import * as React from "react";
import type { ReactNode } from "react";
import { ChevronDown, ChevronUp, ChevronsUpDown } from "lucide-react";
import { cn } from "@/lib/utils";
import { Button } from "@/components/ui/button";
import { GeneLink } from "@/components/gene-link";

export type ColumnAlign = "left" | "center" | "right";
export type SortDir = "asc" | "desc";

export interface DataTableColumn<T> {
  /** Property key on the row (also the value source unless `render`/`sortAccessor` given). */
  key: string;
  header: ReactNode;
  /** Custom cell renderer. Ignored for the `geneColumn`. */
  render?: (row: T) => ReactNode;
  sortable?: boolean;
  align?: ColumnAlign;
  /** Right-aligns and applies tabular numerals. */
  numeric?: boolean;
  /** Value used for sorting (defaults to `row[key]`). */
  sortAccessor?: (row: T) => string | number | null | undefined;
  className?: string;
  /** Fixed column width, e.g. "8rem" or 120. */
  width?: string | number;
}

export interface DataTableProps<T> {
  data: T[];
  columns: DataTableColumn<T>[];
  /** Column key whose cell value is a gene symbol → rendered as a GeneLink. */
  geneColumn?: string;
  /** Rows per page. Pass 0 to render all rows without pagination. Default 25. */
  pageSize?: number;
  initialSort?: { key: string; dir: SortDir };
  rowKey?: (row: T, index: number) => React.Key;
  /** Max height (px or CSS length) for a vertically-scrolling body with sticky header. */
  maxHeight?: number | string;
  emptyMessage?: ReactNode;
  dense?: boolean;
  /** Shade alternating rows for easier row-tracking on wide tables. */
  zebra?: boolean;
  className?: string;
}

function getRaw<T>(row: T, key: string): unknown {
  return (row as Record<string, unknown>)[key];
}

function compare(a: unknown, b: unknown): number {
  const an = a == null || a === "";
  const bn = b == null || b === "";
  if (an && bn) return 0;
  if (an) return 1; // nulls last
  if (bn) return -1;
  if (typeof a === "number" && typeof b === "number") return a - b;
  const na = Number(a);
  const nb = Number(b);
  if (!Number.isNaN(na) && !Number.isNaN(nb)) return na - nb;
  return String(a).localeCompare(String(b));
}

function alignClass(col: DataTableColumn<unknown>): string {
  const align = col.align ?? (col.numeric ? "right" : "left");
  return align === "right"
    ? "text-right"
    : align === "center"
      ? "text-center"
      : "text-left";
}

/**
 * Generic, client-side sortable + paginated table. Click a sortable header to
 * toggle asc/desc. Wide tables scroll horizontally inside their own container
 * so the page never breaks; an optional `maxHeight` adds a sticky-header
 * vertical scroll body.
 */
export function DataTable<T>({
  data,
  columns,
  geneColumn,
  pageSize = 25,
  initialSort,
  rowKey,
  maxHeight,
  emptyMessage = "No rows to display.",
  dense = false,
  zebra = false,
  className,
}: DataTableProps<T>) {
  const [sort, setSort] = React.useState<{ key: string; dir: SortDir } | null>(
    initialSort ?? null
  );
  const [page, setPage] = React.useState(0);

  const sorted = React.useMemo(() => {
    if (!sort) return data;
    const col = columns.find((c) => c.key === sort.key);
    const accessor =
      col?.sortAccessor ?? ((row: T) => getRaw(row, sort.key) as string | number);
    const arr = [...data];
    arr.sort((a, b) => {
      const cmp = compare(accessor(a), accessor(b));
      return sort.dir === "asc" ? cmp : -cmp;
    });
    return arr;
  }, [data, sort, columns]);

  const paginated = pageSize > 0;
  const pageCount = paginated ? Math.max(1, Math.ceil(sorted.length / pageSize)) : 1;
  const safePage = Math.min(page, pageCount - 1);
  const rows = paginated
    ? sorted.slice(safePage * pageSize, safePage * pageSize + pageSize)
    : sorted;

  const handleSort = (col: DataTableColumn<T>) => {
    if (!col.sortable) return;
    setSort((prev) =>
      prev?.key === col.key
        ? { key: col.key, dir: prev.dir === "asc" ? "desc" : "asc" }
        : { key: col.key, dir: "asc" }
    );
    setPage(0);
  };

  const cellPad = dense ? "px-3 py-1.5" : "px-4 py-2";

  return (
    <div className={cn("w-full", className)}>
      <div
        className="overflow-x-auto rounded-lg border border-border"
        style={maxHeight != null ? { maxHeight, overflowY: "auto" } : undefined}
      >
        <table className="w-full text-sm">
          <thead>
            <tr className="border-b border-border bg-muted/50">
              {columns.map((col) => {
                const isSorted = sort?.key === col.key;
                const SortIcon = !col.sortable
                  ? null
                  : !isSorted
                    ? ChevronsUpDown
                    : sort.dir === "asc"
                      ? ChevronUp
                      : ChevronDown;
                return (
                  <th
                    key={col.key}
                    scope="col"
                    aria-sort={
                      isSorted
                        ? sort.dir === "asc"
                          ? "ascending"
                          : "descending"
                        : undefined
                    }
                    style={col.width != null ? { width: col.width } : undefined}
                    className={cn(
                      cellPad,
                      "text-xs font-medium text-muted-foreground",
                      maxHeight != null &&
                        "sticky top-0 z-10 bg-muted shadow-[var(--shadow-elev-1)]",
                      alignClass(col as DataTableColumn<unknown>),
                      col.className
                    )}
                  >
                    {col.sortable ? (
                      <button
                        type="button"
                        onClick={() => handleSort(col)}
                        className={cn(
                          "inline-flex items-center gap-1 transition-colors hover:text-foreground",
                          (col.align ?? (col.numeric ? "right" : "left")) ===
                            "right" && "flex-row-reverse"
                        )}
                      >
                        {col.header}
                        {SortIcon && (
                          <SortIcon
                            className={cn(
                              "size-3",
                              isSorted ? "text-foreground" : "text-muted-foreground/60"
                            )}
                          />
                        )}
                      </button>
                    ) : (
                      col.header
                    )}
                  </th>
                );
              })}
            </tr>
          </thead>
          <tbody>
            {rows.length === 0 ? (
              <tr>
                <td
                  colSpan={columns.length}
                  className={cn(cellPad, "text-center text-sm text-muted-foreground")}
                >
                  {emptyMessage}
                </td>
              </tr>
            ) : (
              rows.map((row, i) => (
                <tr
                  key={rowKey ? rowKey(row, i) : i}
                  className={cn(
                    "border-b border-border/50 transition-colors last:border-0 hover:bg-muted/30",
                    zebra && i % 2 === 1 && "bg-muted/20"
                  )}
                >
                  {columns.map((col) => {
                    const raw = getRaw(row, col.key);
                    let content: ReactNode;
                    if (geneColumn && col.key === geneColumn) {
                      content =
                        raw == null || raw === "" ? (
                          <span className="text-muted-foreground">—</span>
                        ) : (
                          <GeneLink symbol={String(raw)} />
                        );
                    } else if (col.render) {
                      content = col.render(row);
                    } else {
                      content =
                        raw == null || raw === "" ? (
                          <span className="text-muted-foreground">—</span>
                        ) : (
                          String(raw)
                        );
                    }
                    return (
                      <td
                        key={col.key}
                        className={cn(
                          cellPad,
                          "text-sm",
                          col.numeric && "font-numeric",
                          alignClass(col as DataTableColumn<unknown>),
                          col.className
                        )}
                      >
                        {content}
                      </td>
                    );
                  })}
                </tr>
              ))
            )}
          </tbody>
        </table>
      </div>

      {paginated && sorted.length > pageSize && (
        <div className="mt-2 flex items-center justify-between gap-2 text-xs text-muted-foreground">
          <span className="tabular-nums">
            {safePage * pageSize + 1}–
            {Math.min((safePage + 1) * pageSize, sorted.length)} of{" "}
            {sorted.length.toLocaleString("en-US")}
          </span>
          <div className="flex items-center gap-1.5">
            <Button
              variant="outline"
              size="xs"
              disabled={safePage === 0}
              onClick={() => setPage((p) => Math.max(0, p - 1))}
            >
              Previous
            </Button>
            <span className="tabular-nums">
              {safePage + 1} / {pageCount}
            </span>
            <Button
              variant="outline"
              size="xs"
              disabled={safePage >= pageCount - 1}
              onClick={() => setPage((p) => Math.min(pageCount - 1, p + 1))}
            >
              Next
            </Button>
          </div>
        </div>
      )}
    </div>
  );
}
