import type { ReactNode } from "react";
import Link from "next/link";
import { cn } from "@/lib/utils";

/**
 * A gene symbol rendered as an italic monospace link to its gene page.
 *
 * Gene symbols are italicised per publication convention; the link target is
 * the query-param gene route (`/gene?symbol=…`) so it is static-export-safe
 * for arbitrary symbols.
 */
export function GeneLink({
  symbol,
  className,
  children,
}: {
  symbol: string;
  className?: string;
  /** Override the displayed label (defaults to the symbol). */
  children?: ReactNode;
}) {
  return (
    <Link
      href={`/gene?symbol=${encodeURIComponent(symbol)}`}
      className={cn(
        "font-mono text-sm font-medium italic text-primary underline-offset-2 transition-colors hover:underline",
        className
      )}
    >
      {children ?? symbol}
    </Link>
  );
}
