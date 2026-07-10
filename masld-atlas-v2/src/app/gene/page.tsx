"use client";

import { Suspense } from "react";
import { useSearchParams } from "next/navigation";
import { GeneView } from "./gene-view";
import { EmptyState, SkeletonBlock } from "@/components/states";

// Single client route for every gene: /gene?symbol=SYM. Replaces the 27k
// per-symbol static shells — deep links resolve client-side, so any symbol
// works on a static host with no 404 risk. useSearchParams requires a Suspense
// boundary under `output: export`.
function GeneRouteInner() {
  const params = useSearchParams();
  const symbol = params.get("symbol")?.trim();

  if (!symbol) {
    return (
      <div className="mx-auto max-w-6xl px-6 py-10">
        <EmptyState
          title="No gene selected"
          description="Search for a gene with the command palette (⌘K), or follow a gene link from the atlas."
        />
      </div>
    );
  }

  return <GeneView symbol={symbol} />;
}

/** Shape-matched skeleton mirroring the gene card (identity → stats → tabs). */
function GeneSkeleton() {
  return (
    <div className="mx-auto max-w-6xl space-y-6 px-6 py-8">
      <div className="space-y-2">
        <SkeletonBlock className="h-8 w-48" />
        <SkeletonBlock className="h-4 w-72" />
      </div>
      <div className="grid gap-4 lg:grid-cols-[2fr_1fr]">
        <div className="grid grid-cols-2 gap-3 sm:grid-cols-3">
          {Array.from({ length: 6 }).map((_, i) => (
            <SkeletonBlock key={i} className="h-20" />
          ))}
        </div>
        <SkeletonBlock className="h-44" />
      </div>
      <SkeletonBlock className="h-9 w-full max-w-md" />
      <SkeletonBlock className="h-64" />
    </div>
  );
}

export default function GenePage() {
  return (
    <Suspense fallback={<GeneSkeleton />}>
      <GeneRouteInner />
    </Suspense>
  );
}
