"use client";

import { GeneView } from "./gene-view";
import { EmptyState } from "@/components/states";
import { useHashSearchParams } from "@/lib/hash-router";

// Single client route for every gene: #/gene?symbol=SYM. Deep links resolve
// client-side from the URL hash, so any symbol works on a static host with no
// 404 risk (the request always hits `/` -> index.html).
export default function GenePage() {
  const params = useHashSearchParams();
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
