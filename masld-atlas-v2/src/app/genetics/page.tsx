import { GeneticsClient } from "./genetics-client";

export const metadata = {
  title: "Genetics & Ancestry — MASLD Gene Catalog",
};

// The interactive genetics explorer is a client component (DuckDB-WASM queries
// + visx charts). This server wrapper keeps the route's static metadata.
export default function GeneticsPage() {
  return <GeneticsClient />;
}
