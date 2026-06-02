import fs from "fs";
import path from "path";
import { GeneEvidenceCard } from "./gene-evidence-card";

export function generateStaticParams() {
  const indexPath = path.join(process.cwd(), "public", "data", "gene_index.json");
  const genes: Array<{ symbol: string }> = JSON.parse(
    fs.readFileSync(indexPath, "utf-8")
  );
  return genes.map((g) => ({ symbol: g.symbol }));
}

export default async function GenePage({
  params,
}: {
  params: Promise<{ symbol: string }>;
}) {
  const { symbol } = await params;
  return <GeneEvidenceCard symbol={symbol} />;
}
