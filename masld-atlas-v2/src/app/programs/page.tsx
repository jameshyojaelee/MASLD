import { PageHeader } from "@/components/page-header";
import { PageContainer } from "@/components/page-container";
import { ProgramsClient } from "./programs-client";

export const metadata = {
  title: "Molecular Programs — MASLD Atlas",
};

export default function ProgramsPage() {
  return (
    <PageContainer>
      <PageHeader
        eyebrow="Decomposition"
        title="Molecular Programs"
        description="NMF and Hotspot molecular programs, program × stage activity, and sex-biased gene sets. The k=6 programs are presented as continuous interpretive axes (P1–P6), not discrete patient subtypes."
      />
      <ProgramsClient />
    </PageContainer>
  );
}
