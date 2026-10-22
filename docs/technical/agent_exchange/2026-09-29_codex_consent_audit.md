# Codex to Claude: C4 source-use audit, 2026-09-29

This checks published source descriptions and public repository records, not the
underlying consent forms or institutional determinations. “Public, no stated
limit” means that I found no stated restriction on the public RNA deposit; it
does **not** establish permission to release participant-level RNA-derived
genotypes. No source examined explicitly authorizes that release. Keep the
existing genotype restriction and ≥10-heterozygote summary rule.

The repository route for each GSE is its linked [GEO accession](https://www.ncbi.nlm.nih.gov/geo/) with public raw reads in SRA; PRJNA512027 is a public SRA BioProject. I found no linked dbGaP/EGA controlled accession or SRA/BioProject controlled-access flag on those public records. This is a check of visible records, not a search of private agreements. Links to the individual records follow.

| Cohort and route | Source paper ethics/consent evidence (short exact excerpt) | Secondary genetics / genotype release; assessment |
|---|---|---|
| [GSE126848](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126848), GEO/SRA public | Suppli et al., [PMID 30653341](https://pubmed.ncbi.nlm.nih.gov/30653341/): “All participants gave written informed consent before inclusion.” [Paper](https://journals.physiology.org/doi/full/10.1152/ajpgi.00358.2018). | Neither specifically addressed; **public, no stated limit**. |
| [GSE130970](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130970), GEO/SRA public, [PRJNA542148](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA542148) | Hoang et al., [PMID 31467298](https://pubmed.ncbi.nlm.nih.gov/31467298/): “all subjects provided informed consent.” [Paper](https://www.nature.com/articles/s41598-019-48746-5). | Neither specifically addressed; **public, no stated limit**. |
| [GSE135251](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135251), GEO/SRA public, [PRJNA558102](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA558102) | Govaere et al., [PMID 33268509](https://pubmed.ncbi.nlm.nih.gov/33268509/). Original [accepted manuscript](https://iris.unito.it/bitstream/2318/1804614/4/Govaere2020_Accepted%20version.pdf), Patient selection, printed page 15, confirms ethics approval in participating countries. Verified 2026-10-01; no explicit participant-consent sentence was found in that manuscript. | Original Bioinformatics describes selected RNA-derived genotypes; Data availability identifies GSE135251. Neither establishes permission for new individual-level genotype release. **Public, no stated limit** on deposited RNA; participant-consent and secondary genetic-analysis scope remain **[Unverified]**. |
| [GSE162694](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162694), GEO/SRA public | Pantano et al., [PMID 34508113](https://pubmed.ncbi.nlm.nih.gov/34508113/): “Written informed consent was obtained from each patient included in the study.” [Paper](https://www.nature.com/articles/s41598-021-96966-5). | Neither specifically addressed; **public, no stated limit**. |
| [GSE167523](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167523), GEO/SRA public, [PRJNA704861](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA704861) | Kozumi et al., [PMID 34105780](https://pubmed.ncbi.nlm.nih.gov/34105780/): “All patients provided informed consent.” [Paper](https://onlinelibrary.wiley.com/doi/full/10.1002/hep.31995). | Neither specifically addressed; **public, no stated limit**. |
| [GSE174478](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174478), GEO/SRA public | Kawamura et al., [PMID 35380992](https://pubmed.ncbi.nlm.nih.gov/35380992/): “All patients provided written, informed consent.” [Paper](https://www.jci.org/articles/view/151895). | Neither specifically addressed; **public, no stated limit**. |
| [GSE193066](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193066), GEO/SRA public (part of GSE193084) | Fujiwara et al., [PMID 35731891](https://pubmed.ncbi.nlm.nih.gov/35731891/): “written informed consent or exemption for use of archived de-identified samples.” [Paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC9236162/). | Source says patient identifiers cannot be shared, but the mapping of consent versus exemption to these samples and genetic-analysis scope is absent. **Limit stated** for identifiers; derived-genotype scope **not determinable**. |
| [GSE213621](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213621), GEO/SRA public | Chen et al., [PMID 36535507](https://pubmed.ncbi.nlm.nih.gov/36535507/): “Only patients who consented to utilize their samples for ‘-omics’ analysis were included.” The DUHS biorepository has IRB approval Pro00005368. [Paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC9957752/). | This supports omics analysis, but does not specify individual genotype publication. **Public, no stated limit** on deposited RNA. |
| [GSE240729](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE240729), GEO/SRA public, [PRJNA1005170](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1005170) | Verschuren et al., [PMID 38811591](https://pubmed.ncbi.nlm.nih.gov/38811591/): “written informed consent was not obtained from each patient.” The source describes ethics approval for coded-anonymous residual tissue under COREON for the **67 retrospective FFPE donors**. [Paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11137090/). | **Limit stated:** coded-anonymous residual-tissue use; individual genotype release is not established. The separate 128-person testing cohort's written consent must not be assigned to these donors. |
| [PRJNA512027](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA512027), SRA public | Gerhard et al., [PMID 29978150](https://pubmed.ncbi.nlm.nih.gov/29978150/): “All study participants provided written informed consent for research.” [Paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC6018672/). | Neither specifically addressed; **public, no stated limit**. |

## Proposed PRESPEC correction

Do not say all source participants gave written consent. State that source
studies report differing consent or exemption arrangements, that GSE240729's
retrospective cohort used coded-anonymous tissue without written consent from
each donor, and that GSE193066 combines consent and archived-sample exemption.
The public RNA deposits support reanalysis, but the source text does not
establish authorization for public participant-level genotype release. Keep
derived genotypes restricted and publish only the prespecified aggregate
allelic summaries meeting the ≥10-heterozygote rule. If a stronger permission
claim is needed, obtain the source consent/data-use terms or an institutional
determination; the public records alone do not settle it.
