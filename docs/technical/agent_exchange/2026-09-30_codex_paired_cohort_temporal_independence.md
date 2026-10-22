# Temporal independence of the paired RNA–H3 cohorts

2026-09-30. Bounded primary-source procurement and enrollment review. Only public metadata/methods, existing identity receipts and metadata field names were inspected. No RNA/H3 values, clinical outcome values, genotype records, receiving sequence files, author contact, jobs or installations were used. This report is the only file written.

**Result: temporal disjointness is not established.** Penn's receiving cohort was collected August 2014–February 2015; the Vanderbilt paper cites both a 2009-start and a 2018-start trial but supplies no biopsy-date lower bound or trial allocation for all 99 paired training participants. The final paper's supplemental methods also supplied no collection-period or enrollment-allocation passage in the bounded check. An earlier living bariatric biopsy followed by deceased organ donation in 2014–2015 remains compatible with the inspected dates. This is an unresolved possibility, not observed overlap.

## Source-specific known facts

| Question | Direct evidence | What it establishes |
|---|---|---|
| When was receiving Penn Cohort 1 procured? | Primary Study Subjects/Penn Cohort 1 passage: **“between August 2014 and February 2015”**. It separately describes 50 liver donors, 49 deceased and 1 living. [Primary article](https://pmc.ncbi.nlm.nih.gov/articles/PMC6612522/); [official Europe PMC article XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6612522/fullTextXML) | Cohort procurement time, not GEO submission or sequencing date. The paper's other Penn cohort has a different interval and cannot substitute for Cohort 1. |
| Do the 17 paired people belong to that cohort? | The frozen [raw-roster summary](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/paired_h3k27ac_raw_roster_21997742/summary.json) names 17 paired public B identifiers and excludes B50. Rechecked SOFT `individual`/`cohort` fields: all 17 identifiers resolve to Penn Cohort 1. | The paired roster uses the dated cohort. Extra runs/assays do not create additional people. |
| Are those 17 deceased donors? | SOFT has exactly one living-labeled sample title: `B50_living tissue donor_RNA-Seq` (GSM3663220). B50 is RNA-only, not in the paired 17. [Public B50 record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3663220) | **[Inference]** Combining the study's exactly-one-living cohort with the unique B50 living label places the retained 17 among the 49 deceased donors. This is an author-metadata inference, not 17 independently verified death certificates. |
| How were training biopsies collected? | Patients and sample collection: **“at the time of elective bariatric surgery”**. Study approval identifies Vanderbilt IRBs 090657/171845 and NCT00983463/NCT03407833. [Final primary methods](https://insight.jci.org/articles/view/182232) | Living surgical source and trial affiliations. Normal-histology samples are from that surgical cohort. Neither methods section gives an all 99 collection interval or person-to-trial allocation. |
| When do the cited trials start? | [NCT00983463 API](https://clinicaltrials.gov/api/v2/studies/NCT00983463): `startDateStruct.date=2009-10` without a type in that field. [NCT03407833 API](https://clinicaltrials.gov/api/v2/studies/NCT03407833): `2018-02-13`, type `ACTUAL`. Both locate the study at Vanderbilt University Medical Center. | One registry predates Penn procurement; one starts later. Trial start is not an individual biopsy date, and citation of both does not allocate the 99 training people. |

## Additional checks that constrain the date question

The final JCI [supplemental-data PDF](https://insight.jci.org/articles/view/182232/sd/pdf/render/1) was fetched with an 8-MB cap into memory, converted to text through `pdftotext` stdin/stdout, and not saved. Supplemental methods occupy PDF pages 2–19 before supplemental figures begin on page 20. Searches of that methods span for bariatric/wedge biopsy collection, enrollment, recruitment and either trial identifier returned no matching passage. The main paper therefore remains the identified source for patient collection/approval; the supplement did not add an all 99 timing bound. Clinical-characteristic spreadsheet values and supporting-data values were not opened. Absence in this bounded methods search is not proof that no author-held allocation exists.

The already-local primary [preprint XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11071482/fullTextXML) was checked only in Patients and sample collection and supplement descriptions. It also cites both trials with no all 99 biopsy-date lower bound. Its presence and later publication date do not date procurement.

Both frozen GEO family SOFTs contain submission/update date keys. GSE267145 sample-characteristic **key names only** were additionally scanned for date, collection, trial, study or surgery fields; none was found. No clinical characteristic values were used. This specifically avoids confusing `Sample_submission_date` or `Sample_last_update_date` with when the biopsy was taken. The local [GSE267145 identity join](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/gse267145-authoritative-join-21064930/join_evidence.json) documents the training assay join, but an accession/assay join does not establish cross-study person disjointness.

## The exact inference needed, and its missing premise

Let D be the 17 receiving people and T the 99 training people. On the study's deceased-donor roster interpretation, every d in D was deceased by the end of February 2015. Every t in T was alive at their elective bariatric biopsy. If an authoritative source established that **every training biopsy occurred after February 2015**, the same biological person could not belong to both sets, and `D intersect T = empty` would follow from that chronology.

The **missing premise** is the all 99 training-biopsy lower bound. NCT03407833's 2018 start would be useful only with an authoritative allocation demonstrating that all 99 relevant biopsies were prospectively collected within that later study, rather than merely analyzed under its protocol. No such allocation was found. NCT00983463 starts before Penn procurement; the paper does not rule out earlier biopsies from it. IRB numbers, study registration, publication/submission years and analysis dates are not specimen dates.

**[Unverified] Universal biological-person disjointness remains unresolved.** Different centers, donor labels, accessions, clinical source populations and sequencing assays provide source context but cannot replace this missing chronological premise or a direct authoritative no-overlap statement. No overlap is demonstrated, and no overlap probability is estimated. The appropriate present description remains different-source RNA–chromatin transport with qualified participant independence.

The exact metadata that would resolve this bounded question is an authoritative all 99 biopsy-date lower bound after February 2015, a trial allocation explicitly tied to those later biopsies, or an authoritative cross-source no-overlap declaration. This review did not contact source owners or infer identity from molecular data. It does not establish the additional requirements for independent predictive validation, such as the unchanged scorer/target, assay compatibility and receiving performance.

## Inspection receipts and limits

Primary public API/XML/PDF responses were memory-only. Dynamic responses/hashes identify what was inspected, not immutable registry history. The dates were read directly from the official trial JSON; browser rendering alone was insufficient.

| Inspected artifact | Bytes | SHA256 |
|---|---:|---|
| NCT00983463 official API |5,517|`a6fb2d36e24beb6628b4255b36d41fb12f640ce1037c588ada5eaa7872b3a228`|
| NCT03407833 official API |10,710|`19af41547709828609a0d5d5e0f4a0c90636f04d02dc96e5d97f4454add14ef7`|
| Penn primary full-text XML |289,070|`0283da7100fda68b0ce662a48b2b9a0cb06089648e4468fe429ccf158c1165f8`|
| Final JCI supplemental PDF |6,964,033|`0f99e455270a3908727d2cb592d720c867ce58da90dbc56b38d66c25a24ae06b`|
| Already-local GSE128072 SOFT |66,214|`92d47ccb362fec5913e899d765982f413c1b39327ab92b80a40af77536030def`|
| Already-local GSE267145 SOFT |27,272|`4ac26f6e003157919bf19b3efd36f1074c5511a11ca559a33aa6f2d26e06cac5`|

Browser access to the Penn PMC page produced reCAPTCHA, and trial/API rendering was unavailable; public official metadata APIs supplied the stated evidence directly. The final supplemental PDF was inspected through its public article link, without clinical spreadsheets or molecular outcomes. Source enumeration was limited to the two referenced papers, their public methods/supplements, the two cited trials and existing cohort metadata. No new cohort search was performed.

Executed checks: local metadata field/roster census, official trial/date extraction, primary procurement passage inspection, bounded supplemental-methods search and report whitespace check. All completed. Genotype identity testing, source-owner confirmation, individual death/biopsy dates and predictive validation were **not run**. No authority or adoption state changes.
