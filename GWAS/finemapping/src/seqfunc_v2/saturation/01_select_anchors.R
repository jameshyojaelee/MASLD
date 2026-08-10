#!/usr/bin/env Rscript
# Select reliable Tier-1 uniform35 SNVs, lift to GRCh38, orient to the reference,
# exclude calibration labels, and remove variants overlapping GENCODE CDS.
suppressPackageStartupMessages({
  library(data.table); library(jsonlite); library(GenomicRanges)
  library(rtracklayer); library(Rsamtools)
})
ROOT <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
RUN <- Sys.getenv("UNIFORM35_RUN_ROOT", file.path(ROOT,"GWAS/finemapping/runs/uniform35_v2_2026-07-13"))
OUT <- Sys.getenv("SATURATION_V2_ROOT", file.path(ROOT,"GWAS/finemapping/results/seqfunc/haplotype_saturation/v2"))
FASTA <- Sys.getenv("GRCH38_FASTA", "/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa")
CHAIN <- Sys.getenv("HG19_TO_HG38_CHAIN", file.path(ROOT,"data/broadaway_eqtl/hg19ToHg38.over.chain"))
GTF <- Sys.getenv("GENCODE_GTF", "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz")
N <- as.integer(Sys.getenv("SATURATION_N_ANCHORS", "500")); dir.create(file.path(OUT,"anchors"),recursive=TRUE,showWarnings=FALSE)
stopf <- function(...) stop(sprintf(...),call.=FALSE)
for(p in c(FASTA,CHAIN,GTF,file.path(RUN,"aggregate/susie_variant_results.tsv.gz"),file.path(RUN,"aggregate/locus_summary.tsv")))
  if(!file.exists(p) || file.info(p)$size==0) stopf("missing required input: %s",p)
gate <- file.path(OUT,"gates/upstream_gate.json")
if(!file.exists(gate) || !isTRUE(fromJSON(gate)$saturation_submission_allowed)) stopf("upstream gate is not PASS: %s",gate)
uniform_audit <- fromJSON(file.path(RUN,"audit/audit_verdict.json"))

v <- fread(file.path(RUN,"aggregate/susie_variant_results.tsv.gz"))
l <- fread(file.path(RUN,"aggregate/locus_summary.tsv"))
reliable <- unique(l[primary_eligible==TRUE,.(study_name,locus_id)])
v <- reliable[v,on=.(study_name,locus_id),nomatch=0]
# Admit tier 1 (disease NAFLD/PDFF) AND tier 2 (liver-enzyme ALT/AST/GGT).
#
# HISTORY (corrected 2026-08-04): a 2026-07-16 note here claimed the refresh
# yielded ZERO reliable tier-1 loci and that tier-2 was "the only honest reliable
# substrate". That was an artifact, not a finding. The uniform35_v2 run modelled
# 17 of 1,950 loci because palindromic variants were dropped and then the whole
# locus was killed for having dropped them, and because a fixed 1e-3 ridge could
# not factorise this substrate. After those fixes the v3 run completes 1,857/1,950
# loci and yields 43 primary-eligible tier-1 disease loci across 10 MASLD studies.
# The tier widening is retained on its own merits -- liver-enzyme loci are
# legitimate substrate -- but it is no longer a fallback, and the disease-substrate
# gap it was introduced to disclose no longer exists.
#
# The reliability (primary_eligible) and PIP>=0.05 gates are UNCHANGED.
# substrate_phenotype in anchor_contract.json is computed from the data and will
# now report disease_and_enzyme rather than the tier-2-only wording.
v <- v[primary_eligible==TRUE & pip>=0.05 & tier %in% c(1L,2L)]
v[,`:=`(effect_allele=toupper(effect_allele),other_allele=toupper(other_allele))]
v <- v[grepl("^[ACGT]$",effect_allele) & grepl("^[ACGT]$",other_allele) & effect_allele!=other_allele]
if(!nrow(v)) stopf("no reliable biallelic SNVs with PIP >= 0.05 (tier 1 or tier 2)")
n_tier1_reliable <- v[tier==1L,.N]; n_tier2_reliable <- v[tier==2L,.N]

# Single-mapping liftover. Alleles are unoriented genomic alleles until FASTA validation.
gr37 <- GRanges(paste0("chr",v$chromosome),IRanges(v$position,width=1)); mcols(gr37)$row_id <- seq_len(nrow(v))
lo <- liftOver(gr37,import.chain(CHAIN)); keep <- which(lengths(lo)==1L)
gr38 <- unlist(lo[keep]); lifted <- v[mcols(gr38)$row_id]
lifted[,`:=`(chr_hg38=as.character(seqnames(gr38)),pos_hg38=start(gr38))]
valid_chr <- paste0("chr",1:22); lifted <- lifted[chr_hg38 %in% valid_chr]
fa <- FaFile(FASTA); open(fa); on.exit(close(fa),add=TRUE)
obs <- as.character(scanFa(fa,GRanges(lifted$chr_hg38,IRanges(lifted$pos_hg38,width=1))))
comp <- function(x) chartr("ACGT","TGCA",x)
lifted[,fasta_ref:=toupper(obs)]
lifted[,orientation:=fcase(fasta_ref==effect_allele,"effect_is_ref",
                            fasta_ref==other_allele,"other_is_ref",
                            fasta_ref==comp(effect_allele),"complement_effect_is_ref",
                            fasta_ref==comp(other_allele),"complement_other_is_ref",default="mismatch")]
lifted[,ref_hg38:=fasta_ref]
lifted[,alt_hg38:=fcase(orientation=="effect_is_ref",other_allele,
                         orientation=="other_is_ref",effect_allele,
                         orientation=="complement_effect_is_ref",comp(other_allele),
                         orientation=="complement_other_is_ref",comp(effect_allele),default=NA_character_)]
reject <- lifted[orientation=="mismatch" | is.na(alt_hg38)]
lifted <- lifted[orientation!="mismatch" & !is.na(alt_hg38)]

# Hu calibration universe: exclude all released assay rsIDs plus any legacy
# allele-resolved substrate rows. Currin: exclude all source positives/background.
hu_files <- list.files(file.path(ROOT,"GWAS/finemapping/results/seqfunc/mpra_benchmark/v2/source_reproduction"),
  pattern="official_comparison\\.tsv\\.gz$",full.names=TRUE)
if(length(hu_files)!=4L) stopf("expected four Hu official-comparison calibration files; found %d",length(hu_files))
hu_rsid <- unique(unlist(lapply(hu_files,function(p) fread(p,select="rsid")$rsid)))
hu_count_files <- list.files(file.path(ROOT,"GWAS/finemapping/results/seqfunc/mpra_benchmark/v2/source_reproduction"),
  pattern="allele_replicate_counts\\.tsv\\.gz$",full.names=TRUE)
if(length(hu_count_files)!=2L) stopf("expected two complete Hu allele-count tables; found %d",length(hu_count_files))
hu_element <- unique(unlist(lapply(hu_count_files,function(p) fread(p,select="element_id")$element_id)))
# The Hu assay contains genomic oligo intervals "chrN:start-end" plus non-genomic
# negative-control elements ("NegControl_*") that define no window. Keep only the
# genomic-interval elements for the anti-leakage overlap; controls cannot overlap a
# genomic anchor, so dropping them does not weaken the leakage guard.
hu_interval <- hu_element[grepl("^chr[0-9XYM]+:[0-9]+-[0-9]+$",hu_element)]
if(!length(hu_interval)) stopf("no parseable Hu genomic element intervals found")
hu_parts <- tstrsplit(hu_interval,"[:-]")
if(length(hu_parts)!=3L || anyNA(as.integer(hu_parts[[2]])) || anyNA(as.integer(hu_parts[[3]])))
  stopf("could not parse complete Hu element intervals")
hu_gr <- GRanges(hu_parts[[1]],IRanges(as.integer(hu_parts[[2]]),as.integer(hu_parts[[3]])))
legacy <- file.path(ROOT,"GWAS/finemapping/results/seqfunc/mpra_benchmark/mpra_substrate_truth.tsv")
hu_hg19 <- character(); if(file.exists(legacy)){x<-fread(legacy); hu_hg19<-unique(x$variant_id_hg19)}
currin_files <- list.files(file.path(ROOT,"GWAS/finemapping/results/seqfunc/chrombpnet_caqtl/v2_truth"),
  pattern="^currin_.*\\.tsv\\.gz$",full.names=TRUE)
if(length(currin_files)!=3L) stopf("expected three Currin truth/background calibration files; found %d",length(currin_files))
currin_id <- unique(unlist(lapply(currin_files,function(p) fread(p,select="variant_id_hg38")$variant_id_hg38)))
lifted[,variant_id_hg19:=paste(chromosome,position,pmin(effect_allele,other_allele),pmax(effect_allele,other_allele),sep=":")]
lifted[,variant_id_hg38:=paste(chr_hg38,pos_hg38,ref_hg38,alt_hg38,sep=":")]
hu_query <- GRanges(lifted$chr_hg38,IRanges(lifted$pos_hg38,width=1L))
hu_window_rows <- unique(queryHits(findOverlaps(hu_query,hu_gr,ignore.strand=TRUE)))
lifted[,hu_assay_window_excluded:=seq_len(.N)%in%hu_window_rows]
lifted[,calibration_excluded:=hu_assay_window_excluded | (!is.na(rsid)&rsid%in%hu_rsid) |
  variant_id_hg19%in%hu_hg19 | variant_id_hg38%in%currin_id]

# "Noncoding" is operationally defined at the nucleotide level: no overlap
# with a GENCODE v49 CDS feature. Promoters, UTRs, introns and intergenic bases remain.
gt <- import(GTF,format="gtf",feature.type="CDS"); seqlevelsStyle(gt) <- "UCSC"
q <- GRanges(lifted$chr_hg38,IRanges(lifted$pos_hg38,width=1))
coding <- unique(queryHits(findOverlaps(q,gt,ignore.strand=TRUE)))
lifted[,overlaps_gencode_cds:=seq_len(.N)%in%coding]
eligible <- lifted[!calibration_excluded & !overlaps_gencode_cds]

# Deduplicate the physical GRCh38 allele while preserving all supporting
# study/locus rows. Selection uses max uniform PIP only, then deterministic keys.
anchors <- eligible[,.(uniform_max_pip=max(pip),
  n_support_rows=.N,n_studies=uniqueN(study_name),n_loci=uniqueN(paste(study_name,locus_id)),
  studies=paste(sort(unique(study_name)),collapse=";"),
  loci=paste(sort(unique(paste(study_name,locus_id,sep="|"))),collapse=";"),
  ancestries=paste(sort(unique(ancestry)),collapse=";"),
  rsids=paste(sort(unique(na.omit(rsid))),collapse=";"),
  representative_hg19=variant_id_hg19[which.max(pip)]),
  by=.(chr_hg38,pos_hg38,ref_hg38,alt_hg38,variant_id_hg38)]
setorder(anchors,-uniform_max_pip,chr_hg38,pos_hg38,ref_hg38,alt_hg38)
# Truncation to the top N by PIP is a substantive selection decision that was
# made in one unlogged line.  With the v3 fixes the eligible pool is far larger
# than 500, so record exactly what was dropped and where the cut fell.
n_eligible_total <- nrow(anchors)
anchor_truncated <- n_eligible_total > N
pip_cutoff_at_rank_N <- if (anchor_truncated) anchors$uniform_max_pip[[N]] else NA_real_
n_excluded_by_cap <- max(0L, n_eligible_total - N)
# Tier is not a column on the deduplicated allele table -- it collapsed into the
# semicolon-joined `studies` field -- so derive it.  Whether the crown actually
# reaches disease loci is the single most important property of this selection,
# and truncating by PIP rank favours the high-powered liver-enzyme GWAS, so it
# must be recorded rather than assumed.
.tier_tbl <- fread(file.path(ROOT,"GWAS/finemapping/config/gwas_trait_tier.tsv"))
.t1_studies <- .tier_tbl[placement=="main" & tier==1L]$study_name
has_tier1 <- function(s) vapply(s, function(x)
  any(strsplit(x,";",fixed=TRUE)[[1]] %in% .t1_studies), logical(1), USE.NAMES=FALSE)
anchors[, tier1_supported := has_tier1(studies)]
tier_before <- list(tier1_supported=sum(anchors$tier1_supported),
                    tier2_only=sum(!anchors$tier1_supported))
anchors <- anchors[seq_len(min(.N,N))]; anchors[,anchor_rank:=seq_len(.N)]
tier_after <- list(tier1_supported=sum(anchors$tier1_supported),
                   tier2_only=sum(!anchors$tier1_supported))
cat(sprintf("[anchors] tier-1 disease support: %d/%d selected (%.1f%%); %d/%d eligible before the cap\n",
            tier_after$tier1_supported, nrow(anchors),
            100*tier_after$tier1_supported/nrow(anchors),
            tier_before$tier1_supported, n_eligible_total))
if (anchor_truncated) {
  cat(sprintf("[anchors] TRUNCATED %d eligible -> %d (PIP cutoff at rank %d = %.4f; %d excluded)\n",
              n_eligible_total, N, N, pip_cutoff_at_rank_N, n_excluded_by_cap))
}
anchors[,`:=`(window_bp=501L,substitutions_planned=1503L,
  evidence_class="predicted",interpretation="mechanistic_nomination_not_causal_evidence",
  apply_only=TRUE)]
setcolorder(anchors,c("anchor_rank","variant_id_hg38","chr_hg38","pos_hg38","ref_hg38","alt_hg38"))
fwrite(anchors,file.path(OUT,"anchors/anchor_manifest.tsv"),sep="\t")

# Preserve the original fine-mapping unit for linked-PIP aggregation.  PIPs
# from different studies or loci are never added together as though they were
# one posterior.  The downstream integration may sum PIPs only within a
# study+locus+target-gene unit (the expected linked causal-variant count).
support <- merge(
  eligible,
  anchors[,.(anchor_rank,variant_id_hg38)],
  by="variant_id_hg38", all=FALSE, allow.cartesian=TRUE
)
support <- unique(support[,.(
  anchor_rank, variant_id_hg38, study_name, locus_id, trait, tier, ancestry,
  chromosome_hg19=chromosome, position_hg19=position, rsid, pip,
  primary_eligible
)], by=c("variant_id_hg38","study_name","locus_id"))
setorder(support,anchor_rank,study_name,locus_id)
fwrite(support,file.path(OUT,"anchors/anchor_pip_support.tsv"),sep="\t")
fwrite(reject,file.path(OUT,"anchors/liftover_allele_rejections.tsv"),sep="\t")
fwrite(lifted[calibration_excluded | overlaps_gencode_cds,
  .(study_name,locus_id,rsid,variant_id_hg19,variant_id_hg38,hu_assay_window_excluded,calibration_excluded,overlaps_gencode_cds)],
  file.path(OUT,"anchors/eligibility_exclusions.tsv"),sep="\t")
write_json(list(status="PASS",input_reliable_pip005=nrow(v),
  reliable_tier1_disease_pip005=n_tier1_reliable,reliable_tier2_enzyme_pip005=n_tier2_reliable,
  substrate_phenotype=if(n_tier1_reliable==0L)"liver_enzyme_proxy_tier2_only" else "disease_and_enzyme",
  tier1_disease_substrate_gap=if(n_tier1_reliable==0L)
    "ZERO reliable tier-1 (NAFLD/PDFF disease) fine-mapped loci in the uniform35 refresh; all reliable substrate is tier-2 liver-enzyme (BBJ ALT/AST/GGT, East-Asian). Saturation runs on the tier-2 liver-enzyme substrate only; disease-GWAS saturation is NOT substantiated by the honest fine-mapping."
    else "tier-1 disease substrate present",
  single_lifted=nrow(lifted)+nrow(reject),
  allele_rejected=nrow(reject),calibration_excluded=sum(lifted$calibration_excluded),
  coding_excluded=sum(lifted$overlaps_gencode_cds),eligible_unique_alleles=nrow(unique(eligible,by="variant_id_hg38")),
  anchors_selected=nrow(anchors),
  anchor_selection_truncated=anchor_truncated,
  eligible_unique_alleles_total=n_eligible_total,
  pip_cutoff_at_rank_N=pip_cutoff_at_rank_N,
  n_eligible_above_cutoff_excluded=n_excluded_by_cap,
  tier_breakdown_before_cap=tier_before,
  tier_breakdown_selected=tier_after,maximum_anchors=N,build_in="GRCh37",build_out="GRCh38",
  uniform_pipeline_complete=uniform_audit$pipeline_complete,
  uniform_release_gate_pass=uniform_audit$release_gate_pass,
  uniform_release_limitation="global PolyFun coverage gaps retained as a limitation; anchors restricted to primary_eligible reliable loci",
  noncoding_definition="does not overlap a GENCODE v49 CDS feature",
  hu_complete_assay_elements=length(hu_element),
  anchor_pip_support_rows=nrow(support),
  linked_pip_aggregation_unit="study+locus+target_gene; never summed across studies/loci",
  calibration_policy="exclude every anchor within a complete Hu assayed oligo interval, released Hu assay rsIDs/allele-resolved legacy rows, and all Currin source truth/background variants",
  pip_source="uniform-prior SuSiE only",canonical_outputs_mutated=FALSE,apply_only_firewall=TRUE),
  file.path(OUT,"anchors/anchor_contract.json"),pretty=TRUE,auto_unbox=TRUE)
message("Selected ",nrow(anchors)," saturation anchors")
