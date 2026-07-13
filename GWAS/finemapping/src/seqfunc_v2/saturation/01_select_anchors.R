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
v <- v[tier==1L & primary_eligible==TRUE & pip>=0.05]
v[,`:=`(effect_allele=toupper(effect_allele),other_allele=toupper(other_allele))]
v <- v[grepl("^[ACGT]$",effect_allele) & grepl("^[ACGT]$",other_allele) & effect_allele!=other_allele]
if(!nrow(v)) stopf("no Tier-1 reliable biallelic SNVs with PIP >= 0.05")

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
hu_parts <- tstrsplit(hu_element,"[:-]")
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
anchors <- anchors[seq_len(min(.N,N))]; anchors[,anchor_rank:=seq_len(.N)]
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
write_json(list(status="PASS",input_reliable_tier1_pip005=nrow(v),single_lifted=nrow(lifted)+nrow(reject),
  allele_rejected=nrow(reject),calibration_excluded=sum(lifted$calibration_excluded),
  coding_excluded=sum(lifted$overlaps_gencode_cds),eligible_unique_alleles=nrow(unique(eligible,by="variant_id_hg38")),
  anchors_selected=nrow(anchors),maximum_anchors=N,build_in="GRCh37",build_out="GRCh38",
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
