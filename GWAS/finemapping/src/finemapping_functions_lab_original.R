
## defines region around lead SNPs and submits finemapping pipeline for each locus
## sumstat requirements:
## sumstat must have columns c(chromosome, position, allele1 , allele2, beta, se) and be saved to path_df_sumstats 
## we assume allele1 is effect allele 
## leadSNP requirements:
## df.leadSNP has columns BP, CHR, locus (locus is BP:CHR, what we use to label loci)
generate_per_locus_data_and_submit_pipeline <- function(df.leadSNP, path_df_sumstats, sumstats_name, N_tot, N_cases, window_mb,  ld_pop){
  window_bp <- window_mb*1000000
  df.leadSNP <- df.leadSNP %>%
    mutate(LOWER = ifelse((BP-window_bp) < 0, 1, BP-window_bp), #if lead snp is too close to chr beginning
           UPPER = BP+window_bp)
  
  dir.create(paste0("output/", sumstats_name))
  dir.create(paste0("output/", sumstats_name, "/", ld_pop, "_", window_mb, "Mb_window"))
  dir.create(paste0("output/", sumstats_name, "/", ld_pop, "_", window_mb, "Mb_window/ss"))
  dir.create(paste0("output/", sumstats_name, "/", ld_pop, "_", window_mb, "Mb_window/ld"))
  
  for (i in 1:nrow(df.leadSNP)) {
    LOCUS <- df.leadSNP$locus[i]
    CHR <- df.leadSNP$CHR[i]
    START <- df.leadSNP$LOWER[i]
    END <- df.leadSNP$UPPER[i]
    job_name <- paste("run_FM_pipeline", sumstats_name, LOCUS, paste0(window_mb, "Mb"), sep = "_")
    system(paste("sbatch -J", job_name,
                 "src/run_finemapping_per_locus.sh", #submits src/run_finemapping_per_locus.R using cluster -> paralellizing by locus
                 sumstats_name, 
                 path_df_sumstats, 
                 ld_pop,
                 LOCUS,
                 CHR,
                 START,
                 END,
                 N_tot,
                 N_cases,
                 window_mb, sep = " "))
  }
}

## returns sumstat limited to region around lead SNP
get_ss_per_locus <- function(df.sumstats, LOCUS, CHR,  START, END){
  #df.sumstats <- fread(path_df_sumstats)
  df_chrom <- df.sumstats %>% 
    filter(chromosome == CHR) %>% 
    arrange(position)
  ROW_START = which(df_chrom[,"position"]>START)[1]
  ROW_END = which(df_chrom[,"position"]<END)[length(which(df_chrom[,"position"]<END))]
  df_chrom[ROW_START:ROW_END,"locus"] = as.character(LOCUS)
  return(df_chrom[ROW_START:ROW_END,])
}


save_ss_per_locus <- function(sumstats_name, df.sumstats, LOCUS, CHR, START, END, window_mb){
  ss <- get_ss_per_locus(df.sumstats, as.character(LOCUS), CHR,  START, END)
  ss_dir <- paste0("output/", sumstats_name, "/", ld_pop, "_", window_mb, "Mb/ss/")
  ss_path <- paste0(ss_dir, paste(sumstats_name, paste0(window_mb, "Mb"), as.character(LOCUS), sep = "_"), ".txt")
  write_delim(ss, ss_path, delim = "\t")
  
}


## returns LD matrix for each locus, and corresponding sum stats, 
## modified so that beta corresponds to alt SNP from bim file 
## matches bim and sumstats using chr:pos:alt:ref
get_ld_per_locus <- function(ss_per_locus, LOCUS, CHR, START, END){
  path_to_LD <- find_LD_block(CHR, START, END)
  print(path_to_LD)
  bim <- data.frame()
  ld <- matrix(ncol = 0, nrow = 0)
  for(file in path_to_LD){
    BIM <- fread(paste0(file, ".bim"))
    LD <- fread(paste0(file, ".ld"))
    bim <- rbind(bim, BIM)
    ld <- bdiag(ld, as.matrix(LD)) #make sure this is correct
  }
  
  colnames(bim) <- c("chr", "rsid", "dk", "pos", "alt", "ref")
  bim$SNP <- paste(bim$chr, bim$pos, bim$alt, bim$ref, sep = ":")
 
  colnames(ld) <- rownames(ld) <- bim$SNP
  
  ld <- ld[!duplicated(bim$SNP), !duplicated(bim$SNP)]
  bim <- bim[!duplicated(bim$SNP)]

  ss_per_locus$SNP1 <- paste(ss_per_locus$chromosome,  #allele1 = alt, no need to change direction of beta
                             ss_per_locus$position, 
                             ss_per_locus$allele1, 
                             ss_per_locus$allele2, sep = ":")
  ss_per_locus$SNP2 <- paste(ss_per_locus$chromosome, 
                             ss_per_locus$position, 
                             ss_per_locus$allele2, 
                             ss_per_locus$allele1, sep = ":")

  bim$ss_matched <- bim$SNP %in% c(ss_per_locus$SNP2, ss_per_locus$SNP1)
  ld_filtered <- ld[bim$ss_matched, bim$ss_matched]
  bim_filtered <- bim[bim$ss_matched,]
  ss_filtered <- rbind(ss_per_locus %>% 
                         filter(SNP1 %in% bim_filtered$SNP) %>% 
                         rename(SNP = SNP1) %>% 
                         select(!c(SNP2)),
                       ss_per_locus %>% 
                         filter(SNP2 %in% bim_filtered$SNP) %>% 
                         mutate(beta = -1*beta)  %>% #converts beta to be in same direction as ld allele
                         rename(SNP = SNP2, allele1 = allele2, allele2 = allele1) %>%  #then we rename alleles
                         select(!c(SNP1)))
  
  ss_filtered <- inner_join(bim_filtered %>% select(SNP), ss_filtered, by = "SNP") 
  if(sum(is.na(ld_filtered)) != 0){
    print("LD mat contains NAs")
    ld_filtered <- as.data.frame(as.matrix(ld_filtered))
    ld_filtered_noNA <- as.data.frame(as.matrix(ld_filtered)[unname(which(!sapply(ld_filtered, anyNA))), unname(which(!sapply(ld_filtered, anyNA)))]) #ugly as hell but seems to work
    to_remove <- colnames(ld_filtered)[!colnames(ld_filtered) %in% colnames(ld_filtered_noNA)]
    ss_filtered_noNA <- ss_filtered %>% filter(!SNP %in% to_remove)
    print(dim(ss_filtered_noNA))
    print(dim(ld_filtered_noNA))
    if(!isTRUE(all.equal(ss_filtered_noNA$SNP, colnames(ld_filtered_noNA)) )){
      stop("Problem! ss rsid do not agree with LD columns")
    }else{
      return(list(ss_filtered_noNA, ld_filtered_noNA))
    }
  }else if(sum(ss_filtered$SNP != bim_filtered$SNP) != 0){
    stop("Problem! variants out of order")
  }else{
    return(list(ss_filtered, as.data.frame(as.matrix(ld_filtered))))
  }
}


## determines which LD block locus falls in
find_LD_block <- function(CHR, START, END){
  # Legacy UKBB-LD path removed (PolyFun is the production LD reference). Set UKBB_LD_DIR to use this legacy path.
  ukbb_ld_dir <- Sys.getenv("UKBB_LD_DIR")
  approx_blocks <- fread(file.path(ukbb_ld_dir, "approx_LD_blocks.txt"))
  num_blocks = nrow(approx_blocks %>% filter(chr == CHR,  start < END, stop > START))
  print(paste0("Number of blocks: ", num_blocks))
  LD_start <- approx_blocks %>% filter(chr == CHR, start < END, stop > START) %>% pull(start)
  LD_stop <- approx_blocks %>% filter(chr == CHR, start < END, stop > START) %>% pull(stop)
  filepath <- paste0(ukbb_ld_dir, "/chr", CHR, "/", LD_start, ".", LD_stop, "/", LD_start, ".", LD_stop)
  return(filepath)
}

## run CARMA for sumstat, ld pairing, does not need sample size info
## returns credible set for rho = 0.95
run_CARMA <- function(sumstat, ld, sumstats_name, ld_pop, window_mb, locus, LDpanel){
  z.list<-list()
  ld.list<-list()
  lambda.list<-list()
  sumstat$Z <- sumstat$beta/sumstat$se
  z.list[[1]]<-sumstat$Z
  
  print(packageVersion("Matrix"))
  print("Before LD conversion:")
  ld.list[[1]]<-as.matrix(ld)
  print("After LD conversion:")
  
  lambda.list[[1]]<-1 
  CARMA.results<-CARMA(z.list, ld.list, lambda.list=lambda.list, outlier.switch=TRUE, rho.index = 0.95)
  print("After CARMA method:")
  
  dir.create(paste0("output/",sumstats_name,"/", ld_pop, "_", 
                    window_mb, "Mb/CARMA/"))
  
  saveRDS(CARMA.results, file = paste0("output/",sumstats_name,"/", ld_pop, "_", 
                                       window_mb, "Mb/CARMA/", sumstats_name, "_", LDpanel, "_", 
                                       window_mb, "Mb_",locus,".rds"))
  print("Outliers")
  print(sumstat$SNP[c(CARMA.results[[1]]$Outliers$Index)])
  print("After Outliers conversion:")

  ###### Posterior inclusion probability (PIP) and credible set (CS)
  sumstat.result = sumstat %>% mutate(PIP = CARMA.results[[1]]$PIPs, CS = 0) 
  
  if(length(CARMA.results[[1]]$`Credible set`[[2]])!=0){
    for(l in 1:length(CARMA.results[[1]]$`Credible set`[[2]])){ 
      sumstat.result$CS[CARMA.results[[1]]$`Credible set`[[2]][[l]]]=l
    } }
  
  sumstat.result$outlier <- 0
  sumstat.result$outlier[CARMA.results[[1]]$Outliers$Index] <- 1

  print("Before writing file conversion:")
  
  ###### write the GWAS summary statistics with PIP and CS
  fwrite(x = sumstat.result,
         file = paste0("output/",sumstats_name,"/", ld_pop, "_", 
                       window_mb, "Mb/CARMA/", sumstats_name, "_", LDpanel, "_", 
                       window_mb, "Mb_",locus,".txt.gz"), 
         sep = "\t", quote = F, na = "NA", row.names = F, col.names = T, compress = "gzip")
  
}

## run SuSiE for sumstat, ld pairing, using modified variance as suggested for working with binary traits
run_susie <- function(sumstats, LDmat, N_tot, N_cases, sumstats_name, ld_pop,  window_mb, locus, LDpanel, coverage = 0.95){
  window_bp <- window_mb*1000000
  if(length(sumstats$rsid) != nrow(LDmat)){
    print("Problem! SS dimensions do not agree with LD matrix dimensions")
  }
  colnames(LDmat) <- rownames(LDmat) <- paste(sumstats$chromosome,sumstats$position,sumstats$allele1, sumstats$allele2, sep = ":")
  
  dir.create(paste0("output/",sumstats_name,"/",ld_pop,"_", 
                    window_mb, "Mb/susie"))
  if(is.na(N_cases)){
    susie_output <- susie_rss(R = as.matrix(LDmat), n = N_tot, bhat = sumstats$beta, shat = sumstats$se)
  }else{
    phi <- N_cases/N_tot
    susie_output <- susie_rss(R = as.matrix(LDmat), n = N_tot, vary_y = 1/(phi*(1-phi)), bhat = sumstats$beta, shat = sumstats$se)
  }

  saveRDS(susie_output, file = paste0("output/",sumstats_name,"/", ld_pop, "_", 
                                       window_mb, "Mb/susie/", sumstats_name, "_", LDpanel, "_", 
                                       window_mb, "Mb_",locus,".rds"))
  
  if(susie_output$converged == TRUE){
    print("converged!")
    sumstats$converged = susie_output$converged
    sumstats$PIP = susie_output$pip
    cs_out <- susie_output$sets
    if(!is.null(cs_out$cs)){
      print(locus)
      print(cs_out$cs)
      sumstats[,paste("CS")] <- 0
      for (j in 1:length(cs_out$cs)) {
        sumstats[cs_out$cs[[j]], paste("CS")] <- j
      }
    }else{
      print("NUll!")
    }
    ## write finemapped sumstats file
    write_tsv(sumstats, paste0("output/",sumstats_name,"/",ld_pop,"_", 
                               window_mb, "Mb/susie/", LDpanel,"_", locus, "_cov", coverage,"_adjustedvar.tsv"))
  }
}


bring_in_UKBB_data <- function(sumstats_name, ld_pop, window_mb, locus_dot){
  #bring in ss inner join from hail

  ss_matched <- fread(cmd = paste0("zcat output/", sumstats_name, "/", ld_pop, "_", window_mb,
                                  "Mb/ld/", sumstats_name, "_",
                                  window_mb,"Mb_", locus_dot,"_matched.tsv.bgz"), header = T)
  
  ss_matched <- ss_matched %>% separate(alleles, into = c("Hail_1", "Hail_2"), sep = ",") %>%
    mutate(Hail_1 = str_replace_all(Hail_1, regex("\\W+"), ""), 
           Hail_2 = str_replace_all(Hail_2, regex("\\W+"), "")) 
  
  ss_matched <- ss_matched %>% mutate(beta1 = ifelse(Hail_2 == allele1, beta, 
                                                     ifelse(Hail_2 == allele2, -1*beta, NA)))
  #check that all are well matched
  sum(is.na(ss_matched$beta1))
  ss_matched <- ss_matched %>% mutate(beta = beta1) %>% select(!beta1)
  
  #now all betas are in terms of Hail_2 allele
  #write new ss:
  ss_matched <- ss_matched %>% 
                select(rsid, chromosome, position, allele1 = Hail_2, allele2 = Hail_1, beta, se, locus)

  
  LD <- read.table(paste0("output/", sumstats_name, "/", ld_pop, "_", window_mb,
                          "Mb/ld/UKBB_LDmat_",
                          sumstats_name, "_", window_mb,"Mb_", locus_dot,".bgz"))
  LD <- as.matrix(LD)
  #make full matrix
  LD[lower.tri(LD)]  <- t(LD)[lower.tri(LD)]
  return(list(ss_matched, LD))
}