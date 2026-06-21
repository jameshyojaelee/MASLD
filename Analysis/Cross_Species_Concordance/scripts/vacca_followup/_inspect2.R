suppressPackageStartupMessages({library(data.table);library(readxl)})
BASE <- Sys.getenv("MASLD_PROJECT_ROOT","/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
V <- file.path(BASE,"data/external/vacca_2024")
p6 <- as.data.table(read_excel(file.path(V,"42255_2024_1043_MOESM6_ESM.xlsx"),sheet="NES"))
hc <- grep("UCAM|EPoS",names(p6),value=TRUE)
cat("human ref cols:\n"); print(hc)
mc <- setdiff(names(p6)[-1],hc)
cat("\nN model cols:",length(mc),"; first 6:\n"); print(head(mc))
# our mouse fgsea pathways
CS <- file.path(BASE,"Analysis/Cross_Species_Concordance/results")
mr <- fread(file.path(CS,"fgsea_mouse_results.csv"))[grepl("^KEGG_",pathway)]
pnrm <- function(x) toupper(trimws(gsub("_"," ",sub("^KEGG_","",x))))
mr[,pname:=pnrm(pathway)]
v6names <- toupper(trimws(p6[[1]]))
shared <- intersect(unique(mr$pname),v6names)
cat("\nshared KEGG (ours int Vacca):",length(shared),"\n")
metab_kw <- "LIPID|GLUCOS|INSULIN|PPAR|FATTY|GLYCO|CHOLESTEROL|BILE|STEROID|GLUCONEO|PYRUVATE|UNSATURATED|LINOLEIC|KETONE|CARBON METABOLISM|RETINOL|ARACHIDONIC|PROPANOATE|BUTANOATE"
cat("shared that are metabolic:\n")
print(shared[grepl(metab_kw,shared)])
