
suppressPackageStartupMessages(library(dplyr))
scores <- readRDS("/gpfs/commons/groups/sanjana_lab/Cas13/Pathway/results/refactored/all_ssgsea_scores.rds")

# Check InHouse_MCD stats
cat("Summary of InHouse_MCD scores:\n")
scores %>% filter(Dataset == "InHouse_MCD") %>% pull(Score) %>% summary() %>% print()

# Check Ext_MCD_1 stats for comparison
cat("\nSummary of Ext_MCD_1 scores:\n")
scores %>% filter(Dataset == "Ext_MCD_1") %>% pull(Score) %>% summary() %>% print()

# Check a couple of specific pathways
target_pw <- c("HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_FATTY_ACID_METABOLISM")
cat("\nSpecific Pathways (Mean Scores):\n")
scores %>% 
  filter(Pathway %in% target_pw) %>%
  group_by(Dataset, Pathway) %>%
  summarise(Mean_Score = mean(Score), .groups="drop") %>%
  print()
