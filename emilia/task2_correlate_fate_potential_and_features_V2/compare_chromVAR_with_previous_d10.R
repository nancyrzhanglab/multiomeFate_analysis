library(Seurat)
library(multiomeFate)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

# ==============================================================================
# Read data
# ==============================================================================
## Read previous chromVars
load('/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/emiliac/lineage_trace/task_identify_features_predictive_growth_and_lineage_specific/data/chromVar_day10_data.RData')

## Read motif names
motifs <- read.csv('/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/emiliac/lineage_trace/task_identify_features_predictive_growth_and_lineage_specific/data/motif_info.csv')

## Read chromVAR
all_data = multiomeFate::data_loader(which_files = c("chromvar"))

# ==============================================================================
# Wrangle data
# ==============================================================================

chromvar_results_cis <- as.data.frame(chromvar_results_cis)
cell_ids_day10_cis <- colnames(chromvar_results_cis)

cell_ids_to_compare <- intersect(colnames(all_data), cell_ids_day10_cis)

chromvar_results_cis$motif_code <- rownames(chromvar_results_cis)
chromvar_results_cis <- merge(chromvar_results_cis, motifs, by = 'motif_code')
chromvar_results_cis <- chromvar_results_cis[, c(cell_ids_to_compare, "motif_names")]
rownames(chromvar_results_cis) <- chromvar_results_cis$motif_names

chromVAR_day10_CIS <- as.data.frame(all_data[["chromVar.day10_CIS"]]@data)
chromVAR_day10_CIS$motif_names <- rownames(chromVAR_day10_CIS)
chromVAR_day10_CIS <- chromVAR_day10_CIS[, c(cell_ids_to_compare, "motif_names")]

# ==============================================================================
# Compare chromVAR results
# ==============================================================================

result_df <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(result_df) <- c("motif_names", 'rho', 'p_val')

for(m in chromvar_results_cis$motif_name){
  old_chromvar_results <- as.data.frame(t(chromVAR_day10_CIS[m, cell_ids_to_compare]))
  colnames(old_chromvar_results) <- c('old')
  
  new_chromvar_results <- as.data.frame(t(chromvar_results_cis[m, cell_ids_to_compare]))
  colnames(new_chromvar_results) <- c('new')
  
  comp_df <- merge(old_chromvar_results, new_chromvar_results, by = 'row.names')
  res <- stats::cor.test(comp_df$old, comp_df$new,
                         alternative = "two.sided",
                         method = "pearson")
  result_df[nrow(result_df) + 1, ] = c(m, res$estimate, res$p.value)
  
}
quantile(as.numeric(result_df$rho))
write.csv(result_df, '/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/Multiome_fate/out/emilia/task2_correlate_fate_potential_and_features_V2/day10_CIS_chromVar_comp_with_previous.csv', row.names = F)
