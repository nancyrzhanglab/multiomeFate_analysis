library(tidyr)
library(Seurat)
library(multiomeFate)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

out_dir <- '/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/Multiome_fate/out/emilia/task4_differential_winner_loser_d0_V2/'
# ==============================================================================
# Read data
# ==============================================================================
## Read saver and FP
all_data = multiomeFate::data_loader(which_files = c("saver", "fatepotential"))

## Read features to test
load('/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/Multiome_fate/out/emilia/task3_identify_lineage_specific_adapation_features_V2/lineage_specific_adapation_genes_saver.RData')

# ==============================================================================
# Divide into winner and loser cells
# ==============================================================================

# get fate potential
cur_time = 'd0'
fut_time = 'd10'
treatment = 'DABTRAM'
fp_name = paste0("fatepotential_", treatment, '_', cur_time, '_', fut_time)

fp <- as.data.frame(all_data@misc[[fp_name]][["cell_imputed_score"]])
fp$cell_id <- rownames(fp)
colnames(fp) <- c(fp_name, 'cell_id')

# get winner cells
fp_winner <- fp[fp[[fp_name]] > 0, ]

# get loser cells
fp_loser <- fp[fp[[fp_name]] <= 0, ]

# get saver
metadat <- all_data@meta.data
metadat.day0 <- metadat[metadat$dataset == 'day0', ]
metadat.day0$cell_id <- rownames(metadat.day0)

saver_day0 <- t(all_data[["Saver"]]@data)
saver_day0 <- saver_day0[metadat.day0$cell_id, ]
na_df <- saver_day0[which(is.na(as.data.frame(saver_day0)), arr.ind=TRUE), ]
saveRDS(na_df, file = paste0(out_dir, 'na_df_saver.rds'))

saver_day0_winner <- saver_day0[fp_winner$cell_id, ]
saver_day0_loser <- saver_day0[fp_loser$cell_id, ]

# ==============================================================================
# Differential tests
# ==============================================================================
columns <- c('feature', 'mean_winner', 'mean_other', 't_statistic', 'p_val')
t_test_results <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(t_test_results) <- columns

features <- lineage_specific_adapation_TFs[[paste0(tolower(treatment), "_", fut_time, "_saver_cor_vec_top25")]][[1]]

for (f in features) {
  feature_winner <- saver_day0_winner[, f]
  feature_loser <- saver_day0_loser[, f]
  
  feature_winner <- feature_winner[!is.na(feature_winner)]
  feature_loser <- feature_loser[!is.na(feature_loser)]
  
  res <- t.test(feature_winner, feature_loser, alternative = 'two.sided')
  
  t_statistics <- res[["statistic"]][["t"]]
  t_test_p_val <- res[["p.value"]] 
  
  t_test_results[nrow(t_test_results) + 1, ] <- c(
    f, 
    mean(feature_winner), 
    mean(feature_loser), 
    t_statistics, 
    t_test_p_val
  )
}

dim(t_test_results)

# ==============================================================================
# Save results
# ==============================================================================

save(date_of_run, session_info,
     t_test_results, file = paste0(out_dir, 'differential_winner_loser_saver_lineage_specific_gene_', treatment, 't_test_results.RData'))




