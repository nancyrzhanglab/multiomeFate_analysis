rm(list = ls())

set.seed(123)

library(tidyverse)
library(ggpubr)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'

treatment <- 'CIS'

remove_unassigned_cells <- TRUE

date_of_run <- Sys.time()
session_info <- devtools::session_info()

# =============================================================================
# reading data
# =============================================================================

load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_chromVar_day0.RData'))

all_data[['chromVar.day0']] <- all_data_chromVar_day0

# remove cells with no lineage
if(remove_unassigned_cells) {
  print("Removing cells with no assigned lineage")
  all_data$keep <- !is.na(all_data$assigned_lineage)
  if(any(!all_data$keep)){
    print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
    all_data <- subset(all_data, keep == TRUE)
  }
}


nonAdaptingFP <- readRDS(paste0(out_dir, 'final_fit_d0_d10_Nonadpating_thres_0_', treatment, '.rds'))
nonAdaptingFP <- nonAdaptingFP[["cell_imputed_score"]]

adaptingFP <- readRDS(paste0(out_dir, 'final_fit_d0_d10_adpating_thres_0_', treatment, '.rds'))
adaptingFP <- adaptingFP[["cell_imputed_score"]]

# =============================================================================
# Check
# =============================================================================
adaptingFP <- as.data.frame(adaptingFP)
nonAdaptingFP <- as.data.frame(nonAdaptingFP)

df <- merge(adaptingFP, nonAdaptingFP, by = 'row.names')
df$bias <- 10**(df$adaptingFP) / (10**(df$adaptingFP) + 10**(df$nonAdaptingFP))
hist(df$bias, breaks = 50)

ggplot(df, aes(x = nonAdaptingFP, y = bias)) +
  geom_jitter(width = 0.1) +
  # geom_abline(intercept = 0, slope = 1, col = 'red') +
  stat_cor(method = 'spearman') +
  theme_minimal()

colnames(df)[1] <- 'cell_id' 

# =============================================================================
# Correlate
# =============================================================================
chromVAR_cur <- all_data@assays[["chromVar.day0"]]@data
chromVAR_cur <- as.matrix(chromVAR_cur)
chromVAR_cur <- chromVAR_cur[, df$cell_id]
chromVAR_cur <- t(chromVAR_cur)

chromVAR_cor_vec <- sapply(1:ncol(chromVAR_cur), function(j){
  res <- stats::cor.test(df[['bias']], chromVAR_cur[,j],
                         alternative = "two.sided",
                         method = "pearson")
  c(res$estimate, res$p.value)
})

chromVAR_cor_vec <- as.data.frame(t(chromVAR_cor_vec))
colnames(chromVAR_cor_vec) <- c("correlation", "p.value")
rownames(chromVAR_cor_vec) <- colnames(chromVAR_cur)

save(date_of_run, session_info,
     chromVAR_cor_vec, 
     file = paste0(out_dir, 'TFchromVAR_on_day0_cor_vec_', treatment, '.RData'))
