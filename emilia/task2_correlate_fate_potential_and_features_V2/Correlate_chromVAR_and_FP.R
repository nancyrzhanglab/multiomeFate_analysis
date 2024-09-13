library(Seurat)
library(multiomeFate)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

# ==============================================================================
# Read data
# ==============================================================================
## Read chromVAR and FP
all_data = multiomeFate::data_loader(which_files = c("chromvar", "fatepotential"))

# ==============================================================================
# Correlate
# ==============================================================================

# correlating d0-d10 DABTRAM with d0 ChromVAR
cur_time = 'd0'
fut_time = 'd10'
treatment = 'DABTRAM'
fp_name = paste0("fatepotential_", treatment, '_', cur_time, '_', fut_time)

fp <- as.data.frame(all_data@misc[[fp_name]][["cell_imputed_score"]])
fp$cell_id <- rownames(fp)
colnames(fp) <- c(fp_name, 'cell_id')

chromVAR_cur <- t(all_data[["chromVar.day0"]]@data)
chromVAR_cur <- chromVAR_cur[fp$cell_id, ]

dabtram_d0_chromVAR_cor_vec <- sapply(1:ncol(chromVAR_cur), function(j){
  res <- stats::cor.test(fp[[fp_name]], chromVAR_cur[,j],
                         alternative = "two.sided",
                         method = "pearson")
  c(res$estimate, res$p.value)
})

dabtram_d0_chromVAR_cor_vec <- as.data.frame(t(dabtram_d0_chromVAR_cor_vec))
colnames(dabtram_d0_chromVAR_cor_vec) <- c("correlation", "p.value")
rownames(dabtram_d0_chromVAR_cor_vec) <- colnames(chromVAR_cur)

# correlating d0-d10 COCL2 with d0 ChromVAR
cur_time = 'd0'
fut_time = 'd10'
treatment = 'COCL2'
fp_name = paste0("fatepotential_", treatment, '_', cur_time, '_', fut_time)

fp <- as.data.frame(all_data@misc[[fp_name]][["cell_imputed_score"]])
fp$cell_id <- rownames(fp)
colnames(fp) <- c(fp_name, 'cell_id')

chromVAR_cur <- t(all_data[["chromVar.day0"]]@data)
chromVAR_cur <- chromVAR_cur[fp$cell_id, ]

cocl2_d0_chromVAR_cor_vec <- sapply(1:ncol(chromVAR_cur), function(j){
  res <- stats::cor.test(fp[[fp_name]], chromVAR_cur[,j],
                         alternative = "two.sided",
                         method = "pearson")
  c(res$estimate, res$p.value)
})

cocl2_d0_chromVAR_cor_vec <- as.data.frame(t(cocl2_d0_chromVAR_cor_vec))
colnames(cocl2_d0_chromVAR_cor_vec) <- c("correlation", "p.value")
rownames(cocl2_d0_chromVAR_cor_vec) <- colnames(chromVAR_cur)

# correlating d0-d10 CIS with d0 ChromVAR
cur_time = 'd0'
fut_time = 'd10'
treatment = 'CIS'
fp_name = paste0("fatepotential_", treatment, '_', cur_time, '_', fut_time)

fp <- as.data.frame(all_data@misc[[fp_name]][["cell_imputed_score"]])
fp$cell_id <- rownames(fp)
colnames(fp) <- c(fp_name, 'cell_id')

chromVAR_cur <- t(all_data[["chromVar.day0"]]@data)
chromVAR_cur <- chromVAR_cur[fp$cell_id, ]

cis_d0_chromVAR_cor_vec <- sapply(1:ncol(chromVAR_cur), function(j){
  res <- stats::cor.test(fp[[fp_name]], chromVAR_cur[,j],
                         alternative = "two.sided",
                         method = "pearson")
  c(res$estimate, res$p.value)
})

cis_d0_chromVAR_cor_vec <- as.data.frame(t(cis_d0_chromVAR_cor_vec))
colnames(cis_d0_chromVAR_cor_vec) <- c("correlation", "p.value")
rownames(cis_d0_chromVAR_cor_vec) <- colnames(chromVAR_cur)

# correlating d10-wk5 DABTRAM with d10 ChromVAR DABTRAM
cur_time = 'd10'
fut_time = 'w5'
treatment = 'DABTRAM'
fp_name = paste0("fatepotential_", treatment, '_', cur_time, '_', fut_time)

fp <- as.data.frame(all_data@misc[[fp_name]][["cell_imputed_score"]])
fp$cell_id <- rownames(fp)
colnames(fp) <- c(fp_name, 'cell_id')

chromVAR_cur <- t(all_data[[paste0("chromVar.day10_", treatment)]]@data)
chromVAR_cur <- chromVAR_cur[fp$cell_id, ]

dabtram_d10_chromVAR_cor_vec <- sapply(1:ncol(chromVAR_cur), function(j){
  res <- stats::cor.test(fp[[fp_name]], chromVAR_cur[,j],
                         alternative = "two.sided",
                         method = "pearson")
  c(res$estimate, res$p.value)
})

dabtram_d10_chromVAR_cor_vec <- as.data.frame(t(dabtram_d10_chromVAR_cor_vec))
colnames(dabtram_d10_chromVAR_cor_vec) <- c("correlation", "p.value")
rownames(dabtram_d10_chromVAR_cor_vec) <- colnames(chromVAR_cur)

# correlating d10-wk5 COCL2 with d10 ChromVAR COCL2
cur_time = 'd10'
fut_time = 'w5'
treatment = 'COCL2'
fp_name = paste0("fatepotential_", treatment, '_', cur_time, '_', fut_time)

fp <- as.data.frame(all_data@misc[[fp_name]][["cell_imputed_score"]])
fp$cell_id <- rownames(fp)
colnames(fp) <- c(fp_name, 'cell_id')

chromVAR_cur <- t(all_data[[paste0("chromVar.day10_", treatment)]]@data)
chromVAR_cur <- chromVAR_cur[fp$cell_id, ]

cocl2_d10_chromVAR_cor_vec <- sapply(1:ncol(chromVAR_cur), function(j){
  res <- stats::cor.test(fp[[fp_name]], chromVAR_cur[,j],
                         alternative = "two.sided",
                         method = "pearson")
  c(res$estimate, res$p.value)
})

cocl2_d10_chromVAR_cor_vec <- as.data.frame(t(cocl2_d10_chromVAR_cor_vec))
colnames(cocl2_d10_chromVAR_cor_vec) <- c("correlation", "p.value")
rownames(cocl2_d10_chromVAR_cor_vec) <- colnames(chromVAR_cur)

# correlating d10-wk5 CIS with d10 ChromVAR CIS
cur_time = 'd10'
fut_time = 'w5'
treatment = 'CIS'
fp_name = paste0("fatepotential_", treatment, '_', cur_time, '_', fut_time)

fp <- as.data.frame(all_data@misc[[fp_name]][["cell_imputed_score"]])
fp$cell_id <- rownames(fp)
colnames(fp) <- c(fp_name, 'cell_id')

chromVAR_cur <- t(all_data[[paste0("chromVar.day10_", treatment)]]@data)
chromVAR_cur <- chromVAR_cur[fp$cell_id, ]

cis_d10_chromVAR_cor_vec <- sapply(1:ncol(chromVAR_cur), function(j){
  res <- stats::cor.test(fp[[fp_name]], chromVAR_cur[,j],
                         alternative = "two.sided",
                         method = "pearson")
  c(res$estimate, res$p.value)
})

cis_d10_chromVAR_cor_vec <- as.data.frame(t(cis_d10_chromVAR_cor_vec))
colnames(cis_d10_chromVAR_cor_vec) <- c("correlation", "p.value")
rownames(cis_d10_chromVAR_cor_vec) <- colnames(chromVAR_cur)

# ==============================================================================
# Assemble results
# ==============================================================================
chromVAR_cor_vec <- list(
  dabtram_d0_chromVAR_cor_vec,
  cocl2_d0_chromVAR_cor_vec,
  cis_d0_chromVAR_cor_vec,
  dabtram_d10_chromVAR_cor_vec,
  cocl2_d10_chromVAR_cor_vec,
  cis_d10_chromVAR_cor_vec
)
names(chromVAR_cor_vec) <- c(
  'dabtram_d0_chromVAR_cor_vec',
  'cocl2_d0_chromVAR_cor_vec',
  'cis_d0_chromVAR_cor_vec',
  'dabtram_d10_chromVAR_cor_vec',
  'cocl2_d10_chromVAR_cor_vec',
  'cis_d10_chromVAR_cor_vec'
)

save(date_of_run, session_info,
     chromVAR_cor_vec, 
     file = '/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/Multiome_fate/out/emilia/task2_correlate_fate_potential_and_features_V2/chromVAR_cor_vec.RData')




