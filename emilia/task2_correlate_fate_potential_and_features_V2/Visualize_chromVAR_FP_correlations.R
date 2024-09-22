library(tidyverse)
library(ggplot2)
library(GGally)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

# ==============================================================================
# Read data
# ==============================================================================
load('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/chromVAR_cor_vec.RData')

dabtram_d0_chromVAR_cor_vec <- chromVAR_cor_vec[['dabtram_d0_chromVAR_cor_vec']] 
cocl2_d0_chromVAR_cor_vec <- chromVAR_cor_vec[['cocl2_d0_chromVAR_cor_vec']] 
cis_d0_chromVAR_cor_vec <- chromVAR_cor_vec[['cis_d0_chromVAR_cor_vec']] 

dabtram_d10_chromVAR_cor_vec <- chromVAR_cor_vec[['dabtram_d10_chromVAR_cor_vec']] 
cocl2_d10_chromVAR_cor_vec <- chromVAR_cor_vec[['cocl2_d10_chromVAR_cor_vec']] 
cis_d10_chromVAR_cor_vec <- chromVAR_cor_vec[['cis_d10_chromVAR_cor_vec']] 

# ==============================================================================
# Wrangle data
# ==============================================================================
colnames(dabtram_d0_chromVAR_cor_vec) <- paste0(colnames(dabtram_d0_chromVAR_cor_vec), '.DABTRAM_d0')
colnames(cocl2_d0_chromVAR_cor_vec) <- paste0(colnames(cocl2_d0_chromVAR_cor_vec), '.COCL2_d0')
colnames(cis_d0_chromVAR_cor_vec) <- paste0(colnames(cis_d0_chromVAR_cor_vec), '.CIS_d0')

colnames(dabtram_d10_chromVAR_cor_vec) <- paste0(colnames(dabtram_d10_chromVAR_cor_vec), '.DABTRAM_d10')
colnames(cocl2_d10_chromVAR_cor_vec) <- paste0(colnames(cocl2_d10_chromVAR_cor_vec), '.COCL2_d10')
colnames(cis_d10_chromVAR_cor_vec) <- paste0(colnames(cis_d10_chromVAR_cor_vec), '.CIS_d10')

# Day0
d0_cor <- merge(dabtram_d0_chromVAR_cor_vec, cocl2_d0_chromVAR_cor_vec, by = 'row.names')
rownames(d0_cor) <- d0_cor$Row.names
d0_cor <- d0_cor |> select(-Row.names)

d0_cor <- merge(d0_cor, cis_d0_chromVAR_cor_vec, by = 'row.names')
rownames(d0_cor) <- d0_cor$Row.names
d0_cor <- d0_cor |> select(-Row.names)

# Day10
d10_cor <- merge(dabtram_d10_chromVAR_cor_vec, cocl2_d10_chromVAR_cor_vec, by = 'row.names')
rownames(d10_cor) <- d10_cor$Row.names
d10_cor <- d10_cor |> select(-Row.names)

d10_cor <- merge(d10_cor, cis_d10_chromVAR_cor_vec, by = 'row.names')
rownames(d10_cor) <- d10_cor$Row.names
d10_cor <- d10_cor |> select(-Row.names)

# ==============================================================================
# Plotting
# ==============================================================================

ggpairs(d0_cor[, c('correlation.DABTRAM_d0', 'correlation.COCL2_d0', 'correlation.CIS_d0')],
        mapping = aes(alpha = 0.7),
        lower = list(continuous = wrap("points", color = "#9DBDFF")),
        diag = list(continuous = wrap("densityDiag", fill = "#9DBDFF"))) +
  theme_Publication()

ggsave('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/correlation_ATAC_chromVAR_FP_day0.png', 
       width = 8, height = 8, dpi = 300)

ggpairs(d10_cor[, c('correlation.DABTRAM_d10', 'correlation.COCL2_d10', 'correlation.CIS_d10')],
        mapping = aes(alpha = 0.7),
        lower = list(continuous = wrap("points", color = "#9DBDFF")),
        diag = list(continuous = wrap("densityDiag", fill = "#9DBDFF"))) +
  theme_Publication()
ggsave('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/correlation_ATAC_chromVAR_FP_day10.png', 
       width = 8, height = 8, dpi = 300)

# ==============================================================================
# Save session info
# ==============================================================================
session_info <- devtools::session_info()
