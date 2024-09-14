library(tidyverse)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

# ==============================================================================
# Read data
# ==============================================================================
load('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/chromVAR_cor_vec.RData')

dabtram_d10_chromVAR_cor_vec <- chromVAR_cor_vec[['dabtram_d10_chromVAR_cor_vec']] 
cocl2_d10_chromVAR_cor_vec <- chromVAR_cor_vec[['cocl2_d10_chromVAR_cor_vec']] 
cis_d10_chromVAR_cor_vec <- chromVAR_cor_vec[['cis_d10_chromVAR_cor_vec']] 

anova_in_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task1_ANOVA_lineage_specific_features_V2/'
anova_chromVAR_dabtram_d10 <- read.csv(paste0(anova_in_dir, 'day10_DABTRAM_chromVAR_pvals.csv'))
anova_chromVAR_cocl2_d10 <- read.csv(paste0(anova_in_dir, 'day10_COCL2_chromVAR_pvals.csv'))
anova_chromVAR_cis_d10 <- read.csv(paste0(anova_in_dir, 'day10_CIS_chromVAR_pvals.csv'))

# ==============================================================================
# Wrangle data
# ==============================================================================

## Helper
# Take top 25% most correlated / anti-correlated chromVAR features with growth
get_top25_correlated_feature <- function(df) {
  df$abs_correlation <- abs(df$correlation)
  df$p_val_adj <- p.adjust(df$p.value, 'BH')
  df <- df[df$p_val_adj < 0.05, ]
  df <- df %>%
    arrange(desc(abs_correlation))
  
  top25 <- as.integer(nrow(df) * 0.25)
  df_top25 <- df[c(1: top25), ]
  df_top25$feature <- row.names(df_top25)
  
  return(df_top25)
}


## Day10 DABTRAM

# Take lineage specific chromVAR features
anova_chromVAR_dabtram_d10$p_val_adj <- p.adjust(anova_chromVAR_dabtram_d10$p_val, 'BH')
anova_chromVAR_dabtram_d10_sig <- anova_chromVAR_dabtram_d10[anova_chromVAR_dabtram_d10$p_val_adj < 0.05, ]

dabtram_d10_chromVAR_cor_vec_top25 <- get_top25_correlated_feature(dabtram_d10_chromVAR_cor_vec)

# intersect with ANOVA results
dabtram_d10_chromVAR_cor_vec_top25 <- dabtram_d10_chromVAR_cor_vec_top25 %>%
  filter(feature %in% anova_chromVAR_dabtram_d10_sig$feature) %>% 
  select(feature)

## Day10 COCL2

# Take lineage specific chromVAR features
anova_chromVAR_cocl2_d10$p_val_adj <- p.adjust(anova_chromVAR_cocl2_d10$p_val, 'BH')
anova_chromVAR_cocl2_d10_sig <- anova_chromVAR_cocl2_d10[anova_chromVAR_cocl2_d10$p_val_adj < 0.05, ]

cocl2_d10_chromVAR_cor_vec_top25 <- get_top25_correlated_feature(cocl2_d10_chromVAR_cor_vec)

# intersect with ANOVA results
cocl2_d10_chromVAR_cor_vec_top25 <- cocl2_d10_chromVAR_cor_vec_top25 %>%
  filter(feature %in% anova_chromVAR_cocl2_d10_sig$feature) %>% 
  select(feature)

## Day10 CIS

# Take lineage specific chromVAR features
anova_chromVAR_cis_d10$p_val_adj <- p.adjust(anova_chromVAR_cis_d10$p_val, 'BH')
anova_chromVAR_cis_d10_sig <- anova_chromVAR_cis_d10[anova_chromVAR_cis_d10$p_val_adj < 0.05, ]

cis_d10_chromVAR_cor_vec_top25 <- get_top25_correlated_feature(cis_d10_chromVAR_cor_vec)

# intersect with ANOVA results
cis_d10_chromVAR_cor_vec_top25 <- cis_d10_chromVAR_cor_vec_top25 %>%
  filter(feature %in% anova_chromVAR_cis_d10_sig$feature) %>% 
  select(feature)

# ==============================================================================
# Save data
# ==============================================================================
lineage_specific_adapation_TFs <- list(
  dabtram_d10_chromVAR_cor_vec_top25 = dabtram_d10_chromVAR_cor_vec_top25,
  cocl2_d10_chromVAR_cor_vec_top25 = cocl2_d10_chromVAR_cor_vec_top25,
  cis_d10_chromVAR_cor_vec_top25 = cis_d10_chromVAR_cor_vec_top25
)

save(date_of_run, session_info,
     lineage_specific_adapation_TFs, file = '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task3_identify_lineage_specific_adapation_features_V2/lineage_specific_adapation_TFs.RData')

