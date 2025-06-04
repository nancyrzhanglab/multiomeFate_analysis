rm(list = ls())
library(tidyverse)
library(ggplot2)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
data_other_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data_other/Torre_CRISPR_screen_BRAFi/'
result_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/'

# ==============================================================================
# Read data
# ==============================================================================

load(paste0(result_dir, 'saver_cor_vec.RData'))

dabtram_d0_saver_cor_vec <- saver_cor_vec[['dabtram_d0_saver_cor_vec']] 
cocl2_d0_saver_cor_vec <- saver_cor_vec[['cocl2_d0_saver_cor_vec']] 
cis_d0_saver_cor_vec <- saver_cor_vec[['cis_d0_saver_cor_vec']] 


dabtram_d10_saver_cor_vec <- saver_cor_vec[['dabtram_d10_saver_cor_vec']] 
cocl2_d10_saver_cor_vec <- saver_cor_vec[['cocl2_d10_saver_cor_vec']] 
cis_d10_saver_cor_vec <- saver_cor_vec[['cis_d10_saver_cor_vec']] 

crispr_hits <- read.csv(paste0(data_other_dir, 'Screen_hits.csv'))
crispr_hits <- crispr_hits[, c('gene', 'domain', 'library', 'primaryCellularPrimingScreen_medianlFC', 'primaryResistanceScreen_medianlFC')] %>% 
  distinct()

crispr_hits$primaryCellularPrimingScreen_medianlFC <- as.numeric(crispr_hits$primaryCellularPrimingScreen_medianlFC)
crispr_hits$primaryResistanceScreen_medianlFC <- as.numeric(crispr_hits$primaryResistanceScreen_medianlFC)
crispr_hits <- crispr_hits[crispr_hits$library != 'TF', ]
# ==============================================================================
# Format data
# ==============================================================================
colnames(dabtram_d0_saver_cor_vec) <- paste0(colnames(dabtram_d0_saver_cor_vec), '.DABTRAM_d0')
colnames(cocl2_d0_saver_cor_vec) <- paste0(colnames(cocl2_d0_saver_cor_vec), '.COCL2_d0')
colnames(cis_d0_saver_cor_vec) <- paste0(colnames(cis_d0_saver_cor_vec), '.CIS_d0')


colnames(dabtram_d10_saver_cor_vec) <- paste0(colnames(dabtram_d10_saver_cor_vec), '.DABTRAM_d10')
colnames(cocl2_d10_saver_cor_vec) <- paste0(colnames(cocl2_d10_saver_cor_vec), '.COCL2_d10')
colnames(cis_d10_saver_cor_vec) <- paste0(colnames(cis_d10_saver_cor_vec), '.CIS_d10')

# make a master dataframe for Day0
d0_cor.RNA <- merge(dabtram_d0_saver_cor_vec, cocl2_d0_saver_cor_vec, by = 'row.names')
rownames(d0_cor.RNA) <- d0_cor.RNA$Row.names
d0_cor.RNA <- d0_cor.RNA |> select(-Row.names)

d0_cor.RNA <- merge(d0_cor.RNA, cis_d0_saver_cor_vec, by = 'row.names')
rownames(d0_cor.RNA) <- d0_cor.RNA$Row.names
d0_cor.RNA <- d0_cor.RNA |> select(-Row.names)

d0_cor.RNA$gene <- rownames(d0_cor.RNA)
d0_cor.RNA <- merge(d0_cor.RNA, crispr_hits, by = 'gene', all.x = TRUE)

# make a master dataframe for Day10
d10_cor.RNA <- merge(dabtram_d10_saver_cor_vec, cocl2_d10_saver_cor_vec, by = 'row.names')
rownames(d10_cor.RNA) <- d10_cor.RNA$Row.names
d10_cor.RNA <- d10_cor.RNA |> select(-Row.names)

d10_cor.RNA <- merge(d10_cor.RNA, cis_d10_saver_cor_vec, by = 'row.names')
rownames(d10_cor.RNA) <- d10_cor.RNA$Row.names
d10_cor.RNA <- d10_cor.RNA |> select(-Row.names)

d10_cor.RNA$gene <- rownames(d10_cor.RNA)

d10_cor.RNA <- merge(d10_cor.RNA, crispr_hits, by = 'gene', all.x = TRUE)

# ==============================================================================
# Get hits
# ==============================================================================
resistant_hits <- crispr_hits %>% 
  arrange(desc(primaryResistanceScreen_medianlFC)) %>%
  head(100)

priming_hits <- crispr_hits %>% 
  arrange(desc(primaryCellularPrimingScreen_medianlFC)) %>%
  head(100)


# ==============================================================================
# Plot
# ==============================================================================
ggplot(d10_cor.RNA, aes(x = correlation.DABTRAM_d10, y = primaryResistanceScreen_medianlFC)) +
  geom_point() +
  geom_point(data = subset(d10_cor.RNA, gene %in% resistant_hits$gene), aes(x = correlation.DABTRAM_d10, y = primaryResistanceScreen_medianlFC), color = 'red') +
  stat_cor() +
  theme_bw()


ggplot(d0_cor.RNA, aes(x = correlation.DABTRAM_d0, y = primaryCellularPrimingScreen_medianlFC)) +
  geom_point() +
  stat_cor() +
  theme_bw()
