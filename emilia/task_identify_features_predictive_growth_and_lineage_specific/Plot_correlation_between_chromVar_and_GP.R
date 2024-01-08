library(ggplot2)
library(GGally)

load("~/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/day10_chromVar_day10_growth_potential_for_week5_correlation.RData")

colnames(cis_cor_vec) <- c("motif_names", "corr.cis", "p.val.cis")
colnames(cocl2_cor_vec) <- c("motif_names", "corr.cocl2", "p.val.cocl2")
colnames(dabtram_cor_vec) <- c("motif_names", "corr.dabtram","p.val.dabtram")

cor_df <- merge(cis_cor_vec,cocl2_cor_vec, by='motif_names' )
cor_df <- merge(cor_df,dabtram_cor_vec, by='motif_names' )
cor_df <- cor_df[, c("motif_names", "corr.cis", "corr.cocl2", "corr.dabtram",
                  "p.val.cis", "p.val.cocl2", "p.val.dabtram")]
cor_df$corr.cis <- as.numeric(cor_df$corr.cis)
cor_df$corr.cocl2 <- as.numeric(cor_df$corr.cocl2)
cor_df$corr.dabtram <- as.numeric(cor_df$corr.dabtram)
ggpairs(cor_df,columns = 2:4, aes(alpha = 0.3))

# ggsave('~/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/day10_chromVar_day10_growth_potential_for_week5_correlation.png')
