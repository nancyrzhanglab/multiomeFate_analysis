library(ggplot2)
library(GGally)

load("~/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/day10_chromVar_day10_growth_potential_for_week5_correlation_writeup6r.RData")

cis_cor_vec <- correlation_list[['cis_cor_vec']]
cocl2_cor_vec <- correlation_list[['cocl2_cor_vec']]
dabtram_cor_vec <- correlation_list[['dabtram_cor_vec']]

cis_cor_vec$motif_name <- rownames(cis_cor_vec)
cocl2_cor_vec$motif_name <- rownames(cocl2_cor_vec)
dabtram_cor_vec$motif_name <- rownames(dabtram_cor_vec)

colnames(cis_cor_vec) <- c( "corr.cis", "p.val.cis", "motif_name")
colnames(cocl2_cor_vec) <- c("corr.cocl2", "p.val.cocl2", "motif_name")
colnames(dabtram_cor_vec) <- c("corr.dabtram","p.val.dabtram", "motif_name")

cor_df <- merge(cis_cor_vec,cocl2_cor_vec, by='motif_name')
cor_df <- merge(cor_df,dabtram_cor_vec, by='motif_name' )
cor_df <- cor_df[, c("motif_name", "corr.cis", "corr.cocl2", "corr.dabtram",
                  "p.val.cis", "p.val.cocl2", "p.val.dabtram")]
cor_df$corr.cis <- as.numeric(cor_df$corr.cis)
cor_df$corr.cocl2 <- as.numeric(cor_df$corr.cocl2)
cor_df$corr.dabtram <- as.numeric(cor_df$corr.dabtram)
ggpairs(cor_df,columns = 2:4, aes(alpha = 0.3))

ggsave('~/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/day10_chromVar_day10_growth_potential_for_week5_correlation_writeup6r.png')
