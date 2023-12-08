library(ggplot2)
library(ggExtra)
library(GGally)

# ==============================================================================
# Read data
# ==============================================================================

# ChromVAR
load("~/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/day10_chromVar_day10_growth_potential_for_week5_correlation.RData")

colnames(cis_cor_vec) <- c("motif_names", "corr.cis.chromVAR", "p.val.cis.chromVAR")
colnames(cocl2_cor_vec) <- c("motif_names", "corr.cocl2.chromVAR", "p.val.cocl2.chromVAR")
colnames(dabtram_cor_vec) <- c("motif_names", "corr.dabtram.chromVAR","p.val.dabtram.chromVAR")

cor_df <- merge(cis_cor_vec,cocl2_cor_vec, by='motif_names' )
cor_df <- merge(cor_df,dabtram_cor_vec, by='motif_names' )
cor_df <- cor_df[, c("motif_names", "corr.cis.chromVAR", "corr.cocl2.chromVAR", "corr.dabtram.chromVAR",
                     "p.val.cis.chromVAR", "p.val.cocl2.chromVAR", "p.val.dabtram.chromVAR")]
cor_df$corr.cis <- as.numeric(cor_df$corr.cis)
cor_df$corr.cocl2 <- as.numeric(cor_df$corr.cocl2)
cor_df$corr.dabtram <- as.numeric(cor_df$corr.dabtram)

cor_df_chromVAR <- cor_df

# RNA
load("~/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/Writeup6n_lineage-imputation_day0-day10-export.RData")

cis_cor_vec <- as.data.frame(correlation_list[["day10_CIS"]])
cocl2_cor_vec <- as.data.frame(correlation_list[["day10_COCL2"]])
dabtram_cor_vec <- as.data.frame(correlation_list[["day10_DABTRAM"]])

cis_cor_vec$gene <- row.names(cis_cor_vec)
cocl2_cor_vec$gene <- row.names(cocl2_cor_vec)
dabtram_cor_vec$gene <- row.names(dabtram_cor_vec)

colnames(cis_cor_vec) <- c("corr.cis.RNA", "p.val.cis.RNA", "gene")
colnames(cocl2_cor_vec) <- c("corr.cocl2.RNA", "p.val.cocl2.RNA", "gene")
colnames(dabtram_cor_vec) <- c("corr.dabtram.RNA","p.val.dabtram.RNA", "gene")

cor_df <- merge(cis_cor_vec,cocl2_cor_vec, by='gene' )
cor_df <- merge(cor_df,dabtram_cor_vec, by='gene' )
cor_df <- cor_df[, c("gene", "corr.cis.RNA", "corr.cocl2.RNA", "corr.dabtram.RNA",
                     "p.val.cis.RNA", "p.val.cocl2.RNA", "p.val.dabtram.RNA")]
cor_df$corr.cis <- as.numeric(cor_df$corr.cis)
cor_df$corr.cocl2 <- as.numeric(cor_df$corr.cocl2)
cor_df$corr.dabtram <- as.numeric(cor_df$corr.dabtram)

cor_df_RNA <- cor_df

# ==============================================================================
# Identifying targets
# ==============================================================================

# CIS and COCL2
p1 <- ggplot(cor_df_chromVAR) +
  geom_point(aes(x = corr.cis.chromVAR, y = corr.cocl2.chromVAR)) +
  geom_hline(yintercept = -0.2, color='red')  +
  geom_vline(xintercept = 0.1, color='red')
ggMarginal(p1, type="density", alpha=0.5)

p2 <- ggplot(cor_df_RNA) +
  geom_point(aes(x = corr.cis.RNA, y = corr.cocl2.RNA))  +
  geom_hline(yintercept = -0.25, color='red')  +
  geom_vline(xintercept = 0.1, color='red')
ggMarginal(p2, type="density", alpha=0.5)

p1 + p2
neg_cocl2_pos_cis_chromVAR <- cor_df_chromVAR[cor_df_chromVAR$corr.cis.chromVAR > 0.1 & cor_df_chromVAR$corr.cocl2.chromVAR < -0.2, ]
neg_cocl2_pos_cis_RNA <- cor_df_RNA[cor_df_RNA$corr.cis.RNA > 0.1 & cor_df_RNA$corr.cocl2.RNA < -0.25, ]


# DABTRAM and COCL2
p3 <- ggplot(cor_df_chromVAR) +
  geom_point(aes(x = corr.cocl2.chromVAR, y = corr.dabtram.chromVAR)) +
  geom_hline(yintercept = 0.25, color='red')  +
  geom_vline(xintercept = -0.2, color='red')
ggMarginal(p3, type="density", alpha=0.5)

p4 <- ggplot(cor_df_RNA) +
  geom_point(aes(x = corr.cocl2.RNA, y = corr.dabtram.RNA), alpha=0.5)  +
  geom_hline(yintercept = 0.11, color='red')  +
  geom_vline(xintercept = -0.25, color='red')
ggMarginal(p4, type="density", alpha=0.5)

neg_cocl2_pos_dabtram_chromVAR <- cor_df_chromVAR[cor_df_chromVAR$corr.dabtram.chromVAR > 0.25 & cor_df_chromVAR$corr.cocl2.chromVAR < -0.2, ]
neg_cocl2_pos_dabtram_RNA <- cor_df_RNA[cor_df_RNA$corr.dabtram.RNA > 0.11 & cor_df_RNA$corr.cocl2.RNA < -0.25, ]


rna_targets <- merge(neg_cocl2_pos_cis_RNA, neg_cocl2_pos_dabtram_RNA,
                     by = 'gene')
tf_targets <- merge(neg_cocl2_pos_dabtram_chromVAR, neg_cocl2_pos_cis_chromVAR,
                    by = 'motif_names')


write.csv(cor_df_RNA, '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/all_genes_all_corrs.csv', row.names = FALSE)
write.csv(rna_targets, '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/common_genes_in_corr_with_growth_v1.csv', row.names = FALSE)
write.csv(tf_targets, '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/common_motifs_in_corr_with_growth_v1.csv', row.names = FALSE)

p2 + p4
ggsave('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/common_genes_in_corr_with_growth_v1.png')

p1 + p3
ggsave('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/common_motifs_in_corr_with_growth_v1.png')

