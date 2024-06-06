library(dbscan)
library(ggplot2)
library(ggExtra)
library(GGally)
library(gridExtra)
library(grid)

# ==============================================================================
# Read data
# ==============================================================================

# ChromVAR
load("~/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/Results_with_GP_writeup6n/day10_chromVar_day10_growth_potential_for_week5_correlation_writeup6n.RData")

# cis_cor_vec <- correlation_list[['cis_cor_vec']]
# cocl2_cor_vec <- correlation_list[['cocl2_cor_vec']]
# dabtram_cor_vec <- correlation_list[['dabtram_cor_vec']]

# cis_cor_vec$motif_name <- rownames(cis_cor_vec)
# cocl2_cor_vec$motif_name <- rownames(cocl2_cor_vec)
# dabtram_cor_vec$motif_name <- rownames(dabtram_cor_vec)

# colnames(cis_cor_vec) <- c( "corr.cis.chromVAR", "p.val.cis.chromVAR", "motif_name")
# colnames(cocl2_cor_vec) <- c("corr.cocl2.chromVAR", "p.val.cocl2.chromVAR", "motif_name")
# colnames(dabtram_cor_vec) <- c("corr.dabtram.chromVAR","p.val.dabtram.chromVAR", "motif_name")

colnames(cis_cor_vec) <- c( "motif_name", "corr.cis.chromVAR", "p.val.cis.chromVAR")
colnames(cocl2_cor_vec) <- c("motif_name", "corr.cocl2.chromVAR", "p.val.cocl2.chromVAR")
colnames(dabtram_cor_vec) <- c("motif_name", "corr.dabtram.chromVAR","p.val.dabtram.chromVAR")


cor_df <- merge(cis_cor_vec,cocl2_cor_vec, by='motif_name' )
cor_df <- merge(cor_df,dabtram_cor_vec, by='motif_name' )
cor_df <- cor_df[, c("motif_name", "corr.cis.chromVAR", "corr.cocl2.chromVAR", "corr.dabtram.chromVAR",
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
# Identifying targets [ChromVAR]
# ==============================================================================

# CIS and COCL2
cor_mat_chromVAR_cis_cocl2 <- cor_df_chromVAR[, c("motif_name", 
                                                  "corr.cis.chromVAR", 
                                                  "corr.cocl2.chromVAR")]
rownames(cor_mat_chromVAR_cis_cocl2) <- cor_mat_chromVAR_cis_cocl2$motif_name
cor_mat_chromVAR_cis_cocl2 <- subset(cor_mat_chromVAR_cis_cocl2, select=-c(motif_name))
db_chromVAR_cis_cocl2 <- dbscan(cor_mat_chromVAR_cis_cocl2, eps = 0.02, minPts = 2) #eps = 0.03, minPts = 2
db_chromVAR_cis_cocl2
pairs(cor_mat_chromVAR_cis_cocl2, col = db_chromVAR_cis_cocl2$cluster + 1L)

db_chromVAR_cis_cocl2_cluster <- as.data.frame(db_chromVAR_cis_cocl2$cluster)
colnames(db_chromVAR_cis_cocl2_cluster) <- 'cluster_cis_cocl2'
db_chromVAR_cis_cocl2_cluster$motif_name <- rownames(cor_mat_chromVAR_cis_cocl2)

# DABTRAM and COCL2
cor_mat_chromVAR_dabtram_cocl2 <- cor_df_chromVAR[, c("motif_name", 
                                                      "corr.dabtram.chromVAR", 
                                                      "corr.cocl2.chromVAR")]
rownames(cor_mat_chromVAR_dabtram_cocl2) <- cor_mat_chromVAR_dabtram_cocl2$motif_name
cor_mat_chromVAR_dabtram_cocl2 <- subset(cor_mat_chromVAR_dabtram_cocl2, select=-c(motif_name))
db_chromVAR_dabtram_cocl2 <- dbscan(cor_mat_chromVAR_dabtram_cocl2, eps = 0.06, minPts = 1) #eps = 0.06, minPts = 1
db_chromVAR_dabtram_cocl2
pairs(cor_mat_chromVAR_dabtram_cocl2, col = db_chromVAR_dabtram_cocl2$cluster + 1L)

db_chromVAR_dabtram_cocl2_cluster <- as.data.frame(db_chromVAR_dabtram_cocl2$cluster)
colnames(db_chromVAR_dabtram_cocl2_cluster) <- 'cluster_dabtram_cocl2'
db_chromVAR_dabtram_cocl2_cluster$motif_name <- rownames(cor_mat_chromVAR_dabtram_cocl2)

# plotting
chromVAR_cluster <- cor_df_chromVAR[, c("motif_name", "corr.dabtram.chromVAR", "corr.cocl2.chromVAR", "corr.cis.chromVAR")]
chromVAR_cluster <- merge(chromVAR_cluster, db_chromVAR_cis_cocl2_cluster, by='motif_name')
chromVAR_cluster <- merge(chromVAR_cluster, db_chromVAR_dabtram_cocl2_cluster, by='motif_name')

chromVAR_cluster$cluster_dabtram_cocl2 <- as.factor(chromVAR_cluster$cluster_dabtram_cocl2)
chromVAR_cluster$cluster_cis_cocl2 <- as.factor(chromVAR_cluster$cluster_cis_cocl2)

p1 <- ggplot(chromVAR_cluster) +
  geom_point(aes(x = corr.cocl2.chromVAR, y = corr.dabtram.chromVAR, col=cluster_dabtram_cocl2)) +
  theme_bw()
p2 <- ggplot(chromVAR_cluster) +
  geom_point(aes(x = corr.cocl2.chromVAR, y = corr.cis.chromVAR, col=cluster_dabtram_cocl2)) +
  theme_bw()
grid.arrange(p1, p2, nrow=1)

ggplot(chromVAR_cluster) +
  geom_point(aes(x = corr.cocl2.chromVAR, y = corr.dabtram.chromVAR), color='gray') +
  geom_point(data = chromVAR_cluster[chromVAR_cluster$cluster_dabtram_cocl2 == 6, ], aes(x = corr.cocl2.chromVAR, y = corr.dabtram.chromVAR, col = cluster_dabtram_cocl2)) +
  theme_bw()
chromVAR_cluster2 <- chromVAR_cluster[chromVAR_cluster$cluster_dabtram_cocl2 == 2, ]
chromVAR_cluster3 <- chromVAR_cluster[chromVAR_cluster$cluster_dabtram_cocl2 == 3, ]
chromVAR_cluster4 <- chromVAR_cluster[chromVAR_cluster$cluster_dabtram_cocl2 == 4, ]
chromVAR_cluster5 <- chromVAR_cluster[chromVAR_cluster$cluster_dabtram_cocl2 == 5, ]
chromVAR_cluster6 <- chromVAR_cluster[chromVAR_cluster$cluster_dabtram_cocl2 == 6, ]
chromVAR_cluster7 <- chromVAR_cluster[chromVAR_cluster$cluster_dabtram_cocl2 == 7, ]
chromVAR_cluster8 <- chromVAR_cluster[chromVAR_cluster$cluster_dabtram_cocl2 == 8, ]
chromVAR_cluster9 <- chromVAR_cluster[chromVAR_cluster$cluster_dabtram_cocl2 == 9, ]
chromVAR_cluster10 <- chromVAR_cluster[chromVAR_cluster$cluster_dabtram_cocl2 == 10, ]
chromVAR_cluster2_3 <- rbind(chromVAR_cluster2, chromVAR_cluster3)
# write.csv(chromVAR_cluster2_3, '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/common_motifs_in_corr_with_growth_v2.csv', row.names = FALSE)

# ==============================================================================
# Identifying targets [RNA]
# ==============================================================================

# CIS and COCL2
cor_mat_RNA_cis_cocl2 <- cor_df_RNA[, c("gene", "corr.cis.RNA", "corr.cocl2.RNA")]
rownames(cor_mat_RNA_cis_cocl2) <- cor_mat_RNA_cis_cocl2$gene
cor_mat_RNA_cis_cocl2 <- subset(cor_mat_RNA_cis_cocl2, select=-c(gene))
eps_cis_cocl2 <- 0.016
minPts_cis_cocl2 <- 8
db_RNA_cis_cocl2 <- dbscan(cor_mat_RNA_cis_cocl2, eps = eps_cis_cocl2, minPts = minPts_cis_cocl2)
db_RNA_cis_cocl2
pairs(cor_mat_RNA_cis_cocl2, col = db_RNA_cis_cocl2$cluster + 1L)


db_RNA_cis_cocl2_cluster <- as.data.frame(db_RNA_cis_cocl2$cluster)
colnames(db_RNA_cis_cocl2_cluster) <- 'cluster_cis_cocl2'
db_RNA_cis_cocl2_cluster$gene <- rownames(cor_mat_RNA_cis_cocl2)

# DABTRAM and COCL2
cor_mat_RNA_dabtram_cocl2 <- cor_df_RNA[, c("gene", "corr.dabtram.RNA", "corr.cocl2.RNA")]
rownames(cor_mat_RNA_dabtram_cocl2) <- cor_mat_RNA_dabtram_cocl2$gene
cor_mat_RNA_dabtram_cocl2 <- subset(cor_mat_RNA_dabtram_cocl2, select=-c(gene))
eps_dabtram_cocl2 <- 0.019
minPts_dabtram_cocl2 <- 8
db_RNA_dabtram_cocl2 <- dbscan(cor_mat_RNA_dabtram_cocl2, eps = eps_dabtram_cocl2, minPts = minPts_dabtram_cocl2)
db_RNA_dabtram_cocl2
pairs(cor_mat_RNA_dabtram_cocl2, col = db_RNA_dabtram_cocl2$cluster + 1L)

db_RNA_dabtram_cocl2_cluster <- as.data.frame(db_RNA_dabtram_cocl2$cluster)
colnames(db_RNA_dabtram_cocl2_cluster) <- 'cluster_dabtram_cocl2'
db_RNA_dabtram_cocl2_cluster$gene <- rownames(cor_mat_RNA_dabtram_cocl2)

# plotting
RNA_cluster <- cor_df_RNA[, c("gene", "corr.dabtram.RNA", "corr.cocl2.RNA", "corr.cis.RNA")]
RNA_cluster <- merge(RNA_cluster, db_RNA_cis_cocl2_cluster, by='gene')
RNA_cluster <- merge(RNA_cluster, db_RNA_dabtram_cocl2_cluster, by='gene')

RNA_cluster$cluster_dabtram_cocl2 <- as.factor(RNA_cluster$cluster_dabtram_cocl2)
RNA_cluster$cluster_cis_cocl2 <- as.factor(RNA_cluster$cluster_cis_cocl2)

RNA_cluster_no0s <- RNA_cluster[RNA_cluster$cluster_dabtram_cocl2 != 0 & RNA_cluster$cluster_cis_cocl2 != 0, ]

p1 <- ggplot(RNA_cluster) +
  geom_point(aes(x = corr.cocl2.RNA, y = corr.dabtram.RNA, col=cluster_dabtram_cocl2), alpha=0.3, color='gray') +
  geom_point(data = RNA_cluster_no0s, aes(x = corr.cocl2.RNA, y = corr.dabtram.RNA, col=cluster_cis_cocl2)) +
  theme_bw() + theme(legend.position = "none")
p2 <- ggplot(RNA_cluster) +
  geom_point(aes(x = corr.cocl2.RNA, y = corr.cis.RNA, col=cluster_dabtram_cocl2), alpha=0.3, color='gray') +
  geom_point(data = RNA_cluster_no0s, aes(x = corr.cocl2.RNA, y = corr.cis.RNA, col=cluster_cis_cocl2)) +
  theme_bw()  #+ theme(legend.position = "none")
grid.arrange(p1, p2, nrow=1)

RNA_cluster4 <- RNA_cluster[RNA_cluster$cluster_dabtram_cocl2 == 4 & RNA_cluster$cluster_cis_cocl2 == 4, ]
# write.csv(RNA_cluster4, '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/common_genes_in_corr_with_growth_v2.csv', row.names = FALSE)
RNA_cluster3 <- RNA_cluster[RNA_cluster$cluster_cis_cocl2 == 3, ]
RNA_cluster5 <- RNA_cluster[RNA_cluster$cluster_cis_cocl2 == 5, ]

write.csv(RNA_cluster3, '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/common_genes_in_corr_with_growth_v3_cis_cocl2_clust3.csv', row.names = FALSE)
write.csv(RNA_cluster5, '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/common_genes_in_corr_with_growth_v3_cis_cocl2_clust5.csv', row.names = FALSE)


ggplot(RNA_cluster) +
  geom_point(aes(x = corr.cocl2.RNA, y = corr.cis.RNA, col=cluster_cis_cocl2), alpha=0.3, color='gray') +
  geom_point(data = RNA_cluster[RNA_cluster$cluster_cis_cocl2 == 5, ], aes(x = corr.cocl2.RNA, y = corr.cis.RNA, col=cluster_cis_cocl2)) +
  theme_bw()

