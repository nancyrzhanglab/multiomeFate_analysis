library(Seurat)
library(tidyverse)
library(data.table)
library(ggplot2)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
result_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/'

remove_unassigned_cells <- TRUE

# ==============================================================================
# Read data general
# ==============================================================================

load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_saver.RData'))

all_data[["Saver"]] <- all_data_saver
# remove cells with no lineage
if(remove_unassigned_cells) {
  print("Removing cells with no assigned lineage")
  all_data$keep <- !is.na(all_data$assigned_lineage)
  if(any(!all_data$keep)){
    print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
    all_data <- subset(all_data, keep == TRUE)
  }
}


# ==============================================================================
# Calculate mean expression
# ==============================================================================

all_data.day0 <- subset(all_data, dataset == 'day0')
all_data.day0.saver <- all_data.day0$Saver
all_data.day0.saver.mean <- rowMeans(all_data.day0.saver)
all_data.day0.saver.mean <- as.data.frame(all_data.day0.saver.mean)
all_data.day0.saver.mean$gene <- rownames(all_data.day0.saver.mean)

# ==============================================================================
# Read correlation results
# ==============================================================================
load(paste0(result_dir, 'saver_cor_vec.RData'))

# ==============================================================================
# Get and sort table
# ==============================================================================

getAndSortTable <- function(name) {
  cor_vec <- as.data.frame(saver_cor_vec[[name]]) %>% drop_na()
  cor_vec$p.value <- as.numeric(cor_vec$p.value)
  cor_vec$p_adj <- p.adjust(cor_vec$p.value, method = 'BH')
  cor_vec <- cor_vec[cor_vec$p_adj < 0.05, ]
  
  cor_vec$gene <- rownames(cor_vec)
  cor_vec$correlation <- as.numeric(cor_vec$correlation)
  cor_vec <- cor_vec[order(cor_vec$correlation, decreasing = TRUE),]
  
  return(cor_vec)
}

day0_dabtram <- getAndSortTable('dabtram_d0_saver_cor_vec')
day0_cocl2 <- getAndSortTable('cocl2_d0_saver_cor_vec')
day0_cis <- getAndSortTable('cis_d0_saver_cor_vec')

day10_dabtram <- getAndSortTable('dabtram_d10_saver_cor_vec')
day10_cocl2 <- getAndSortTable('cocl2_d10_saver_cor_vec')
day10_cis <- getAndSortTable('cis_d10_saver_cor_vec')

# ==============================================================================
# Plot correlations
# ==============================================================================
day0_dabtram.plot <- day0_dabtram[, c('gene', 'correlation')]
colnames(day0_dabtram.plot) <- c('gene', 'dabtram_correlation')
day0_dabtram.plot <- day0_dabtram.plot[day0_dabtram.plot$dabtram_correlation > 0, ]

day0_cocl2.plot <- day0_cocl2[, c('gene', 'correlation')]
colnames(day0_cocl2.plot) <- c('gene', 'cocl2_correlation')
day0_cocl2.plot <- day0_cocl2.plot[day0_cocl2.plot$cocl2_correlation > 0, ]

day0_cis.plot <- day0_cis[, c('gene', 'correlation')]
colnames(day0_cis.plot) <- c('gene', 'cis_correlation')
day0_cis.plot <- day0_cis.plot[day0_cis.plot$cis_correlation > 0, ]

day0_correlations <- merge(day0_dabtram.plot, day0_cocl2.plot, by = 'gene')
day0_correlations <- merge(day0_correlations, day0_cis.plot, by = 'gene')

ggpairs(day0_correlations, columns = 2:4, 
        lower = list(continuous = wrap("points", alpha = 0.5, size = 0.5)))

day0_correlations <- day0_correlations[!grepl('\\.', day0_correlations$gene), ]

eg = bitr(day0_correlations$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
x <- enrichPathway(gene=eg$ENTREZID, pvalueCutoff = 0.05, readable=TRUE)
x_res <- x@result
x_res <- x_res[x_res$qvalue < 0.1, ]

# Performing GO
ego <- clusterProfiler::enrichGO(gene          = eg$SYMBOL, #TODO: determine a list of background genes
                                 OrgDb         = org.Hs.eg.db::org.Hs.eg.db,
                                 keyType       = "SYMBOL",
                                 ont           = "MF",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 0.05,
                                 qvalueCutoff  = 0.05,
                                 readable      = TRUE)
head(ego)
ego_res <- ego@result
ego_res <- ego_res[ego_res$qvalue < 0.1, ]

# ==============================================================================
# Take day0 intersect
# ==============================================================================

day0_dabtram <- day0_dabtram[!grepl('\\.', day0_dabtram$gene), ]
day0_dabtram <- merge(day0_dabtram, all_data.day0.saver.mean, by = 'gene', all.x = TRUE)
day0_dabtram <- day0_dabtram[day0_dabtram$all_data.day0.saver.mean > 0.1, ]
day0_dabtram <- day0_dabtram[order(day0_dabtram$correlation, decreasing = TRUE),]

day0_cocl2 <- day0_cocl2[!grepl('\\.', day0_cocl2$gene), ]
day0_cocl2 <- merge(day0_cocl2, all_data.day0.saver.mean, by = 'gene', all.x = TRUE)
day0_cocl2 <- day0_cocl2[day0_cocl2$all_data.day0.saver.mean > 0.1, ]
day0_cocl2 <- day0_cocl2[order(day0_cocl2$correlation, decreasing = TRUE),]

day0_cis <- day0_cis[!grepl('\\.', day0_cis$gene), ]
day0_cis <- merge(day0_cis, all_data.day0.saver.mean, by = 'gene', all.x = TRUE)
day0_cis <- day0_cis[day0_cis$all_data.day0.saver.mean > 0.1, ]
day0_cis <- day0_cis[order(day0_cis$correlation, decreasing = TRUE),]

n_pct <- 0.2
day0_dabtram.topN <- day0_dabtram[1:as.integer(nrow(day0_dabtram) * n_pct), ]$gene
day0_cocl2.topN <- day0_cocl2[1:as.integer(nrow(day0_cocl2) * n_pct), ]$gene
day0_cis.topN <- day0_cis[1:as.integer(nrow(day0_cis) * n_pct), ]$gene

day0_intersect <- intersect(intersect(day0_dabtram.topN, day0_cocl2.topN), day0_cis.topN)
intersect <- as.data.frame(day0_intersect)
colnames(intersect) <- c('gene')

intersect <- merge(intersect, all_data.day0.saver.mean, by = 'gene', all.x = TRUE)
intersect <- intersect[!grepl('\\.', intersect$gene), ]

write.csv(intersect, paste0('~/Downloads/day0_intersect.csv'), row.names = FALSE)
