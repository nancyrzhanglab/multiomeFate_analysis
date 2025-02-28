rm(list = ls())
library(tidyverse)
library(GSEABase)
library(Biobase) #base functions for bioconductor; required by GSEABase
library(clusterProfiler) # provides a suite of tools for functional enrichment analysis
library(enrichplot) # great for making the standard GSEA enrichment plots
library(ggplot2)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
result_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/'
result_dir2 <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_differential_winner_loser_d0_V2/'

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
load(paste0(result_dir2, 'differential_winner_loser_saver_lineage_specific_adaptation_gene_DABTRAM_t_test_results.RData'))
diff_res_dabtram <- t_test_results


hallmark <- read.gmt("~/Downloads/h.all.v2024.1.Hs.symbols.gmt")
reactome <- read.gmt("~/Downloads/c2.cp.reactome.v2024.1.Hs.symbols.gmt")
threeCA <- read.gmt("~/Downloads/c4.3ca.v2024.1.Hs.symbols.gmt")

# ==============================================================================
# GSEA
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
# Wrangle diff test results
# ==============================================================================
wrangle_data <- function(diff_res) {
  diff_res$p_val <- as.numeric(diff_res$p_val)
  diff_res$p_adj <- p.adjust(diff_res$p_val, method = 'BH')
  diff_res$neg_log10_p_val <- -log10(diff_res$p_val)
  
  diff_res$mean_winner <- as.numeric(diff_res$mean_winner)
  diff_res$mean_other <- as.numeric(diff_res$mean_other)
  diff_res$fold_change <- log2(diff_res$mean_winner / diff_res$mean_other)
  
  # order by fold_change in decreasing order
  diff_res <- diff_res[order(diff_res$fold_change, decreasing = TRUE),]
  
  return(diff_res)
}
diff_res_dabtram <- wrangle_data(diff_res_dabtram)
diff_res_dabtram <- diff_res_dabtram[diff_res_dabtram$p_adj < 0.05, ]
diff_res_dabtram <- diff_res_dabtram[order(diff_res_dabtram$fold_change, decreasing = TRUE), ]

day0_dabtram <- day0_dabtram[day0_dabtram$gene %in% diff_res_dabtram$feature,]


# gsea_input <- diff_res_dabtram$fold_change
# names(gsea_input) <- diff_res_dabtram$feature

day0_dabtram <- merge(day0_dabtram, all_data.day0.saver.mean, by = 'gene')
day0_dabtram <- day0_dabtram[!grepl('\\.', day0_dabtram$gene),]
day0_dabtram <- day0_dabtram[order(day0_dabtram$correlation, decreasing = TRUE),]
day0_dabtram <- day0_dabtram[day0_dabtram$all_data.day0.saver.mean > 0.1, ]
day0_dabtram.topN <- day0_dabtram[1:as.integer(nrow(day0_dabtram) * 0.25), ]
day0_dabtram.topN <- day0_dabtram.topN[day0_dabtram.topN$correlation > 0 ,]

day0_dabtram.bottomN <- tail(day0_dabtram, as.integer(nrow(day0_dabtram) * 0.25))
day0_dabtram.bottomN <- day0_dabtram.bottomN[day0_dabtram.bottomN$correlation < 0 ,]

day0_cocl2 <- merge(day0_cocl2, all_data.day0.saver.mean, by = 'gene')
day0_cocl2 <- day0_cocl2[!grepl('\\.', day0_cocl2$gene),]
day0_cocl2 <- day0_cocl2[order(day0_cocl2$correlation, decreasing = TRUE),]
day0_cocl2 <- day0_cocl2[day0_cocl2$all_data.day0.saver.mean > 0.1, ]
day0_cocl2.topN <- day0_cocl2[1:as.integer(nrow(day0_cocl2) * 0.25), ]
day0_cocl2.topN <- day0_cocl2.topN[day0_cocl2.topN$correlation > 0 ,]

day0_cocl2.bottomN <- tail(day0_cocl2, as.integer(nrow(day0_cocl2) * 0.25))
day0_cocl2.bottomN <- day0_cocl2.bottomN[day0_cocl2.bottomN$correlation < 0 ,]

day0_cis <- merge(day0_cis, all_data.day0.saver.mean, by = 'gene')
day0_cis <- day0_cis[!grepl('\\.', day0_cis$gene),]
day0_cis <- day0_cis[order(day0_cis$correlation, decreasing = TRUE),]
day0_cis <- day0_cis[day0_cis$all_data.day0.saver.mean > 0.1, ]
day0_cis.topN <- day0_cis[1:as.integer(nrow(day0_cis) * 0.25), ]
day0_cis.topN <- day0_cis.topN[day0_cis.topN$correlation > 0 ,]

day0_cis.bottomN <- tail(day0_cis, as.integer(nrow(day0_cis) * 0.25))
day0_cis.bottomN <- day0_cis.bottomN[day0_cis.bottomN$correlation < 0 ,]

# write.csv(day0_dabtram, paste0('~/Downloads/day0_dabtram.csv'))

# GSEA DABTRAM
gsea_input_dabtram <- day0_dabtram$correlation
names(gsea_input_dabtram) <- day0_dabtram$gene

set.seed(123)
GSEA_res_dabtram <- GSEA(geneList = gsea_input_dabtram, 
                 TERM2GENE = threeCA, 
                 pvalueCutoff = 0.2,
                 seed = T,
                 verbose = F)
GSEA_res_dabtram.df <- as_tibble(GSEA_res_dabtram@result)
GSEA_res_dabtram.df <- GSEA_res_dabtram.df[, c('ID', 'NES', 'p.adjust', 'qvalue')]
colnames(GSEA_res_dabtram.df) <- c('ID', 'NES.DABTRAM', 'p.adjust.DABTRAM', 'qvalue.DABTRAM')

# GSEA COCL2
gsea_input_cocl2 <- day0_cocl2$correlation
names(gsea_input_cocl2) <- day0_cocl2$gene

set.seed(123)
GSEA_res_cocl2 <- GSEA(geneList = gsea_input_cocl2, 
                         TERM2GENE = reactome, 
                         seed = T,
                         verbose = F)
GSEA_res_cocl2.df <- as_tibble(GSEA_res_cocl2@result)
GSEA_res_cocl2.df <- GSEA_res_cocl2.df[, c('ID', 'NES', 'p.adjust', 'qvalue')]
colnames(GSEA_res_cocl2.df) <- c('ID', 'NES.COCL2', 'p.adjust.COCL2', 'qvalue.COCL2')

# GSEA CIS
gsea_input_cis <- day0_cis$correlation
names(gsea_input_cis) <- day0_cis$gene

set.seed(123)
GSEA_res_cis <- GSEA(geneList = gsea_input_cis, 
                       TERM2GENE = reactome, 
                       seed = T,
                       verbose = F)
GSEA_res_cis.df <- as_tibble(GSEA_res_cis@result)
GSEA_res_cis.df <- GSEA_res_cis.df[, c('ID', 'NES', 'p.adjust', 'qvalue')]
colnames(GSEA_res_cis.df) <- c('ID', 'NES.CIS', 'p.adjust.CIS', 'qvalue.CIS')

GSEA_res.df <- merge(GSEA_res_dabtram.df, GSEA_res_cocl2.df, by = c('ID'), all = T)
GSEA_res.df <- merge(GSEA_res.df, GSEA_res_cis.df, by = c('ID'), all = T)

# create enrichment plots using the enrichplot package
gseaplot2(GSEA_res, 
          geneSetID = c(1, 2, 3), #can choose multiple signatures to overlay in this plot
          pvalue_table = FALSE, #can set this to FALSE for a cleaner plot
          title = "Top 3 Significant Pathways") #can also turn off this title

GSEA_res.df <- GSEA_res.df %>%
  mutate(phenotype = case_when(
    NES > 0 ~ "Positive",
    NES < 0 ~ "Negative"))

# create 'bubble plot' to summarize y signatures across x phenotypes
ggplot(GSEA_res.df[1:14,], aes(x=phenotype, y=ID)) + 
  geom_point(aes(size=setSize, color = NES, alpha=-log10(p.adjust))) +
  scale_color_gradient(low="blue", high="red") +
  theme_bw()

library(ReactomePA)
eg = bitr(day0_dabtram.topN$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
x <- enrichPathway(gene=eg$ENTREZID, pvalueCutoff = 0.05, readable=TRUE)
x_dabtram <- x@result
x_dabtram <- x_dabtram[x_dabtram$p.adjust < 0.05, ]
x_dabtram <- x_dabtram[, c('ID', 'Description', 'GeneRatio', 'p.adjust')]
colnames(x_dabtram) <- c('ID', 'Description', 'GeneRatio.DABTRAM', 'p.adjust.DABTRAM')

eg = bitr(day0_cocl2.topN$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
x <- enrichPathway(gene=eg$ENTREZID, pvalueCutoff = 0.05, readable=TRUE)
x_cocl2 <- x@result
x_cocl2 <- x_cocl2[x_cocl2$p.adjust < 0.05, ]
x_cocl2 <- x_cocl2[, c('ID', 'Description', 'GeneRatio', 'p.adjust')]
colnames(x_cocl2) <- c('ID', 'Description', 'GeneRatio.COCL2', 'p.adjust.COCL2')

eg = bitr(day0_cis.topN$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
x <- enrichPathway(gene=eg$ENTREZID, pvalueCutoff = 0.05, readable=TRUE)
x_cis <- x@result
x_cis <- x_cis[x_cis$p.adjust < 0.05, ]
x_cis <- x_cis[, c('ID', 'Description', 'GeneRatio', 'p.adjust')]
colnames(x_cis) <- c('ID', 'Description', 'GeneRatio.CIS', 'p.adjust.CIS')

x_all <- merge(x_dabtram, x_cocl2, by = c('ID', 'Description'), all = T)
x_all <- merge(x_all, x_cis, by = c('ID', 'Description'), all = T)


gene_intersect1 <- intersect(day0_dabtram.topN$gene, day0_cocl2.topN$gene)
gene_intersect2 <- intersect(day0_dabtram.topN$gene, day0_cis.topN$gene)
gene_intersect3 <- intersect(day0_cocl2.topN$gene, day0_cis.topN$gene)
gene_intersect <- union(gene_intersect1, gene_intersect2)
gene_intersect <- union(gene_intersect, gene_intersect3)

gene_intersect <- intersect(day0_dabtram.topN$gene, day0_cocl2.topN$gene)
gene_intersect <- intersect(gene_intersect, day0_cis.topN$gene)

gene_intersect.bottom <- intersect(day0_dabtram.bottomN$gene, day0_cocl2.bottomN$gene)
gene_intersect.bottom <- intersect(gene_intersect.bottom, day0_cis.bottomN$gene)

eg = bitr(gene_intersect.bottom, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
x <- enrichPathway(gene=eg$ENTREZID, pvalueCutoff = 0.05, readable=TRUE)
x_res <- x@result
x_res <- x_res[x_res$p.adjust < 0.2, ]
x_res <- x_res[, c('ID', 'Description', 'GeneRatio', 'p.adjust')]

# Performing GO
ego <- clusterProfiler::enrichGO(gene          = eg$SYMBOL, #TODO: determine a list of background genes
                                 OrgDb         = org.Hs.eg.db::org.Hs.eg.db,
                                 keyType       = "SYMBOL",
                                 ont           = "BP",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 0.05,
                                 qvalueCutoff  = 0.05,
                                 readable      = TRUE)
head(ego)
ego_res <- ego@result
ego_res <- ego_res[ego_res$p.adjust < 0.1, ]
