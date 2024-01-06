library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(rrvgo)

# ==============================================================================
# Read data
# ==============================================================================
# input_dir <- '~/Dropbox/Thesis/Lineage_trace/outputs/task0_explore_lineage_variability/features_associated_with_plasticity/'
# genes_to_test <- read.csv(paste0(input_dir, 'day10_DABTRAM_gene_exp_mean_day10_plasticity_correlation_TOP1000.csv'))
# genes_to_test <- genes_to_test$features_to_label

input_dir <- '~/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/'
genes_to_test <- read.csv(paste0(input_dir, 'common_genes_in_corr_with_growth_v2.csv'))
# genes_to_test <- read.csv(paste0(input_dir, 'common_genes_in_corr_with_growth_v3_cis_cocl2_clust3.csv'))
genes_to_test <- genes_to_test$gene
# ==============================================================================
# GO term
# ==============================================================================

# Performing GO
ego <- clusterProfiler::enrichGO(gene          = genes_to_test, #TODO: determine a list of background genes
                                 OrgDb         = org.Hs.eg.db::org.Hs.eg.db,
                                 keyType       = "SYMBOL",
                                 ont           = "BP",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 0.05,
                                 qvalueCutoff  = 0.05,
                                 readable      = TRUE)
head(ego)

if(length(ego$ID) > 5) {
  
  ## see https://www.bioconductor.org/packages/release/bioc/vignettes/rrvgo/inst/doc/rrvgo.html
  simMatrix <- rrvgo::calculateSimMatrix(ego$ID,
                                         orgdb="org.Hs.eg.db",
                                         ont="BP",
                                         method="Rel")
  
  scores <- setNames(-log10(ego$qvalue), ego$ID)
  reducedTerms <- rrvgo::reduceSimMatrix(simMatrix,
                                         scores,
                                         threshold=0.7,
                                         orgdb="org.Hs.eg.db")
  
  # png(paste0("../../../out/fig/main/sns_", names(file_vec)[kk], "_revigo-treemap.png"),
  #     height = 2000, width = 2000,
  #     units = "px", res = 500)
  rrvgo::treemapPlot(reducedTerms)
  # graphics.off()
}
