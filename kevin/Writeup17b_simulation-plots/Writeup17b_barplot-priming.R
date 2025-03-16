rm(list=ls())
library(Seurat)
library(multiomeFate)
library(ggplot2)
library(Ckmeans.1d.dp)

# version on HPC3
# out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup14/Writeup14_priming-setting_"
# plot_folder <- "~/project/Multiome_fate/git/multiomeFate_analysis_kevin/fig/kevin/Writeup14/Writeup14_priming_v2-setting_"
# func_folder <- "~/project/Multiome_fate/git/multiomeFate_analysis_kevin/kevin/Writeup14_simulation/"
# load(paste0(out_folder, "simulation_v2.RData"))
# load(paste0(out_folder, "v2_seurat_CoSpar_prepared.RData"))

# version locally
out_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/out/Writeup14/Writeup14_priming-setting_"
plot_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/git/multiomeFate_analysis/fig/kevin/Writeup17b/Writeup17b_priming_"
csv_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/out/Writeup14/"

func_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/git/multiomeFate_analysis/kevin/Writeup14_simulation/"
load(paste0(out_folder, "simulation_v2.RData"))
load(paste0(out_folder, "v2_seurat_CoSpar_prepared.RData"))

source("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/git/multiomeFate_analysis/kevin/Writeup17b_simulation-plots/plot_signed_logpvalue.R")
source("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/git/multiomeFate_analysis/kevin/Writeup17b_simulation-plots/multtest_correction.R")

####


source(paste0(func_folder, "func_seurat.R"))
names(simulation_res$lineage_assignment) <- rownames(simulation_res$embedding_mat)
simulation_data <- .form_simulation_seurat_fate(final_fit = final_fit,
                                                simulation_res = simulation_res)

################
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
tab_mat <- tab_mat[order(tab_mat[,"week5_COCL2"], decreasing = TRUE),]

# compute the "denoised" gene expression profiles
all_data <- subset(all_data, dataset == "day10_COCL2")
cell_embedding <- all_data[["fasttopic.COCL2"]]@cell.embeddings
gene_embedding <- all_data[["fasttopic.COCL2"]]@feature.loadings
denoised_mat <- tcrossprod(cell_embedding, gene_embedding)
denoised_mat <- denoised_mat[Seurat::Cells(simulation_data),]
p <- ncol(denoised_mat)

#############

# other way to define the truth:
# Separate cells by true fate potential positive or negative
winner_cells <- which(simulation_data$fatepotential_true >= 0)
loser_cells <- which(simulation_data$fatepotential_true < 0)
true_cell_labels <- rep("Loser", length(simulation_data$assigned_lineage))
true_cell_labels[winner_cells] <- "Winner"
wilcox_results <- sapply(1:p, function(j){
  tmp <- stats::wilcox.test(
    x = denoised_mat[winner_cells,j],
    y = denoised_mat[loser_cells,j]
  )
  logfc <- log2(mean(denoised_mat[winner_cells,j])) - log2(mean(denoised_mat[loser_cells,j]))
  
  c(logFC = logfc,
    pval = tmp$p.value)
})
wilcox_results <- t(wilcox_results)
gene_df <- as.data.frame(wilcox_results)
rownames(gene_df) <- colnames(denoised_mat)
colnames(gene_df) <- c("true_logFC", "true_pval")
gene_df$Gene <- rownames(gene_df)

##########

# add naive
tmp <- log10(simulation_data@misc[["lineage_observed_count"]]+1)
winner_lineages <- names(sort(tmp, decreasing = TRUE)[1:10])
tab_mat[winner_lineages,] # verify the sizes
loser_lineages <- setdiff(names(tmp), winner_lineages)

winner_cells <- which(simulation_data$assigned_lineage %in% winner_lineages)
loser_cells <- which(simulation_data$assigned_lineage %in% loser_lineages)

wilcox_results <- sapply(1:p, function(j){
  tmp <- stats::wilcox.test(
    x = denoised_mat[winner_cells,j],
    y = denoised_mat[loser_cells,j]
  )
  logfc <- log2(mean(denoised_mat[winner_cells,j])) - log2(mean(denoised_mat[loser_cells,j]))
  
  c(logFC = logfc,
    pval = tmp$p.value)
})
wilcox_results <- t(wilcox_results)
rownames(wilcox_results) <- colnames(denoised_mat)
gene_df$naive_logFC <- wilcox_results[rownames(gene_df),"logFC"]
gene_df$naive_pval <- wilcox_results[rownames(gene_df),"pval"]

#############

# add cospar
df <- read.csv(paste0(csv_folder, "simulation_priming_cospar_obs_v2.csv"))
rownames(df) <- df$X

cospar_fate <- df$fate_bias_intraclone_transition_map_High.Low
names(cospar_fate) <- rownames(df)
cospar_fate <- cospar_fate[rownames(denoised_mat)]

winner_cells <- names(cospar_fate)[cospar_fate > 0.6]
loser_cells <- names(cospar_fate)[cospar_fate < 0.4]

cospar_results <- sapply(1:p, function(j){
  tmp <- stats::wilcox.test(
    x = denoised_mat[winner_cells,j],
    y = denoised_mat[loser_cells,j]
  )
  logfc <- log2(mean(denoised_mat[winner_cells,j])) - log2(mean(denoised_mat[loser_cells,j]))
  
  c(logFC = logfc,
    pval = tmp$p.value)
})
cospar_results <- t(cospar_results)
gene_df$cospar_logFC <- cospar_results[,"logFC"]
gene_df$cospar_pval <- cospar_results[,"pval"]

#############

# add cypher
cypher_fate <- simulation_data$fatepotential[rownames(denoised_mat)]
cypher_fate <- cypher_fate[rownames(denoised_mat)]

winner_cells <- names(cypher_fate)[cypher_fate > 0]
loser_cells <- names(cypher_fate)[cypher_fate < -0.1]

cypher_results <- sapply(1:p, function(j){
  tmp <- stats::wilcox.test(
    x = denoised_mat[winner_cells,j],
    y = denoised_mat[loser_cells,j]
  )
  logfc <- log2(mean(denoised_mat[winner_cells,j])) - log2(mean(denoised_mat[loser_cells,j]))
  
  c(logFC = logfc,
    pval = tmp$p.value)
})
cypher_results <- t(cypher_results)
gene_df$cypher_logFC <- cypher_results[,"logFC"]
gene_df$cypher_pval <- cypher_results[,"pval"]

#########

gene_list <- vector("list", length = 4)
names(gene_list) <- c("true", "naive", "cospar", "cypher")

tmp <- gene_df$true_pval; names(tmp) <- rownames(gene_df)
gene_df$true_pval <- .multtest_correction(
  pvalue_vec = tmp
)

tmp <- gene_df$true_pval; tmp <- stats::p.adjust(tmp, method = "BH"); names(tmp) <- rownames(gene_df)
gene_list[["true"]] <- names(tmp)[which(tmp <= 0.35)]

tmp <- gene_df$naive_pval; names(tmp) <- rownames(gene_df)
gene_list[["naive"]] <- names(tmp)[order(tmp, decreasing = FALSE)[1:length(gene_list[["true"]])]]

tmp <- gene_df$cospar_pval; names(tmp) <- rownames(gene_df)
gene_list[["cospar"]] <- names(tmp)[order(tmp, decreasing = FALSE)[1:length(gene_list[["true"]])]]

tmp <- gene_df$cypher_pval; names(tmp) <- rownames(gene_df)
gene_list[["cypher"]] <- names(tmp)[order(tmp, decreasing = FALSE)[1:length(gene_list[["true"]])]]

jaccard_vec <- sapply(2:4, function(i){
  numerator <- length(intersect(gene_list[["true"]], gene_list[[i]]))
  denominator <- length(unique(c(gene_list[["true"]], gene_list[[i]])))
  numerator/denominator
})
names(jaccard_vec) <- names(gene_list)[2:4]
jaccard_vec

