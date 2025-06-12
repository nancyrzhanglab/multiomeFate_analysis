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
plot_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/git/multiomeFate_analysis/fig/kevin/Writeup14/tmp_Writeup14_priming-setting_"
func_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/git/multiomeFate_analysis/kevin/Writeup14_simulation/"
load(paste0(out_folder, "simulation_v2.RData"))
load(paste0(out_folder, "v2_seurat_CoSpar_prepared.RData"))

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
  
  c(logfc = logfc,
    p.value = tmp$p.value)
})
wilcox_results <- t(wilcox_results)
gene_df <- as.data.frame(wilcox_results)
rownames(gene_df) <- colnames(denoised_mat)
colnames(gene_df) <- c("true_logfc", "true_p.value")
gene_df$true_p.value.adj <- stats::p.adjust(gene_df$true_p.value, method = "BH")

tmp <- Ckmeans.1d.dp::Ckmeans.1d.dp(
  x = abs(gene_df$true_logfc),
  k = 2
)
true_gene_vec <- rep(TRUE, nrow(gene_df))
true_gene_vec[which(gene_df$true_p.value.adj >= 0.05)] <- FALSE
true_gene_vec[which(tmp$cluster == which.min(tmp$centers))] <- FALSE
gene_df$true_gene <- true_gene_vec

##########
# Now do this for the estimated fate potential
winner_cells <- which(simulation_data$fatepotential >= 0)
loser_cells <- which(simulation_data$fatepotential < 0)
wilcox_results <- sapply(1:p, function(j){
  tmp <- stats::wilcox.test(
    x = denoised_mat[winner_cells,j],
    y = denoised_mat[loser_cells,j]
  )
  logfc <- log2(mean(denoised_mat[winner_cells,j])) - log2(mean(denoised_mat[loser_cells,j]))
  
  c(logfc = logfc,
    p.value = tmp$p.value)
})
wilcox_results <- t(wilcox_results)
rownames(wilcox_results) <- colnames(denoised_mat)
gene_df$est_logfc <- wilcox_results[rownames(gene_df),"logfc"]
gene_df$est_p.value <- wilcox_results[rownames(gene_df),"p.value"]
gene_df$est_p.value.adj <- stats::p.adjust(gene_df$est_p.value, method = "BH")
tmp <- Ckmeans.1d.dp::Ckmeans.1d.dp(
  x = abs(gene_df$est_logfc),
  k = 2
)
est_gene_vec <- rep(TRUE, nrow(gene_df))
est_gene_vec[which(gene_df$est_p.value.adj >= 0.05)] <- FALSE
est_gene_vec[which(tmp$cluster == which.min(tmp$centers))] <- FALSE
gene_df$est_gene <- est_gene_vec
# table(gene_df$true_gene, gene_df$est_gene)
gene_df$est_log10p.value <- -log10(gene_df$est_p.value)

label_vec <- sapply(1:p, function(j){
  if(gene_df$est_gene[j] == TRUE & gene_df$true_gene[j] == TRUE) {
    return("Both")
  } else if(gene_df$est_gene[j] == FALSE & gene_df$true_gene[j] == TRUE) {
    return("True")
  } else if(gene_df$est_gene[j] == TRUE & gene_df$true_gene[j] == FALSE) {
    return("Estimate")
  } else{
    return("Neither")
  }
})
gene_df$est_label <- label_vec
gene_df <- gene_df[c(which(gene_df$est_label == "Neither"),
                     which(gene_df$est_label == "Both"),
                     which(gene_df$est_label == "Estimate"),
                     which(gene_df$est_label == "True")),
]

colors <- setNames(c("#e0e0e0", "#cc0967", "#109163", "#bc6a17"), 
                   c("Neither", "True", "Both", "Estimate"))

cor_val <- stats::cor(gene_df$est_logfc, gene_df$true_logfc)
ymult_test <- min(-log10(gene_df$est_p.value[which(gene_df$est_p.value.adj <= 0.05)]))

# make a correlation logfc scatterplot
plot_1 <- ggplot(gene_df, aes(x = true_logfc, y = est_logfc, color = est_label)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = colors) +
  geom_vline(xintercept = 0, color = "#2f2f2f") +
  geom_hline(yintercept = 0, color = "#2f2f2f") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#2f2f2f") +
  labs(title = paste0("Winner-vs-loser cells based on estimated fate:",
                      "\nCorrelation in true-and-estimated logfc: ", 
                      round(cor_val,2),
                      "\nBoth: ", length(which(gene_df$est_label == "Both")),
                      ", Only true: ", length(which(gene_df$est_label == "True")),
                      ", Only estimate: ", length(which(gene_df$est_label == "Estimate")),
                      "\nNum. genes: ", nrow(gene_df)),
       x = "True LogFC",
       y = "Estimated LogFC") +
  theme_minimal() 
ggsave(plot_1, 
       file = paste0(plot_folder, "true-vs-estimated-logfc.png"),
       width = 5, height = 5, units = "in")

##########

# do a DE analysis between the winner lineages vs loser lineages
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
  
  c(logfc = logfc,
    p.value = tmp$p.value)
})
wilcox_results <- t(wilcox_results)
rownames(wilcox_results) <- colnames(denoised_mat)
gene_df$naive_logfc <- wilcox_results[rownames(gene_df),"logfc"]
gene_df$naive_p.value <- wilcox_results[rownames(gene_df),"p.value"]
gene_df$naive_p.value.adj <- stats::p.adjust(gene_df$naive_p.value, method = "BH")

tmp <- Ckmeans.1d.dp::Ckmeans.1d.dp(
  x = abs(gene_df$naive_logfc),
  k = 2
)
est_gene_vec <- rep(TRUE, nrow(gene_df))
est_gene_vec[which(gene_df$naive_p.value.adj >= 0.05)] <- FALSE
est_gene_vec[which(tmp$cluster == which.min(tmp$centers))] <- FALSE
gene_df$naive_gene <- est_gene_vec
# table(gene_df$true_gene, gene_df$naive_gene)
gene_df$naive_log10p.value <- -log10(gene_df$naive_p.value)

label_vec <- sapply(1:p, function(j){
  if(gene_df$naive_gene[j] == TRUE & gene_df$true_gene[j] == TRUE) {
    return("Both")
  } else if(gene_df$naive_gene[j] == FALSE & gene_df$true_gene[j] == TRUE) {
    return("True")
  } else if(gene_df$naive_gene[j] == TRUE & gene_df$true_gene[j] == FALSE) {
    return("Naive")
  } else{
    return("Neither")
  }
})
gene_df$naive_label <- label_vec
gene_df <- gene_df[c(which(gene_df$naive_label == "Neither"),
                     which(gene_df$naive_label == "Both"),
                     which(gene_df$naive_label == "Naive"),
                     which(gene_df$naive_label == "True")),
]

cor_val <- stats::cor(gene_df$naive_logfc, gene_df$true_logfc)

colors <- setNames(c("#e0e0e0", "#cc0967", "#109163", "#bc6a17"), 
                   c("Neither", "True", "Both", "Naive"))

# make a correlation logfc scatterplot
plot_1 <- ggplot(gene_df, aes(x = true_logfc, y = naive_logfc, color = naive_label)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = colors) +
  geom_vline(xintercept = 0, color = "#2f2f2f") +
  geom_hline(yintercept = 0, color = "#2f2f2f") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#2f2f2f") +
  labs(title = paste0("Winner-vs-loser lineages (naive):",
                      "\nCorrelation in true-and-naive logfc: ", 
                      round(cor_val,2),
                      "\nBoth: ", length(which(gene_df$naive_label == "Both")),
                      ", Only true: ", length(which(gene_df$naive_label == "True")),
                      ", Only naive: ", length(which(gene_df$naive_label == "Naive")),
                      "\nNum. genes: ", nrow(gene_df)),
       x = "True LogFC",
       y = "Estimated LogFC") +
  theme_minimal() 
ggsave(plot_1, 
       file = paste0(plot_folder, "true-vs-naive-logfc.png"),
       width = 5, height = 5, units = "in")