rm(list=ls())
library(Seurat)
library(Signac)
library(multiomeFate)
library(ggplot2)
library(Ckmeans.1d.dp)

out_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/out/Writeup14/Writeup14_priming-setting_"
plot_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/git/multiomeFate_analysis/fig/kevin/Writeup14/tmp_Writeup14_priming-setting_"
func_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/git/multiomeFate_analysis/kevin/Writeup14_simulation/"

source(paste0(func_folder, "func_seurat.R"))

load(paste0(out_folder, "simulation_v2.RData"))
load(paste0(out_folder, "all_data_fasttopics.RData"))

names(simulation_res$lineage_assignment) <- rownames(simulation_res$embedding_mat)
simulation_data <- .form_simulation_seurat_fate(final_fit = final_fit,
                                                simulation_res = simulation_res)

#########

print("Plotting violin plots")

cell_imputed_score <- simulation_data@meta.data[,"fatepotential"]
names(cell_imputed_score) <- Seurat::Cells(simulation_data)

plot1 <- multiomeFate:::.plot_anova_helper(seurat_object = simulation_data,
                                           cell_imputed_score = cell_imputed_score,
                                           lineage_future_size = simulation_data@misc$lineage_observed_count,
                                           assigned_lineage_variable = "assigned_lineage",
                                           ylab = paste0("Growth potential"),
                                           ylim = c(max(stats::quantile(cell_imputed_score, na.rm = TRUE, prob = 0.05), -5), 
                                                    NA))
plot1

################

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

hist(abs(gene_df$true_logfc), breaks = 25)
rug(abs(gene_df$true_logfc)[which(true_gene_vec)], col = 2, lwd = 2)

##########

# cluster_res <- Ckmeans.1d.dp::Ckmeans.1d.dp(
#   x = log10(simulation_data@misc[["lineage_observed_count"]]+1),
#   k = 4
# )
# 
# tmp <- log10(simulation_data@misc[["lineage_observed_count"]]+1)
# hist(tmp)
# rug(tmp[which(cluster_res$cluster == which.max(cluster_res$centers))], col = 2, lwd = 2)

# do a DE analysis between the winner lineages vs loser lineages
lineage_names <- names(simulation_data@misc[["lineage_observed_count"]])
winner_lineages <- lineage_names[which(cluster_res$cluster == which.max(cluster_res$centers))]
loser_lineages <- lineage_names[which(cluster_res$cluster == which.min(cluster_res$centers))]

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
plot_1
