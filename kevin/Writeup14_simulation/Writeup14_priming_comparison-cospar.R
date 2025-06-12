rm(list=ls())

library(Seurat)
library(multiomeFate)
library(ggplot2)
library(Ckmeans.1d.dp)

setting <- "priming"
out_folder <- paste0("~/project/Multiome_fate/out/kevin/Writeup14/Writeup14_", setting, "-setting_")
plot_folder <- paste0("~/project/Multiome_fate/git/multiomeFate_analysis_kevin/fig/kevin/Writeup14/Writeup14_", setting, "-setting_")
func_folder <- "~/project/Multiome_fate/git/multiomeFate_analysis_kevin/kevin/Writeup14_simulation/"
csv_folder <- "~/project/Multiome_fate/git/multiomeFate_analysis_kevin/csv/kevin/Writeup14/"

df <- read.csv(paste0(csv_folder, "sijia_simulation_priming_cospar_obs.csv"))
rownames(df) <- df$X

day10_idx <- which(df$dataset == "day10_COCL2")
fatepotential_true <- df$fatepotential_true[day10_idx]
cospar_bias <- df$fate_bias_intraclone_transition_map_High.Low[day10_idx]
cor_val <- stats::cor(fatepotential_true, cospar_bias)

png(paste0(plot_folder, "cell_true-vs-cospar.png"))
plot(x = fatepotential_true,
     y = cospar_bias,
     xlab = "Cell's true fate potential",
     ylab = "CoSpar's estimated fate bias for High vs. Low",
     main = paste0("Setting ", setting, ": Corr=", round(cor_val,2)))
graphics.off()

##################

source(paste0(func_folder, "func_seurat.R"))

all_data <- multiomeFate:::data_loader(which_files = "fasttopics")

load(paste0(out_folder, "simulation.RData"))
names(simulation_res$lineage_assignment) <- rownames(simulation_res$embedding_mat)
simulation_data <- .form_simulation_seurat_fate(final_fit = final_fit,
                                                simulation_res = simulation_res)

# compute the "denoised" gene expression profiles
all_data <- subset(all_data, dataset == "day10_COCL2")
cell_embedding <- all_data[["fasttopic.COCL2"]]@cell.embeddings
gene_embedding <- all_data[["fasttopic.COCL2"]]@feature.loadings
denoised_mat <- tcrossprod(cell_embedding, gene_embedding)
denoised_mat <- denoised_mat[Seurat::Cells(simulation_data),]
p <- ncol(denoised_mat)

# compute the A: {true correlation between each gene and the fate potential}, with p-values
fatepotential_true <- simulation_data$fatepotential_true
cor_mat <- t(sapply(1:ncol(denoised_mat), function(j){
  tmp <- stats::cor.test(x = denoised_mat[,j],
                         y = fatepotential_true)
  
  c(correlation = as.numeric(tmp$estimate), p.value = as.numeric(tmp$p.value))
}))
rownames(cor_mat) <- colnames(denoised_mat)
cor_df <- as.data.frame(cor_mat)
colnames(cor_df) <- c("true_correlation", "true_p.value")
cor_df$true_p.value.adj <- stats::p.adjust(cor_df$true_p.value, method = "BH")
true_gene_vec <- rep(TRUE, nrow(cor_df))
true_gene_vec[which(cor_df$true_p.value.adj >= 0.05)] <- FALSE
true_gene_vec[which(abs(cor_df$true_correlation) <= 0.3)] <- FALSE
cor_df$true_gene <- true_gene_vec

# compute the B: {estimated correlation between each gene and the estimated cospar bias}, with p-values
fatepotential <- df$fate_bias_intraclone_transition_map_High.Low
names(fatepotential) <- rownames(df)
fatepotential <- fatepotential[rownames(denoised_mat)]
cor_mat_est <- t(sapply(1:ncol(denoised_mat), function(j){
  tmp <- stats::cor.test(x = denoised_mat[,j],
                         y = fatepotential)
  
  c(correlation = as.numeric(tmp$estimate), p.value = tmp$p.value)
}))
cor_df$est_correlation <- cor_mat_est[,"correlation"]
cor_df$est_p.value <- cor_mat_est[,"p.value"]
cor_df$est_p.value.adj <- stats::p.adjust(cor_df$est_p.value, method = "BH")
est_gene_vec <- rep(TRUE, nrow(cor_df))
est_gene_vec[which(cor_df$est_p.value.adj >= 0.05)] <- FALSE
est_gene_vec[which(abs(cor_df$est_correlation) <= 0.15)] <- FALSE
cor_df$est_gene <- est_gene_vec

# compute the correlation between {A} and {B}, and make the scatterplot
cor_val <- stats::cor(cor_df$true_correlation, 
                      cor_df$est_correlation)

plot_1 <- ggplot(cor_df, aes(x = true_correlation, y = est_correlation)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = "gray") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#2f2f2f") +
  labs(title = sprintf("True-vs-estimated correlation: \nCorrelation: %.2f", cor_val),
       x = "Correlation with true fate potential",
       y = "Correlation with estimated fate potential") +
  theme_minimal() + coord_equal()
ggsave(plot_1, 
       file = paste0(plot_folder, "true-vs-cospar-correlation.png"),
       width = 6, height = 6, units = "in")

###############

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

# Now do this for the estimated fate potential
day10_idx <- which(df$dataset == "day10_COCL2")
cospar_bias <- df$fate_bias_intraclone_transition_map_High.Low[day10_idx]

winner_cells <- which(cospar_bias < 0.5) # this is just a parameterization thing
loser_cells <- which(cospar_bias > 0.5)
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
gene_df$est_logfc <- wilcox_results[,"logfc"]
gene_df$est_p.value <- wilcox_results[,"p.value"]
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

# make a correlation logfc scatterplot
plot_1 <- ggplot(gene_df, aes(x = true_logfc, y = est_logfc, color = est_label)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = colors) +
  geom_vline(xintercept = 0, color = "#2f2f2f") +
  geom_hline(yintercept = 0, color = "#2f2f2f") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#2f2f2f") +
  labs(title = paste0("Winner-vs-loser cells based on estimated cospar bias:",
                      "\nCorrelation in true-and-cospar logfc: ", 
                      round(cor_val,2),
                      "\nBoth: ", length(which(gene_df$est_label == "Both")),
                      ", Only true: ", length(which(gene_df$est_label == "True")),
                      ", Only estimate: ", length(which(gene_df$est_label == "Estimate")),
                      "\nNum. genes: ", nrow(gene_df)),
       x = "True LogFC",
       y = "Estimated LogFC") +
  theme_minimal() 
ggsave(plot_1, 
       file = paste0(plot_folder, "true-vs-cospar-logfc.png"),
       width = 5, height = 5, units = "in")


