rm(list=ls())
library(Seurat)
library(multiomeFate)

# version locally
out_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/out/Writeup14/Writeup14_priming-setting_"
plot_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/git/multiomeFate_analysis/fig/kevin/Writeup17b/Writeup17b_priming-setting_"
csv_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/out/Writeup14/"

func_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/git/multiomeFate_analysis/kevin/Writeup14_simulation/"
load(paste0(out_folder, "simulation_v2.RData"))
load(paste0(out_folder, "v2_seurat_CoSpar_prepared.RData"))

############

# first the "true" fates
cell_features <- cbind(1, simulation_res$embedding_mat)
true_coefficient_vec <- c(simulation_res$coefficient_intercept, simulation_res$coefficient_vec)
names(true_coefficient_vec) <- c("Intercept", colnames(simulation_res$embedding_mat))
fatepotential_true <- as.numeric(log10(exp(cell_features %*% true_coefficient_vec)))
names(fatepotential_true) <- rownames(cell_features)

cell_imputed_score <- rep(NA, length(Seurat::Cells(all_data)))
names(cell_imputed_score) <- Seurat::Cells(all_data)
cell_imputed_score[names(fatepotential_true)] <- fatepotential_true

plot1 <- multiomeFate:::plot_cellGrowthUmap(
  seurat_object = all_data,
  cell_imputed_score = cell_imputed_score,
  colors_use = rev(all_data@misc$fatepotential_colors),
  na_color = all_data@misc$fatepotential_na_colors,
  reduction = "ft.COCL2.umap",
  title = "Priming setting: True fate potential"
)

ggplot2::ggsave(filename = paste0(plot_folder, "Writeup17b_priming_fatepotential-true.png"),
                plot1, width = 5, height = 4.5, units = "in")
ggplot2::ggsave(filename = paste0(plot_folder, "Writeup17b_priming_fatepotential-true.pdf"),
                plot1, width = 5, height = 4.5, units = "in")


# then the "estimated" fates
cell_imputed_score <- rep(NA, length(Seurat::Cells(all_data)))
names(cell_imputed_score) <- Seurat::Cells(all_data)
cell_imputed_score[names(final_fit$cell_imputed_score)] <- final_fit$cell_imputed_score

plot1 <- multiomeFate:::plot_cellGrowthUmap(
  seurat_object = all_data,
  cell_imputed_score = cell_imputed_score,
  colors_use = rev(all_data@misc$fatepotential_colors),
  na_color = all_data@misc$fatepotential_na_colors,
  reduction = "ft.COCL2.umap",
  title = "Priming setting: Estimated fate potential"
)

ggplot2::ggsave(filename = paste0(plot_folder, "Writeup17b_priming_fatepotential-estimate.png"),
                plot1, width = 5, height = 4.5, units = "in")
ggplot2::ggsave(filename = paste0(plot_folder, "Writeup17b_priming_fatepotential-estimate.pdf"),
                plot1, width = 5, height = 4.5, units = "in")


# tne the cospar bias
# add cospar
df <- read.csv(paste0(csv_folder, "simulation_priming_cospar_obs_v2.csv"))
rownames(df) <- df$X
cospar_fate <- df$fate_bias_intraclone_transition_map_High.Low
names(cospar_fate) <- rownames(df)
cell_imputed_score <- rep(NA, length(Seurat::Cells(all_data)))
names(cell_imputed_score) <- Seurat::Cells(all_data)
cell_imputed_score[names(cospar_fate)] <- cospar_fate
cell_imputed_score[which(all_data$dataset != "day10_COCL2")] <- NA

plot1 <- multiomeFate:::plot_cellGrowthUmap(
  seurat_object = all_data,
  cell_imputed_score = cell_imputed_score,
  colors_use = rev(all_data@misc$fatepotential_colors),
  na_color = all_data@misc$fatepotential_na_colors,
  reduction = "ft.COCL2.umap",
  title = "Priming setting: CoSPAR"
)

ggplot2::ggsave(filename = paste0(plot_folder, "Writeup17b_priming_cospar.png"),
                plot1, width = 5, height = 4.5, units = "in")
ggplot2::ggsave(filename = paste0(plot_folder, "Writeup17b_priming_cospar.pdf"),
                plot1, width = 5, height = 4.5, units = "in")

# Finally, case-control
winner_lineages <- names(sort(simulation_res$lineage_future_size, decreasing = TRUE)[1:15])
loser_lineages <- setdiff(names(simulation_res$lineage_future_size), winner_lineages)

winner_cells <- rownames(simulation_res$embedding_mat)[which(simulation_res$lineage_assignment %in% winner_lineages)]
loser_cells <- rownames(simulation_res$embedding_mat)[which(simulation_res$lineage_assignment %in% loser_lineages)]

status_vec <- rep("None", length(Seurat::Cells(all_data)))
names(status_vec) <- Seurat::Cells(all_data)
status_vec[winner_cells] <- "Winner"
status_vec[loser_cells] <- "Loser"
all_data$status <- status_vec

plot1 <- scCustomize::DimPlot_scCustom(
  seurat_object = all_data, 
  group.by = "status",
  colors_use = c(Winner = "darkorange", Loser = "chartreuse4", None = "gray80"),
  reduction = "ft.COCL2.umap",
)

ggplot2::ggsave(filename = paste0(plot_folder, "Writeup17b_priming_true-status.png"),
                plot1, width = 5, height = 4.5, units = "in")
ggplot2::ggsave(filename = paste0(plot_folder, "Writeup17b_priming_true-status.pdf"),
                plot1, width = 5, height = 4.5, units = "in")

