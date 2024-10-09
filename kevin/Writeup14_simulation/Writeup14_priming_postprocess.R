rm(list=ls())
library(Seurat)
library(multiomeFate)

out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup14/Writeup14_priming-setting_"
plot_folder <- "~/project/Multiome_fate/git/multiomeFate_analysis_kevin/fig/kevin/Writeup14/Writeup14_priming-setting_"
func_folder <- "~/project/Multiome_fate/git/multiomeFate_analysis_kevin/kevin/Writeup14_simulation/"

source(paste0(func_folder, "func_seurat.R"))

load(paste0(out_folder, "simulation.RData"))

names(simulation_res$lineage_assignment) <- rownames(simulation_res$embedding_mat)
all_data <- .form_simulation_seurat_fate(final_fit = final_fit,
                                         simulation_res = simulation_res)

#########

print("Plotting violin plots")

cell_imputed_score <- all_data@meta.data[,"fatepotential"]
names(cell_imputed_score) <- Seurat::Cells(all_data)

plot1 <- multiomeFate:::.plot_anova_helper(seurat_object = all_data,
                            cell_imputed_score = cell_imputed_score,
                            lineage_future_size = all_data@misc$lineage_observed_count,
                            assigned_lineage_variable = "assigned_lineage",
                            ylab = paste0("Growth potential"),
                            ylim = c(max(stats::quantile(cell_imputed_score, na.rm = TRUE, prob = 0.05), -5), 
                                     NA))
ggplot2::ggsave(filename = paste0(plot_folder, "fatepotential-violinplot.png"),
                plot1, device = "png", width = 6, height = 3, units = "in")


#######

print("Plotting lineage scatterplot")

lineage_imputed_count <- all_data@misc$lineage_observed_count
lineage_future_count <- all_data@misc$lineage_imputed_count
lineage_future_count <- lineage_future_count[names(lineage_imputed_count)]

plot1 <- multiomeFate:::plot_lineageScatterplot(
  lineage_future_count = lineage_future_count,
  lineage_imputed_count = lineage_imputed_count,
  title = "Lineage fate potential, (Log-scale)"
)

ggplot2::ggsave(filename = paste0(plot_folder, "fatepotential-lineage_prediction.png"),
                plot1, device = "png", width = 10, height = 10, units = "in")


#########

fatepotential <- all_data$fatepotential
fatepotential_true <- all_data$fatepotential_true

plot1 <-.plot_cellFateScatterplot(
  fatepotential = fatepotential,
  fatepotential_true = fatepotential_true,
  title = "Cell fate potential, (Log-scale)"
)

ggplot2::ggsave(filename = paste0(plot_folder, "fatepotential-cell_prediction.png"),
                plot1, device = "png", width = 10, height = 10, units = "in")

