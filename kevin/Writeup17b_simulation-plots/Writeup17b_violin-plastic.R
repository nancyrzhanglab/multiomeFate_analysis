rm(list=ls())
library(Seurat)
library(multiomeFate)
library(ggplot2)

# version on HPC3
# out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup14/Writeup14_plastic-setting_"
# plot_folder <- "~/project/Multiome_fate/git/multiomeFate_analysis_kevin/fig/kevin/Writeup14/Writeup14_plastic-setting_"
# func_folder <- "~/project/Multiome_fate/git/multiomeFate_analysis_kevin/kevin/Writeup14_simulation/"
# load(paste0(out_folder, "simulation_v2.RData"))

# version locally
out_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/out/Writeup14/Writeup14_plastic-setting_"
plot_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/git/multiomeFate_analysis/fig/kevin/Writeup17b/Writeup17b_plastic-setting_"
func_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/git/multiomeFate_analysis/kevin/Writeup14_simulation/"

load(paste0(out_folder, "simulation.RData"))
source(paste0(func_folder, "func_seurat.R"))

cell_imputed_score <- as.numeric(simulation_res$embedding_mat %*% simulation_res$coefficient_vec) + simulation_res$coefficient_intercept
cell_imputed_score <- log10(exp(cell_imputed_score))
names(cell_imputed_score) <- rownames(simulation_res$embedding_mat)

assigned_lineage <- simulation_res$lineage_assignment
names(assigned_lineage) <- rownames(simulation_res$embedding_mat)

lineage_vec <- assigned_lineage[names(cell_imputed_score)]
tab_vec <- table(assigned_lineage)

# compute the median for each lineage
median_vec <- sapply(unique(lineage_vec), function(lineage){
  idx <- which(lineage_vec == lineage)
  stats::median(cell_imputed_score[idx])
})
names(median_vec) <- unique(lineage_vec)

lineage_sizes <- simulation_res$lineage_future_size
lineage_names <- names(lineage_sizes)[order(lineage_sizes, decreasing = TRUE)[1:20]]

cell_idx <- which(lineage_vec %in% lineage_names)
median_vec <- median_vec[lineage_names]
lineage_sizes <- lineage_sizes[lineage_names]

# Loop through and adjust duplicates minimally
seen_values <- c()  # Store already seen values
for (i in length(lineage_sizes):1) {
  while (lineage_sizes[i] %in% seen_values) {
    lineage_sizes[i] <- lineage_sizes[i] + 1  # Smallest possible increment
  }
  seen_values <- c(seen_values, lineage_sizes[i])  # Update seen values
}

lineage_labels <- setNames(paste0("(", lineage_sizes, ")"), names(lineage_sizes))

# form data frame
df <- data.frame(lineage = lineage_labels[as.character(lineage_vec[cell_idx])],
                 imputed_count = cell_imputed_score[cell_idx])
df

col_vec <- rep("lightgray", length(lineage_names))
names(col_vec) <- lineage_names

plot1 <- ggplot2::ggplot(df, ggplot2::aes(x=lineage, y=imputed_count))
plot1 <- plot1 + ggplot2::geom_violin(trim=T, scale = "width", ggplot2::aes(fill=lineage))
plot1 <- plot1 + ggplot2::scale_fill_manual(values = col_vec) 
plot1 <- plot1 + ggplot2::geom_jitter(shape=16, 
                                      position=ggplot2::position_jitter(0.2), alpha = 0.1, size = 0.5)
plot1 <- plot1 + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::scale_x_discrete(limits = lineage_labels,
                                           guide = ggplot2::guide_axis(angle = 45))
plot1 <- plot1 + ggplot2::stat_summary(fun = median, geom = "crossbar", 
                                       width = 0.75, color = "#0D8242")

ggplot2::ggsave(filename = paste0(plot_folder, "fatepotential_true-violinplot.png"),
                plot1, device = "png", width = 6, height = 3, units = "in")

########

plot1 <- ggplot2::ggplot(df, ggplot2::aes(x = lineage, y = imputed_count, fill = lineage)) +
  geom_violin(trim = TRUE, scale = "width", color = "black") +  # Black border for clarity
  scale_fill_manual(values = col_vec) +
  stat_summary(fun = median, geom = "crossbar", width = 0.75, color = "#0D8242") +  # Median line
  scale_x_discrete(limits = lineage_labels, guide = guide_axis(angle = 45)) +  # Rotated x-axis labels
  labs(y = NULL, x = NULL) +  # No axis titles
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    panel.grid.major.x = element_blank(),  # Remove vertical grid lines
    panel.grid.minor = element_blank(),    # Remove minor grid lines
    panel.grid.major.y = element_line(color = "gray80", linetype = "dashed")
  ) +
  coord_cartesian(ylim = c(-1, 2.7))  # Keep y-axis range

plot1 <- plot1 + ggplot2::labs(y = NULL) + ggplot2::theme(axis.title = ggplot2::element_blank()) 

ggplot2::ggsave(filename = paste0(plot_folder, "fatepotential_true-violinplot_cleaned.png"),
                plot1, device = "png", width = 4, height = 1.5, units = "in", 
                dpi = 450)
