# based on https://github.com/nancyrzhanglab/multiomeFate_analysis/blob/emilia/emilia/Final/Fig2/Show_in_vitro_heterogneity.R#L242

rm(list = ls())

library(multiomeFate)
library(Seurat)
library(UCell)
library(tidyverse)
library(data.table)
library(ggplot2)
library(ggpubr)
library(gridExtra)

plot_folder <- "~/project/Multiome_fate/git/multiomeFate_analysis_kevin/fig/kevin/Writeup18/"
results_dir <- "~/project/Multiome_fate/out/emilia/task2_correlate_fate_potential_and_features_V2/"

all_data <- multiomeFate:::data_loader(which_files = c("fasttopics"))
scores <- read.csv(paste0(results_dir, 'MITF_AXL_UCell_scores.csv'))
rownames(scores) <- scores$X; scores <- scores[,-1]
stopifnot(all(rownames(scores) == Seurat::Cells(all_data)))

all_data$MITF_Program <- scores[,"MITF_Program"]
all_data$AXL_Program <- scores[,"AXL_Program"]

############

treatment_vec <- c("CIS", "COCL2", "DABTRAM")
for(treatment in treatment_vec){
  all_data_subset <- subset(all_data, dataset %in% paste0("day10_", treatment))
  df <- all_data_subset@meta.data[,c("dataset", "assigned_lineage", "AXL_Program")]
  df$AXL_Program <- scale(df$AXL_Program)
  df_lin <- df %>% 
    group_by(assigned_lineage) %>% 
    summarise(mean = mean(AXL_Program),
              var = var(AXL_Program),
              lin_size = n()) %>% 
    filter(lin_size > 10)
  df_lin <- df_lin[order(df_lin$var), ]
  lins_to_plot <- c(head(df_lin, 4)$assigned_lineage,
                    tail(df_lin, 4)$assigned_lineage)
  
  
  df_to_plot <- df[df$assigned_lineage %in% lins_to_plot, ]
  plot1 <- ggplot2::ggplot(df_to_plot, ggplot2::aes(x = reorder(assigned_lineage, -AXL_Program, median), y = AXL_Program)) +
    ggplot2::geom_violin(scale = 'width')+
    ggplot2::geom_boxplot(width = 0.2, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.1, alpha = 0.6) +
    ggplot2::ylab('AXL Program') +
    ggplot2::theme_classic() +
    ggplot2::xlab('') +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   legend.position = 'none') +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 6))
  ggplot2::ggsave(plot1, filename = paste0(plot_folder, "Writeup18_violin_axl_", treatment, ".png"),
                  height = 8, width = 5)
  
  ###############
  
  # Step 1: Compute median per lineage
  medians <- df_to_plot %>%
    group_by(assigned_lineage) %>%
    summarize(median_axl = median(AXL_Program))
  medians$median_axl <- pmin(pmax(medians$median_axl, -2), 2)
  
  # Step 2: Merge into original data
  df_to_plot <- df_to_plot %>%
    left_join(medians, by = "assigned_lineage") %>%
    mutate(assigned_lineage = factor(
      assigned_lineage,
      levels = medians$assigned_lineage[order(-medians$median_axl)]
    ))
  
  # Step 3: Create a continuous color scale using a custom gradient palette
  # We define a gradient from -2 (purple) to 0 (bisque) to 2 (orange)
  gradient_palette <- scales::col_numeric(
    palette = c("#604CC3", "bisque", "#FFA500"),
    domain = c(-2, 2),
    na.color = "grey80"
  )
  
  # Step 4: Apply violin fill colors based on the lineage median
  fill_colors <- setNames(gradient_palette(medians$median_axl), medians$assigned_lineage)
  
  # Step 5: Plot
  plot1 <- ggplot(df_to_plot, aes(x = assigned_lineage, y = AXL_Program, fill = assigned_lineage)) +
    geom_violin(scale = 'width', width = 0.9, color = "gray80", linewidth = 0.2) +
    geom_boxplot(width = 0.5, outlier.shape = NA, color = "black", fill = "white", size = 0.3) +
    geom_jitter(width = 0.1, size = 0.8, alpha = 0.5, color = "black") +
    scale_fill_manual(values = fill_colors) +
    ylab('AXL Program') +
    xlab('') +
    theme_classic(base_size = 8) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_text(size = 7),
      axis.title.y = element_text(size = 8),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none",
      plot.margin = margin(5, 5, 5, 5)
    )
  plot1 <- plot1 + stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black", fatten = 0)
  
  ggplot2::ggsave(plot1, filename = paste0(plot_folder, "Writeup18_violin_axl_", treatment, "_cleaned.png"),
                  height = 1.5, width = 3)
}

