# https://github.com/nancyrzhanglab/multiomeFate_analysis/blob/emilia/emilia/Final/Fig2/Plot_lin_var_scores.R#L199
rm(list = ls())

set.seed(123)

library(tidyverse)
library(ggpubr)
library(ggplot2)
library(grid)
library(ggthemes)

results_dir <- '~/project/Multiome_fate/out/emilia/task0_explore_lineage_variability_V2/'
plot_folder <- "~/project/Multiome_fate/git/multiomeFate_analysis_kevin/fig/kevin/Writeup18/"

# read lineage variability scores
lin_var.day10_COCL2 <- read.csv(paste0(results_dir, 'day10_COCL2/lineage_variability_day10_COCL2_saver_sample.csv'))
lin_var.day10_COCL2_Shuffled <- read.csv(paste0(results_dir, 'day10_COCL2/lineage_variability_shuffledday10_COCL2_saver_sample.csv'))

theme_Publication<- function(base_size=12, base_family="sans") {

  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            plot.subtitle = element_text(face = "bold", hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold"),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line.x = element_line(colour="black"),
            axis.line.y = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="white"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.box.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "mm"),
            legend.direction = "horizontal",
            legend.key.size= unit(0.5, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"),
            strip.background=element_rect(colour="gray90",fill="gray90"),
            strip.text = element_text(face="bold")
    ))
}


dataset_colors <- c(day0 = "gray",
                    day10_CIS = "#FBD08C",
                    day10_COCL2 = "#6DC49C",
                    day10_DABTRAM = "#9D85BE",
                    week5_CIS = "#C96D29",
                    week5_COCL2 = "#0F8241",
                    week5_DABTRAM = "#623594",
                    shuffled.day0 = 'gray90',
                    shuffled.day10_DABTRAM = 'gray900',
                    shuffled.day10_COCL2 = 'gray90',
                    shuffled.day10_CIS = 'gray90',
                    shuffled.week5_DABTRAM = 'gray90',
                    shuffled.week5_COCL2 = 'gray90',
                    shuffled.week5_CIS = 'gray90')

lin_var.day10_COCL2$dataset <- 'day10_COCL2'
lin_var.day10_COCL2$category <- 'day10_COCL2'
lin_var.day10_COCL2_Shuffled$dataset <- 'day10_COCL2'
lin_var.day10_COCL2_Shuffled$category <- 'shuffled.day10_COCL2'


df.day10_COCL2 <- rbind(lin_var.day10_COCL2, lin_var.day10_COCL2_Shuffled)
df.day10_COCL2$category <- factor(df.day10_COCL2$category, levels = c('shuffled', 'day10_COCL2'))

df1 <- rbind(lin_var.day10_COCL2, lin_var.day10_COCL2_Shuffled)
pvalue_res <- stats::wilcox.test(
  x = df1$normalized_avg_eud_dist_by_shuffle[which(df1$category == "day10_COCL2")],
  y = df1$normalized_avg_eud_dist_by_shuffle[which(df1$category == "shuffled.day10_COCL2")]
)

# plot COCL2 violin plot
df1 <- rbind(lin_var.day10_COCL2, lin_var.day10_COCL2_Shuffled)
plot1 <- ggplot(df1, aes(x = dataset, y= normalized_avg_eud_dist_by_shuffle)) +
  geom_violin(aes(fill = category), scale = 'width', width = 0.8) +
  geom_boxplot(aes(group = category), width = 0.2, position = position_dodge(0.8), fill = 'white', outlier.shape = NA) +
  scale_fill_manual(values = dataset_colors) +
  ylab('Lineage variability (RNA)') +
  theme_Publication() +
  ggplot2::ggtitle(paste0("Wilcoxon p-value: ", pvalue_res$p.value))
ggsave(plot1, file = paste0(plot_folder, 'Writeup18_plasticity_violin_scores.png'), width = 8, height = 3)

##########

# Step 5: Plot
plot1 <- ggplot(df1, aes(x = category, y = normalized_avg_eud_dist_by_shuffle, fill = category)) +
  geom_violin(scale = 'width', width = 0.9, color = "gray70", linewidth = 0.2) +
  geom_boxplot(width = 0.5, outlier.shape = NA, color = "black", fill = "white", size = 0.3) +
  scale_fill_manual(values = dataset_colors) +
  ylab('Lineage variability (RNA)') +
  xlab('') +
  theme_classic(base_size = 8) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 7),
    axis.title.y = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5)
  ) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black", fatten = 0)

ggplot2::ggsave(plot1, filename = paste0(plot_folder, "Writeup18_plasticity_violin_scores_cleaned.png"),
                width = 2, height = 1)
