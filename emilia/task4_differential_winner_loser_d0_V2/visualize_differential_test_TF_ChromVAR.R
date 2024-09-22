library(ggplot2)

in_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_differential_winner_loser_d0_V2/'

theme_Publication<- function(base_size=12, base_family="sans") {
  library(grid)
  library(ggthemes)
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
            strip.background=element_rect(colour="#F0F0F0",fill="#F0F0F0"),
            strip.text = element_text(face="bold")
    ))
}
# ==============================================================================
# Read data
# ==============================================================================

load(paste0(in_dir, 'differential_winner_loser_chromVAR_lineage_specific_adapation_TF_DABTRAM_t_test_results.RData'))
diff_res_dabtram <- t_test_results

load(paste0(in_dir, 'differential_winner_loser_chromVAR_lineage_specific_adaptation_TF_COCL2_t_test_results.RData'))
diff_res_cocl2 <- t_test_results

load(paste0(in_dir, 'differential_winner_loser_chromVAR_lineage_specific_adaptation_TF_CIS_t_test_results.RData'))
diff_res_cis <- t_test_results

# ==============================================================================
# Wrangle data
# ==============================================================================

## helper
wrangle_data <- function(diff_res) {
  diff_res$p_val <- as.numeric(diff_res$p_val)
  diff_res$p_adj <- p.adjust(diff_res$p_val, method = 'BH')
  diff_res$neg_log10_p_val <- -log10(diff_res$p_val)
  
  diff_res$mean_winner <- as.numeric(diff_res$mean_winner)
  diff_res$mean_other <- as.numeric(diff_res$mean_other)
  # diff_res$fold_change <- log2(diff_res$mean_winner / diff_res$mean_other)
  diff_res$fold_change <- diff_res$mean_winner - diff_res$mean_other
  
  # order by fold_change in decreasing order
  diff_res <- diff_res[order(diff_res$fold_change, decreasing = TRUE),]
  
  return(diff_res)
}


diff_res_dabtram <- wrangle_data(diff_res_dabtram)
diff_res_cocl2 <- wrangle_data(diff_res_cocl2)
diff_res_cis <- wrangle_data(diff_res_cis)

# ==============================================================================
# Plot
# ==============================================================================
# Volcano plot for DABTRAM
diff_res_dabtram$neg_log10_p_val <- ifelse(diff_res_dabtram$neg_log10_p_val > 100, 100, diff_res_dabtram$neg_log10_p_val)

up_10_dabtram <- head(diff_res_dabtram, 10)
bottom_10_dabtram <- tail(diff_res_dabtram, 10)
thres_dabtram <- min(diff_res_dabtram[diff_res_dabtram$p_adj < 0.05,]$neg_log10_p_val)

ggplot(diff_res_dabtram, aes(x = fold_change, y = neg_log10_p_val)) +
  geom_point() +
  geom_point(data = up_10_dabtram, aes(x = fold_change, y = neg_log10_p_val), color = 'red') +
  geom_point(data = bottom_10_dabtram, aes(x = fold_change, y = neg_log10_p_val), color = 'blue') +
  ggrepel::geom_text_repel(data = up_10_dabtram, aes(label = feature), color = 'red', nudge_x = 0.2, nudge_y = 0.5) +
  ggrepel::geom_text_repel(data = bottom_10_dabtram, aes(label = feature), color = 'blue',nudge_x = -0.2,  nudge_y = 0.5, max.overlaps = 50) +
  geom_hline(yintercept = thres_dabtram, linetype = 'dashed', color = 'gray') +
  labs(x = 'diff (winner - loser)', y = '-log10 p-value', title = 'day 0 DABTRAM') +
  xlim(-3, 3) +
  theme_Publication()

ggsave(paste0(in_dir, 'differential_winner_loser_saver_lineage_specific_adaptation_TF_DABTRAM_volcano_plot.png'),dpi = 300, width = 6, height = 6)

# Volcano plot for COCL2
diff_res_cocl2$neg_log10_p_val <- ifelse(diff_res_cocl2$neg_log10_p_val > 100, 100, diff_res_cocl2$neg_log10_p_val)

up_10_cocl2<- head(diff_res_cocl2, 10)
bottom_10_cocl2 <- tail(diff_res_cocl2, 10)
thres_cocl2 <- min(diff_res_cocl2[diff_res_cocl2$p_adj < 0.05,]$neg_log10_p_val)

ggplot(diff_res_cocl2, aes(x = fold_change, y = neg_log10_p_val)) +
  geom_point() +
  geom_point(data = up_10_cocl2, aes(x = fold_change, y = neg_log10_p_val), color = 'red') +
  geom_point(data = bottom_10_cocl2, aes(x = fold_change, y = neg_log10_p_val), color = 'blue') +
  ggrepel::geom_text_repel(data = up_10_cocl2, aes(label = feature), color = 'red', nudge_x = 0.05, nudge_y = 0.5) +
  ggrepel::geom_text_repel(data = bottom_10_cocl2, aes(label = feature), color = 'blue',nudge_x = -0.05,  nudge_y = 0.5, max.overlaps = 50) +
  geom_hline(yintercept = thres_cocl2, linetype = 'dashed', color = 'gray') +
  labs(x = 'diff (winner - loser)', y = '-log10 p-value', title = 'day 0 COCL2') +
  xlim(-1.8, 1.8) +
  theme_Publication()

ggsave(paste0(in_dir, 'differential_winner_loser_saver_lineage_specific_adaptation_TF_COCL2_volcano_plot.png'),dpi = 300, width = 6, height = 6)

# Volcano plot for CIS
diff_res_cis$neg_log10_p_val <- ifelse(diff_res_cis$neg_log10_p_val > 100, 100, diff_res_cis$neg_log10_p_val)

up_10_cis<- head(diff_res_cis, 10)
bottom_10_cis <- tail(diff_res_cis, 10)
thres_cis <- min(diff_res_cis[diff_res_cis$p_adj < 0.05,]$neg_log10_p_val)

ggplot(diff_res_cis, aes(x = fold_change, y = neg_log10_p_val)) +
  geom_point() +
  geom_point(data = up_10_cis, aes(x = fold_change, y = neg_log10_p_val), color = 'red') +
  geom_point(data = bottom_10_cis, aes(x = fold_change, y = neg_log10_p_val), color = 'blue') +
  ggrepel::geom_text_repel(data = up_10_cis, aes(label = feature), color = 'red', nudge_x = 0.05, nudge_y = 0.5) +
  ggrepel::geom_text_repel(data = bottom_10_cis, aes(label = feature), color = 'blue',nudge_x = -0.05,  nudge_y = 0.5) +
  geom_hline(yintercept = thres_cis, linetype = 'dashed', color = 'gray') +
  labs(x = 'diff (winner - loser)', y = '-log10 p-value', title = 'day 0 CIS') +
  # xlim(-0.8, 0.8) +
  theme_Publication()

ggsave(paste0(in_dir, 'differential_winner_loser_saver_lineage_specific_adaptation_TF_CIS_volcano_plot.png'),dpi = 300, width = 6, height = 6)
