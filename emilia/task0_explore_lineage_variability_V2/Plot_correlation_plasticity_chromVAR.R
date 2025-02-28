library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

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
results_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/Archive_Sydney/task0_explore_lineage_variability/features_associated_with_plasticity/'
sample_name <- 'day10_DABTRAM'
load(paste0(results_dir, sample_name, '_motif_mean_chromVAR_day10_plasticity_byRNA_correlation.RData'))

# ==============================================================================
# Wrangle data
# ==============================================================================

cor_vec$neg_log10_pval <- (-1) * log10(cor_vec$p.value)
cor_vec <- cor_vec %>%
  arrange(desc(correlation))
cor_vec$order <- seq(1, nrow(cor_vec))
cor_vec <- cor_vec %>% drop_na()

# ==============================================================================
# Plotting
# ==============================================================================

cor_vec$feature <- rownames(cor_vec)
cor_vec[, c('Num', 'feature_short')] <- str_split_fixed(cor_vec$feature, '::', 2)
ap1_related <- cor_vec %>% filter(str_detect(feature, 'JUN')) %>% select(feature) %>% pull()
ap1_related <- c(ap1_related, cor_vec %>% filter(str_detect(feature, 'FOS')) %>% select(feature) %>% pull())
sox_related <- cor_vec %>% filter(str_detect(feature, 'SOX')) %>% select(feature) %>% pull()
tead_related <- cor_vec %>% filter(str_detect(feature, 'TEAD')) %>% select(feature) %>% pull()
mitf <- cor_vec %>% filter(str_detect(feature, 'MITF')) %>% select(feature) %>% pull()
snai_related <- cor_vec %>% filter(str_detect(feature, 'SNAI')) %>% select(feature) %>% pull()

features_to_label_top <- head(cor_vec, 4) %>% select(feature) %>% pull()
features_to_label_top <- c(features_to_label_top, 'TEAD1')

features_to_label_bottom <- c('MITF')
features_to_label_bottom <- c(features_to_label_bottom, 'SOX10')
features_to_label_bottom <- c(features_to_label_bottom, 'SNAI2')

ggplot(cor_vec, aes(x = order, y = correlation)) +
  geom_point(size=0.5) +
  geom_point(data = subset(cor_vec, feature %in% ap1_related), size=2, color = '#0079FF') +
  geom_point(data = subset(cor_vec, feature %in% sox_related), size=2, color = '#5DB996') +
  geom_point(data = subset(cor_vec, feature %in% tead_related), size=2, color = '#F0BB78') +
  geom_point(data = subset(cor_vec, feature %in% mitf), size=2, color = '#F72C5B') +
  geom_point(data = subset(cor_vec, feature %in% snai_related), size=2, color = '#FFA07A') +
  ggrepel::geom_text_repel(data = subset(cor_vec, feature %in% features_to_label_top),
                           ggplot2::aes(label = Num),
                           box.padding = ggplot2::unit(0.4, 'lines'),
                           point.padding = ggplot2::unit(0.2, 'lines'),
                           nudge_y = 0.3,
                           nudge_x = 0.1,
                           color = 'black',
                           size=4,
                           max.overlaps = 80) +
  ggrepel::geom_text_repel(data = subset(cor_vec, feature %in% features_to_label_bottom),
                           ggplot2::aes(label = Num),
                           box.padding = ggplot2::unit(0.4, 'lines'),
                           point.padding = ggplot2::unit(0.2, 'lines'),
                           nudge_y = 0.5,
                           nudge_x = -2,
                           color = 'black',
                           size=4,
                           max.overlaps = 80) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  ylab('Correlation') +
  ggtitle(sample_name) +
  ylim(-0.5, 1.1) +
  xlim(-0.5, 645) +
  theme_Publication()
ggsave(paste0(results_dir, sample_name, '_motif_mean_chromVAR_day10_plasticity_byRNA_correlation.png'), width = 4, height = 5, units = 'in', dpi = 300)
