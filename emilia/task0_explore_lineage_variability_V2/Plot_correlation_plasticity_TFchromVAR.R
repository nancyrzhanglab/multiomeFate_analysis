
library(tidyverse)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
results_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task0_explore_lineage_variability_V2/Feature_correlations/'
figure_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/figures/MultiomeFate_V3/Fig4/'

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
            legend.position = "right",
            legend.box.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "mm"),
            legend.direction = "vertical",
            legend.key.size= unit(0.5, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"),
            strip.background=element_rect(colour="#F0F0F0",fill="#F0F0F0"),
            strip.text = element_text(face="bold")
    ))
}

dataset_colors <- c(day0 = "gray",
                    day10_CIS = "#FBD08C",
                    day10_COCL2 = "#6DC49C",
                    day10_DABTRAM = "#9D85BE",
                    week5_CIS = "#C96D29",
                    week5_COCL2 = "#0F8241",
                    week5_DABTRAM = "#623594")



# =============================================================================
# reading data
# =============================================================================
cor.obj = readRDS(paste0(results_dir, "day10_linvar_TFchromVAR_cor_df.rds"))

df.day10_DABTRAM <- cor.obj[['day10_DABTRAM_linvar_cor_df']]
df.day10_COCL2 <- cor.obj[['day10_COCL2_linvar_cor_df']]
df.day10_CIS <- cor.obj[['day10_CIS_linvar_cor_df']]


# ==============================================================================
# plotting
# ==============================================================================
df.day10_DABTRAM$p.value.adj <- p.adjust(df.day10_DABTRAM$p.value, method = "BH")
df.day10_DABTRAM$significant <- ifelse(df.day10_DABTRAM$p.value.adj < 0.05, "yes", "no")
df.day10_DABTRAM <- df.day10_DABTRAM[order(df.day10_DABTRAM$correlation, decreasing = F), ]
df.day10_DABTRAM$order <- 1:nrow(df.day10_DABTRAM)

thres.lower <- df.day10_DABTRAM %>% 
  filter(significant == 'yes') %>%
  filter(correlation < 0) %>%
  summarise(lower = max(correlation)) %>% 
  pull(lower)
  
thres.upper <- df.day10_DABTRAM %>%
  filter(significant == 'yes') %>%
  filter(correlation > 0) %>%
  summarise(upper = min(correlation)) %>% 
  pull(upper)

df.day10_DABTRAM$TF <- rownames(df.day10_DABTRAM)
tf.anno <- df.day10_DABTRAM %>% 
  tail(8)

ggplot(df.day10_DABTRAM, aes(x = order, y = correlation)) +
  geom_point(size = 0.5, color = 'gray') +
  geom_point(data = tf.anno, aes(x = order, y = correlation), size = 0.5, color = 'red') +
  ggrepel::geom_text_repel(data = tf.anno, 
                           aes(label = TF), 
                           size = 4, 
                           nudge_y = 0.1,
                           max.overlaps = 20,
                           color = 'red') +
  geom_hline(yintercept = thres.lower, linetype = "dashed", color = "black") +
  geom_hline(yintercept = thres.upper, linetype = "dashed", color = "black") +
  ylim(-0.4, 1) +
  xlim(0, 660) +
  theme_Publication() +
  labs(x = "TFs", y = "Correlation with lineage variability") +
  theme(legend.position = "none") +
  ggtitle("Day10 DABTRAM") 
ggsave(paste0(figure_dir, "Supp_day10_DABTRAM_linvar_cor_TFchromVAR.pdf"), width = 3, height = 4)


df.day10_COCL2$p.value.adj <- p.adjust(df.day10_COCL2$p.value, method = "BH")
df.day10_COCL2$significant <- ifelse(df.day10_COCL2$p.value.adj < 0.05, "yes", "no")
df.day10_COCL2 <- df.day10_COCL2[order(df.day10_COCL2$correlation, decreasing = F), ]
df.day10_COCL2$order <- 1:nrow(df.day10_COCL2)
thres.lower <- df.day10_COCL2 %>% 
  filter(significant == 'yes') %>%
  filter(correlation < 0) %>%
  summarise(lower = max(correlation)) %>% 
  pull(lower)
thres.upper <- df.day10_COCL2 %>%
  filter(significant == 'yes') %>%
  filter(correlation > 0) %>%
  summarise(upper = min(correlation)) %>% 
  pull(upper)  
df.day10_COCL2$TF <- rownames(df.day10_COCL2)
tf.anno <- df.day10_COCL2 %>% 
  tail(8)
ggplot(df.day10_COCL2, aes(x = order, y = correlation)) +
  geom_point(size = 0.5, color = 'gray') +
  geom_point(data = tf.anno, aes(x = order, y = correlation), size = 0.5, color = 'red') +
  ggrepel::geom_text_repel(data = tf.anno, 
                           aes(label = TF), 
                           size = 4, 
                           nudge_y = 0.1,
                           max.overlaps = 20,
                           color = 'red') +
  geom_hline(yintercept = thres.lower, linetype = "dashed", color = "black") +
  geom_hline(yintercept = thres.upper, linetype = "dashed", color = "black") +
  ylim(-0.4, 1) +
  xlim(0, 660) +
  theme_Publication() +
  labs(x = "TFs", y = "Correlation with lineage variability") +
  theme(legend.position = "none") +
  ggtitle("Day10 COCL2")
ggsave(paste0(figure_dir, "Supp_day10_COCL2_linvar_cor_TFchromVAR.pdf"), width = 3, height = 4)

df.day10_CIS$p.value.adj <- p.adjust(df.day10_CIS$p.value, method = "BH")
df.day10_CIS$significant <- ifelse(df.day10_CIS$p.value.adj < 0.05, "yes", "no")
df.day10_CIS <- df.day10_CIS[order(df.day10_CIS$correlation, decreasing = F), ]
df.day10_CIS$order <- 1:nrow(df.day10_CIS)
thres.lower <- df.day10_CIS %>% 
  filter(significant == 'yes') %>%
  filter(correlation < 0) %>%
  summarise(lower = max(correlation)) %>% 
  pull(lower)
