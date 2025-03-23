rm(list=ls())
library(Seurat)
library(tidyverse)
library(ggpubr)
library(ggdensity)
library(ggplot2)
library(RColorBrewer)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'

remove_unassigned_cells <- TRUE
treatment <- 'DABTRAM'

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
            axis.text = element_text(size = 16),
            # axis.ticks = element_blank(),
            axis.line.x = element_line(colour="black", linewidth = 1),
            axis.line.y = element_line(colour="black", linewidth = 1),
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

# ==============================================================================
# Read data
# ==============================================================================
load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))
fp.d0_d10 <- as.data.frame(all_data_fatepotential[[paste0("fatepotential_", treatment, "_d0_d10")]][["cell_imputed_score"]])

# fate bias from d0 to week5 adapting
df.bias <- read.csv(paste0(out_dir, 'adapting_bias_thres_0_', treatment, '.csv'))

# ==============================================================================
# Wrangle
# ==============================================================================
metadat <- all_data@meta.data
metadat$cell_id <- rownames(metadat)

colnames(fp.d0_d10) <- c("fatepotential_d0_d10")
fp.d0_d10$cell_id <- rownames(fp.d0_d10)

df <- merge(df.bias, fp.d0_d10, by = 'cell_id')
df <- merge(df, metadat[, c('cell_id', 'assigned_lineage')], by = 'cell_id')

ggplot(df, aes(x = fatepotential_d0_d10, y = bias)) +
  geom_point() +
  geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 0.5, linetype = 'dashed') +
  scale_fill_manual(values = brewer.pal(5, "YlGn")) +
  labs(title = treatment,
       x = 'Fate potential from d0 to d10',
       y = 'Fate bias from d0 to w5-adapting') +
  theme_Publication()

ggplot(df, aes(x = adaptingFP, y = nonAdaptingFP)) +
  geom_point() +
  stat_cor() +
  # geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  # scale_fill_manual(values = brewer.pal(5, "YlGn")) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'red', linewidth = 2) +
  labs(title = treatment,
       x = 'Fate potential from d0 to d10-adapting',
       y = 'Fate potential from d0 to d10-non-adapting') +
  theme_Publication()
ggsave(paste0(out_dir, 'adaptingFP_thres_0_vs_nonAdaptingFP_d0_d10_', treatment, '.png'), width = 4.5, height = 5, dpi = 300)

ggplot(df, aes(x = fatepotential_d0_d10, y = bias)) +
  geom_point(aes(color = assigned_lineage)) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 0.5, linetype = 'dashed') +
  scale_fill_manual(values = brewer.pal(5, "YlGn")) +
  labs(title = treatment,
       x = 'Fate potential from d0 to d10',
       y = 'Fate bias from d0 to w5-adapting') +
  theme_Publication()
# ggsave(paste0(out_dir, 'adapting_bias_thres_0_vs_fp_d0_d10_', treatment, '.png'), width = 4.5, height = 5, dpi = 300)

df$adapting_cell_count <- 10**df$adaptingFP
df <- df[order(df$adapting_cell_count),]
ggplot(df, aes(x = fatepotential_d0_d10, y = bias)) +
  geom_point(aes(color = adapting_cell_count)) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 0.5, linetype = 'dashed') +
  scale_size_continuous(breaks = c(0.001614558, 0.030358710, 0.051100146, 0.087332899, 1.366026551)) +
  # scale_color_gradient(low = "blue", high = "red") +
  labs(title = treatment,
       x = 'Fate potential from d0 to d10',
       y = 'Fate bias from d0 to w5-adapting') +
  theme_Publication()

df <- df[order(df$adaptingFP),]
ggplot(df, aes(x = fatepotential_d0_d10, y = adaptingFP)) +
  geom_point(aes(fill = frac_adapting), size = 2, shape = 21, color = 'darkgray') +
  scale_fill_gradient(low = "gray", high = "red") +
  stat_cor() +
  labs(title = treatment,
       x = 'Fate potential from d0 to d10',
       y = 'Fate potential from d0 to d10-adapting') +
  theme_Publication()

df$imputed_cell_count <- 10**df$fatepotential_d0_d10
df$adapting_cell_count <- 10**df$adaptingFP
df$frac_adapting <- df$adapting_cell_count / df$imputed_cell_count * 100
df$frac_adapting <- ifelse(df$frac_adapting > 100, 100, df$frac_adapting)

df <- df[order(df$frac_adapting, decreasing = F),]
ggplot(df, aes(x = fatepotential_d0_d10, y = bias)) +
  geom_point(aes(fill = adapting_cell_count, size = adapting_cell_count),  shape = 21, color = '#888888', stroke = 0.5) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 0.5, linetype = 'dashed') +
  scale_fill_gradient(low = "#DCDCDC", high = "#7C00FE") +
  # scale_size_continuous(breaks = c(0, 0.1, 0.2, 0.5)) +
  scale_size(range = c(0, 5)) +
  labs(title = treatment,
       y = 'Day0 fate bias to Day10-adapting',
       x = 'Fate potential from Day0 to Day10') +
  theme_Publication()

ggsave(paste0(figure_dir, 'adapting_bias_thres_0_vs_fp_d0_d10_DABTRAM.pdf'), width = 5, height = 3.5)

lineage.size <- df %>% 
  group_by(assigned_lineage) %>% 
  summarise(d0.lin.size = n()) %>% 
  arrange(desc(d0.lin.size))

df.1 <- df[df$assigned_lineage == 'Lin71654', ]
df.1 <- merge(df, lineage.size, by = 'assigned_lineage')
ggplot(df.1, aes(x = fatepotential_d0_d10, y = bias)) +
  geom_point(aes(size = d0.lin.size), shape = 21) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 0.5, linetype = 'dashed') +
  scale_size_continuous(breaks = c(2, 4, 6, 8, 10)) +
  labs(title = treatment,
       x = 'Fate potential from d0 to d10',
       y = 'Fate bias from d0 to w5-adapting') +
  theme_Publication()

high.bias <- df.1[df.1$bias > 0.5, ] %>% pull(assigned_lineage)
more.than.1.cell <- lineage.size[lineage.size$d0.lin.size > 1, ] %>% pull(assigned_lineage)
df.2 <- df.1[df.1$assigned_lineage %in% high.bias, ]
df.2$assigned_lineage_2 <- ifelse(df.2$assigned_lineage %in% more.than.1.cell, df.2$assigned_lineage, 'Singleton')

ggplot(df.2, aes(x = fatepotential_d0_d10, y = bias)) +
  geom_point(aes(color = assigned_lineage_2), shape = 21) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 0.5, linetype = 'dashed') +
  facet_wrap(~assigned_lineage_2) +
  # scale_size_continuous(breaks = c(2, 4, 6, 8, 10)) +
  labs(title = treatment,
       x = 'Fate potential from d0 to d10',
       y = 'Fate bias from d0 to w5-adapting') +
  theme_Publication()


# ===============================================================================================
# ===============================================================================================
# ===============================================================================================

# =============================================================================
# Read data
# =============================================================================
out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'

load(paste0(out_dir, 'TFchromVAR_on_day0_cor_vec_DABTRAM.RData'))
cor_vec.DABTRAM <- chromVAR_cor_vec

load('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/chromVAR_cor_vec.RData')

dabtram_d0_chromVAR_cor_vec <- chromVAR_cor_vec[['dabtram_d0_chromVAR_cor_vec']] 

# =============================================================================
# Wrangle
# =============================================================================
cor_vec.DABTRAM$TF <- rownames(cor_vec.DABTRAM)
colnames(cor_vec.DABTRAM) <- c('correlation.d0.fb', 'p.value.d0.fb', 'TF')

dabtram_d0_chromVAR_cor_vec$TF <- rownames(dabtram_d0_chromVAR_cor_vec)
colnames(dabtram_d0_chromVAR_cor_vec) <- c('correlation.fp.d0.d10', 'p.value.fp.d0.d10', 'TF')

comp.df <- merge(cor_vec.DABTRAM, dabtram_d0_chromVAR_cor_vec, by = 'TF')
tf.top.fb <- comp.df %>% arrange(desc(correlation.d0.fb)) %>% head(1)
tf.top.fp <- comp.df %>% arrange(desc(correlation.fp.d0.d10)) %>% head(1)

ggplot(comp.df, aes(x = correlation.fp.d0.d10, y = correlation.d0.fb)) +
  geom_point(color= 'gray') +
  geom_point(data = tf.top.fb, color = 'red', size = 3) +
  geom_point(data = tf.top.fp, color = 'red', size = 3) +
  ggrepel::geom_text_repel(data = tf.top.fb, aes(label = TF), color = 'red') +
  ggrepel::geom_text_repel(data = tf.top.fp, aes(label = TF), color = 'red') +
  labs(title = treatment,
       x = 'Cor with fate potential from d0 to d10',
       y = 'Cor with fate bias from d0 to w5-adapting') +
  theme_Publication()

ggsave(paste0(out_dir, 'comparison_TFchromVAR_fp_d0_d10_vs_fb_d0_', treatment, '.png'), width = 4, height = 4, dpi = 300)


# =============================================================================
rm(list = ls())
# =============================================================================
# Read data
# =============================================================================
out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'

load(paste0(out_dir, 'geneSaver_on_day0_cor_vec_DABTRAM.RData'))
cor_vec.DABTRAM <- cor_vec

load('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/saver_cor_vec.RData')

dabtram_d0_cor_vec <- saver_cor_vec[['dabtram_d0_saver_cor_vec']] 

# =============================================================================
# Wrangle
# =============================================================================
colnames(cor_vec.DABTRAM) <- c('gene', 'correlation.d0.fb', 'p.value.d0.fb')

dabtram_d0_cor_vec$gene <- rownames(dabtram_d0_cor_vec)
colnames(dabtram_d0_cor_vec) <- c('correlation.fp.d0.d10', 'p.value.fp.d0.d10', 'gene')

comp.df <- merge(cor_vec.DABTRAM, dabtram_d0_cor_vec, by = 'gene')
comp.df <- comp.df[!grepl('\\.', comp.df$gene),]

gene.top.fb <- comp.df %>% arrange(desc(correlation.d0.fb)) %>% head(10)
gene.top.fp <- comp.df %>% arrange(desc(correlation.fp.d0.d10)) %>% head(10)

ggplot(comp.df, aes(x = correlation.fp.d0.d10, y = correlation.d0.fb)) +
  geom_point(color= 'gray') +
  geom_point(data = gene.top.fb, color = 'red', size = 2) +
  geom_point(data = gene.top.fp, color = 'red', size = 2) +
  ggrepel::geom_text_repel(data = gene.top.fb, aes(label = gene), color = 'red',
                           min.segment.length = 0.1, max.overlaps = 20) +
  ggrepel::geom_text_repel(data = gene.top.fp, aes(label = gene), color = 'red',
                           min.segment.length = 0.1, max.overlaps = 20) +
  labs(title = treatment,
       x = 'Cor with fate potential from d0 to d10',
       y = 'Cor with fate bias from d0 to w5-adapting') +
  theme_Publication()

ggsave(paste0(out_dir, 'comparison_geneSaver_fp_d0_d10_vs_fb_d0_', treatment, '.png'), width = 4, height = 4, dpi = 300)
