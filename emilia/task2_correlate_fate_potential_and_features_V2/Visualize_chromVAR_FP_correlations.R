rm(list = ls())
library(tidyverse)
library(ggplot2)
library(GGally)
library(gridExtra)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

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
load('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/chromVAR_cor_vec.RData')

dabtram_d0_chromVAR_cor_vec <- chromVAR_cor_vec[['dabtram_d0_chromVAR_cor_vec']] 
cocl2_d0_chromVAR_cor_vec <- chromVAR_cor_vec[['cocl2_d0_chromVAR_cor_vec']] 
cis_d0_chromVAR_cor_vec <- chromVAR_cor_vec[['cis_d0_chromVAR_cor_vec']] 

dabtram_d10_chromVAR_cor_vec <- chromVAR_cor_vec[['dabtram_d10_chromVAR_cor_vec']] 
cocl2_d10_chromVAR_cor_vec <- chromVAR_cor_vec[['cocl2_d10_chromVAR_cor_vec']] 
cis_d10_chromVAR_cor_vec <- chromVAR_cor_vec[['cis_d10_chromVAR_cor_vec']] 

# ==============================================================================
# Wrangle data
# ==============================================================================
colnames(dabtram_d0_chromVAR_cor_vec) <- paste0(colnames(dabtram_d0_chromVAR_cor_vec), '.DABTRAM_d0')
colnames(cocl2_d0_chromVAR_cor_vec) <- paste0(colnames(cocl2_d0_chromVAR_cor_vec), '.COCL2_d0')
colnames(cis_d0_chromVAR_cor_vec) <- paste0(colnames(cis_d0_chromVAR_cor_vec), '.CIS_d0')

colnames(dabtram_d10_chromVAR_cor_vec) <- paste0(colnames(dabtram_d10_chromVAR_cor_vec), '.DABTRAM_d10')
colnames(cocl2_d10_chromVAR_cor_vec) <- paste0(colnames(cocl2_d10_chromVAR_cor_vec), '.COCL2_d10')
colnames(cis_d10_chromVAR_cor_vec) <- paste0(colnames(cis_d10_chromVAR_cor_vec), '.CIS_d10')

# Day0
d0_cor <- merge(dabtram_d0_chromVAR_cor_vec, cocl2_d0_chromVAR_cor_vec, by = 'row.names')
rownames(d0_cor) <- d0_cor$Row.names
d0_cor <- d0_cor |> dplyr::select(-`Row.names`)

d0_cor <- merge(d0_cor, cis_d0_chromVAR_cor_vec, by = 'row.names')
rownames(d0_cor) <- d0_cor$Row.names
d0_cor <- d0_cor |> dplyr::select(-`Row.names`)


# Day10
d10_cor <- merge(dabtram_d10_chromVAR_cor_vec, cocl2_d10_chromVAR_cor_vec, by = 'row.names')
rownames(d10_cor) <- d10_cor$Row.names
d10_cor <- d10_cor |> dplyr::select(-Row.names)

d10_cor <- merge(d10_cor, cis_d10_chromVAR_cor_vec, by = 'row.names')
rownames(d10_cor) <- d10_cor$Row.names
d10_cor <- d10_cor |> dplyr::select(-Row.names)

# ==============================================================================
# Plotting
# ==============================================================================

ggpairs(d0_cor[, c('correlation.DABTRAM_d0', 'correlation.COCL2_d0', 'correlation.CIS_d0')],
        mapping = aes(alpha = 0.7),
        lower = list(continuous = wrap("points", color = "#9DBDFF")),
        diag = list(continuous = wrap("densityDiag", fill = "#9DBDFF"))) +
  theme_Publication()

ggsave('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/correlation_ATAC_chromVAR_FP_day0.png', 
       width = 8, height = 8, dpi = 300)

ggpairs(d10_cor[, c('correlation.DABTRAM_d10', 'correlation.COCL2_d10', 'correlation.CIS_d10')],
        mapping = aes(alpha = 0.7),
        lower = list(continuous = wrap("points", color = "#9DBDFF")),
        diag = list(continuous = wrap("densityDiag", fill = "#9DBDFF"))) +
  theme_Publication()
ggsave('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/correlation_ATAC_chromVAR_FP_day10.png', 
       width = 8, height = 8, dpi = 300)

# ==============================================================================
# Plotting individual ones
# ==============================================================================


p1 <- ggplot(d0_cor, aes(x = correlation.DABTRAM_d0, y = correlation.COCL2_d0)) +
  # geom_point(color = 'pink') +
  geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "pearson", label.x = -0.5, label.y = 0.45) +
  xlim(c(-0.505, 0.505)) +
  ylim(c(-0.505, 0.505)) +
  scale_fill_manual(values = brewer.pal(5, "GnBu")) +
  theme_Publication()

p2 <- ggplot(d0_cor, aes(x = correlation.CIS_d0, y = correlation.COCL2_d0)) +
  # geom_point(color = 'pink') +
  geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "pearson", label.x = -0.5, label.y = 0.45) +
  xlim(c(-0.505, 0.505)) +
  ylim(c(-0.505, 0.505)) +
  scale_fill_manual(values = brewer.pal(5, "GnBu")) +
  theme_Publication()


p3 <- ggplot(d0_cor, aes(x = correlation.DABTRAM_d0, y = correlation.CIS_d0)) +
  # geom_point(color = 'pink') +
  geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "pearson", label.x = -0.5, label.y = 0.45) +
  xlim(c(-0.505, 0.505)) +
  ylim(c(-0.505, 0.505)) +
  scale_fill_manual(values = brewer.pal(5, "GnBu")) +
  theme_Publication()

p4 <- grid.arrange(p1, p2, p3, ncol = 3)
ggsave('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/correlation_ChromVAR_FP_day0_panel.pdf', 
       p4, 
       width = 8, height = 3, dpi = 300)

d10_cor <- d10_cor %>% drop_na()

p4 <- ggplot(d10_cor, aes(x = correlation.DABTRAM_d10, y = correlation.COCL2_d10)) +
  # geom_point(color = 'pink') +
  geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "pearson", label.x = -0.95, label.y = 0.8) +
  xlim(c(-1, 1)) +
  ylim(c(-1, 1)) +
  scale_fill_manual(values = brewer.pal(5, "GnBu")) +
  theme_Publication()

p5 <- ggplot(d10_cor, aes(x = correlation.CIS_d10, y = correlation.COCL2_d10)) +
  # geom_point(color = 'pink') +
  geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "pearson", label.x = -0.2, label.y = 0.8) +
  xlim(c(-0.2, 0.2)) +
  ylim(c(-1, 1)) +
  scale_fill_manual(values = brewer.pal(5, "GnBu")) +
  theme_Publication()


p6 <- ggplot(d10_cor, aes(x = correlation.DABTRAM_d10, y = correlation.CIS_d10)) +
  # geom_point(color = 'pink') +
  geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "pearson", label.x = -0.95, label.y = 0.16) +
  xlim(c(-1, 1)) +
  ylim(c(-0.2, 0.2)) +
  scale_fill_manual(values = brewer.pal(5, "GnBu")) +
  theme_Publication()

p7 <- grid.arrange(p4, p5, p6, ncol = 3)

ggsave('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/correlation_ChromVAR_FP_day10_panel.pdf', 
       p7, 
       width = 8, height = 3, dpi = 300)
# ==============================================================================
# Save session info
# ==============================================================================
session_info <- devtools::session_info()



# =============================================================================
# Plot
# =============================================================================
d0_cor$TF <- rownames(d0_cor)
d10_cor$TF <- rownames(d10_cor)
tfs <- d0_cor$TF
tfs_toplot <- tfs[grepl('JUN', tfs) | grepl('FOS', tfs) | grepl('SOX10', tfs) | grepl('MITF', tfs) | grepl('TEAD', tfs) | grepl('STAT', tfs) | grepl('IRF3', tfs) ]
tfs_toplot <- tfs_toplot[!grepl('(var.2)', tfs_toplot)]

keyTFs.dabtram.d0 <- d0_cor[d0_cor$TF %in% tfs_toplot, 'correlation.DABTRAM_d0']
keyTFs.cocl2.d0 <- d0_cor[d0_cor$TF %in% tfs_toplot, 'correlation.COCL2_d0']
keyTFs.cis.d0 <- d0_cor[d0_cor$TF %in% tfs_toplot, 'correlation.CIS_d0']

keyTFs.dabtram.d10 <- d10_cor[d10_cor$TF %in% tfs_toplot, 'correlation.DABTRAM_d10']
keyTFs.cocl2.d10 <- d10_cor[d10_cor$TF %in% tfs_toplot, 'correlation.COCL2_d10']
keyTFs.cis.d10 <- d10_cor[d10_cor$TF %in% tfs_toplot, 'correlation.CIS_d10']

x_min <- min(c(d0_cor$correlation.DABTRAM_d0, d0_cor$correlation.COCL2_d0, 
               d0_cor$correlation.CIS_d0, d10_cor$correlation.DABTRAM_d10, 
               d10_cor$correlation.COCL2_d10, d10_cor$correlation.CIS_d10), rm.na = T)
x_min <- -1
x_max <- max(c(d0_cor$correlation.DABTRAM_d0, d0_cor$correlation.COCL2_d0, 
               d0_cor$correlation.CIS_d0, d10_cor$correlation.DABTRAM_d10, 
               d10_cor$correlation.COCL2_d10, d10_cor$correlation.CIS_d10), rm.na = T)

p1 <- ggplot(d0_cor, aes(x = correlation.DABTRAM_d0)) +
  geom_density(fill = 'grey') +
  geom_vline(xintercept = keyTFs.dabtram.d0, color = 'red') +
  ggrepel::geom_text_repel(data = subset(d0_cor, TF %in% tfs_toplot), 
                           aes(label = TF, y = 3)) +
  theme_Publication() +
  xlim(x_min, x_max) +
  labs(ylab = 'Density', xlab = 'Correlation (DABTRAM)')

p2 <- ggplot(d0_cor, aes(x = correlation.COCL2_d0)) +
  geom_density(fill = 'grey') +
  geom_vline(xintercept = keyTFs.cocl2.d0, color = 'red') +
  ggrepel::geom_text_repel(data = subset(d0_cor, TF %in% tfs_toplot), 
                           aes(label = TF, y = 3)) +
  theme_Publication() +
  xlim(-1, x_max) +
  labs(ylab = 'Density', xlab = 'Correlation (COCL2)')

p3 <- ggplot(d0_cor, aes(x = correlation.CIS_d0)) +
  geom_density(fill = 'grey') +
  geom_vline(xintercept = keyTFs.cis.d0, color = 'red') +
  ggrepel::geom_text_repel(data = subset(d0_cor, TF %in% tfs_toplot), 
                           aes(label = TF, y = 3)) +
  theme_Publication() +
  xlim(-1, x_max) +
  labs(ylab = 'Density', xlab = 'Correlation (CIS)')

p4 <- ggplot(d10_cor, aes(x = correlation.DABTRAM_d10)) +
  geom_density(fill = 'grey') +
  geom_vline(xintercept = keyTFs.dabtram.d10, color = 'red') +
  ggrepel::geom_text_repel(data = subset(d10_cor, TF %in% tfs_toplot), 
                           aes(label = TF, y = 3)) +
  theme_Publication() +
  xlim(-1, x_max) +
  labs(ylab = 'Density', xlab = 'Correlation (DABTRAM)')

p5 <- ggplot(d10_cor, aes(x = correlation.COCL2_d10)) +
  geom_density(fill = 'grey') +
  geom_vline(xintercept = keyTFs.cocl2.d10, color = 'red') +
  ggrepel::geom_text_repel(data = subset(d10_cor, TF %in% tfs_toplot), 
                           aes(label = TF, y = 3)) +
  theme_Publication() +
  xlim(-1, x_max) +
  labs(ylab = 'Density', xlab = 'Correlation (COCL2)')

p6 <- ggplot(d10_cor, aes(x = correlation.CIS_d10)) +
  geom_density(fill = 'grey') +
  geom_vline(xintercept = keyTFs.cis.d10, color = 'red') +
  ggrepel::geom_text_repel(data = subset(d10_cor, TF %in% tfs_toplot), 
                           aes(label = TF, y = 3)) +
  theme_Publication() +
  xlim(-1, x_max) +
  labs(ylab = 'Density', xlab = 'Correlation (CIS)')

p7 <- grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 1)

ggsave('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/TFchromVAR_correlation_density.pdf', 
       plot = p7, width = 6, height = 8)


d10_cor <- d10_cor[order(d10_cor$correlation.DABTRAM_d10, decreasing = T),]
d10_cor$order <- 1:nrow(d10_cor)
tfs_toAnno <- c('BATF::JUN', 'FOS','FOSL1::JUN', 'FOSB::JUNB', 'JUN', 'IRF3', 'MITF', 'SOX10', 'STAT1','STAT3', 'STAT1::STAT2', 'TEAD1')
top <- d10_cor %>% 
  mutate(padjust = p.adjust(p.value.DABTRAM_d10, method = 'BH')) %>%
  filter(correlation.DABTRAM_d10 > 0) %>% 
  filter(padjust < 0.05) %>%
  arrange(correlation.DABTRAM_d10) %>%
  head(1) %>% 
  pull(correlation.DABTRAM_d10)

bottom <- d10_cor %>% 
  mutate(padjust = p.adjust(p.value.DABTRAM_d10, method = 'BH')) %>%
  filter(correlation.DABTRAM_d10 < 0) %>% 
  filter(padjust < 0.05) %>%
  arrange(correlation.DABTRAM_d10) %>%
  tail(1) %>% 
  pull(correlation.DABTRAM_d10)

ggplot(d10_cor, aes(x = order, y = correlation.DABTRAM_d10)) +
  geom_point() +
  geom_point(data = subset(d10_cor, TF %in% tfs_toAnno), color = 'red') +
  ggrepel::geom_text_repel(data = subset(d10_cor, TF %in% tfs_toAnno), color = 'red',
                           aes(label = TF), box.padding = unit(0.3, 'lines'),
                           point.padding = unit(1, 'lines'),
                           max.overlaps = 50) +
  geom_hline(yintercept = top, color = 'gray', linetype = 'dashed') +
  geom_hline(yintercept = bottom, color = 'gray', linetype = 'dashed') +
  ggtitle('TF (chromVAR, D10-to-Week5)') +
  theme_Publication()

d0_cor <- d0_cor[order(d0_cor$correlation.DABTRAM_d0, decreasing = T),]
d0_cor$order <- 1:nrow(d0_cor)
top <- d0_cor %>% 
  mutate(padjust = p.adjust(p.value.DABTRAM_d0, method = 'BH')) %>%
  filter(correlation.DABTRAM_d0 > 0) %>% 
  filter(padjust < 0.05) %>%
  arrange(correlation.DABTRAM_d0) %>%
  head(1) %>% 
  pull(correlation.DABTRAM_d0)

bottom <- d0_cor %>% 
  mutate(padjust = p.adjust(p.value.DABTRAM_d0, method = 'BH')) %>%
  filter(correlation.DABTRAM_d0 < 0) %>% 
  filter(padjust < 0.05) %>%
  arrange(correlation.DABTRAM_d0) %>%
  tail(1) %>% 
  pull(correlation.DABTRAM_d0)
ggplot(d0_cor, aes(x = order, y = correlation.DABTRAM_d0)) +
  geom_point() +
  geom_point(data = subset(d0_cor, TF %in% tfs_toAnno), color = 'red') +
  ggrepel::geom_text_repel(data = subset(d0_cor, TF %in% tfs_toAnno), color = 'red',
                           aes(label = TF), box.padding = unit(0.3, 'lines'),
                           point.padding = unit(1, 'lines'),
                           max.overlaps = 50) +
  geom_hline(yintercept = top, color = 'gray', linetype = 'dashed') +
  geom_hline(yintercept = bottom, color = 'gray', linetype = 'dashed') +
  ggtitle('TF (chromVAR, D0-to-D10)') +
  theme_Publication()

