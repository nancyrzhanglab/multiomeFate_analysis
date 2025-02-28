library(tidyverse)
library(ggplot2)
library(GGally)

keygenes <- list(
  jackpot = sort(c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E",
                   "C1S", "FRZB", "SERPINB2", "SERPINE1", "NGFR",
                   "SERPINE2", "NDRG1", "FEZF1", "EGR3", "VGF",
                   "WNT5A", "POSTN", "PDGFRB", "NRG1", "VEGFC", "FOSL1",
                   "RUNX2", "LOXL2", "JUN", "PDGFRC", "CD44", "ID3")),
  DABTRAM = sort(c("AXL", "EGFR", "NGFR", "IGFBP5", "ANXA1",
                   "IGFBP7", "JUNB", "BASP1", "IER2", "JUN",
                   "CXCL12", "ANXA2", "FOS", "MMP2", "GLRX",
                   "IL6ST", "PRNP", "FOSB", "CTSL", "SLC12A8",
                   "TFPI2", "MYL6", "IFITM3", "CAV1", "CD44")),
  COCL2 = sort(c("CD44", "FN1", "HPCAL1", "SLC16A3", "IGFBP5",
                 "COL6A2", "MPC2", "PLIN2", "HLA-A", "IGFBP7",
                 "CAV1"))
)
keygenes <- unlist(keygenes)

# =============================================================================
# Read results
# =============================================================================

dabtram_cor_vec <- read.csv('~/Downloads/cor_vec_DABTRAM.csv')
cocl2_cor_vec <- read.csv('~/Downloads/cor_vec_COCL2.csv')
cis_cor_vec <- read.csv('~/Downloads/cor_vec_CIS.csv')

# =============================================================================
# Wrangle data
# =============================================================================
dabtram_cor_vec <- dabtram_cor_vec[, c('gene', 'cor')]
colnames(dabtram_cor_vec) <- c('gene', 'cor.DABTRAM')

cocl2_cor_vec <- cocl2_cor_vec[, c('gene', 'cor')]
colnames(cocl2_cor_vec) <- c('gene', 'cor.COCL2')

cis_cor_vec <- cis_cor_vec[, c('gene', 'cor')]
colnames(cis_cor_vec) <- c('gene', 'cor.CIS')

df <- merge(dabtram_cor_vec, cocl2_cor_vec, by = 'gene')
df <- merge(df, cis_cor_vec, by = 'gene')
df$keygene <- ifelse(df$gene %in% keygenes, 'keygenes', 'other')

# =============================================================================
# Plot
# =============================================================================
# df <- df[!grepl('\\.', df$gene), ]

# use ggpairs and highligh specific data points and put the keygenes above others
df <- df %>% arrange(desc(keygene))
ggpairs(df, 
        columns = c(2, 3, 4), 
        aes(color = keygene),
        lower = list(continuous = 'points'),
        diag = list(continuous = "densityDiag", discrete = "barDiag", na = "naDiag")) +
  theme_bw()

p1 <- ggplot(df, aes(x = cor.DABTRAM, y = cor.COCL2, color = keygene)) +
  geom_point() +
  stat_cor() +
  theme_bw() +
  theme(panel.background = element_rect(colour = NA),
        plot.background = element_rect(colour = NA),
        panel.grid.major = element_line(colour="white"),
        panel.grid.minor = element_blank(),
        legend.position = 'none') +
  labs(x = 'corr.DABTRAM', y = 'corr.COCL2')
p2 <- ggplot(df, aes(x = cor.DABTRAM, y = cor.CIS, color = keygene)) +
  geom_point() +
  stat_cor() +
  theme_bw() +
  theme(panel.background = element_rect(colour = NA),
        plot.background = element_rect(colour = NA),
        legend.position = 'none',
        panel.grid.major = element_line(colour="white"),
        panel.grid.minor = element_blank()) +
  labs(x = 'corr.DABTRAM', y = 'corr.CIS')
p3 <- ggplot(df, aes(x = cor.CIS, y = cor.COCL2, color = keygene)) +
  geom_point() +
  stat_cor() +
  theme_bw() +
  theme(panel.background = element_rect(colour = NA),
        plot.background = element_rect(colour = NA),
        panel.grid.major = element_line(colour="white"),
        panel.grid.minor = element_blank(),
        legend.position = 'none') +
  labs(x = 'corr.CIS', y = 'corr.COCL2')
p4 <- grid.arrange(p1, p2, p3, ncol = 3, nrow = 1)


ggsave('~/Downloads/cor_d0_w5_RNA.pdf', p4, width = 9, height = 3)

# =============================================================================
# Find shared positively regulated
# =============================================================================

df.pos <- df %>% 
  filter(cor.DABTRAM > 0,  cor.CIS > 0)

# =============================================================================
# Compare COCL2 with original
# =============================================================================

out_dir2 <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/'
load(paste0(out_dir2, 'saver_cor_vec.RData'))

treatment <- 'cis'
tm_d0_d10 <- saver_cor_vec[[paste0(treatment, "_d0_saver_cor_vec")]]
tm_d10_w5 <- saver_cor_vec[[paste0(treatment, "_d10_saver_cor_vec")]]

tm_d0_d10$gene <- rownames(tm_d0_d10)
tm_d10_w5$gene <- rownames(tm_d10_w5)

colnames(tm_d0_d10) <- c('cor.d0.d10', 'pval.d0.d10', 'gene')
colnames(tm_d10_w5) <- c('cor.d10.w5', 'pval.d10.w5', 'gene')

df2 <- merge(tm_d0_d10, tm_d10_w5, by = 'gene')
df2 <- merge(df2, cis_cor_vec, by = 'gene')

df2$keygene <- ifelse(df2$gene %in% keygenes, 'keygenes', 'other')

df2 <- df2[, c('gene', 'cor.d0.d10', 'cor.d10.w5', paste0('cor.', toupper(treatment)), 
               'keygene', 'pval.d0.d10', 'pval.d10.w5')]

# df2 <- df2[!grepl('\\.', df2$gene), ]
df2 <- df2 %>% arrange(desc(keygene))

ggpairs(df2, 
        columns = c(2, 3, 4), 
        aes(color = keygene),
        lower = list(continuous = 'points')) +
  theme_bw()

