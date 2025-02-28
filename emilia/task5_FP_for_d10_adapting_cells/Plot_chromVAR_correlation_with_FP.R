library(tidyverse)
library(ggplot2)
library(GGally)

keyTFs <- c('JUN', 'FOS::JUN', 'FOSB::JUN', 'FOSL1::JUN', 'FOSL2::JUN', 'JUN::JUNB', 
            'FOS::JUNB', 'FOSB::JUNB', 'FOSL1::JUNB', 'FOSL2::JUNB', 'FOS::JUND',
            'FOSL1::JUND', 'FOSL2::JUND', 'JUNB', 'JUND',
            'TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 
            'SOX10', 'MITF', 
            'SNAI1', 'SNAI2', 'SNAI3')

# =============================================================================
# Read results
# =============================================================================

dabtram_cor_vec <- read.csv('~/Downloads/cor_vec_TF_DABTRAM.csv')
cocl2_cor_vec <- read.csv('~/Downloads/cor_vec_TF_COCL2.csv')
cis_cor_vec <- read.csv('~/Downloads/cor_vec_TF_CIS.csv')

# =============================================================================
# Wrangle data
# =============================================================================
dabtram_cor_vec <- dabtram_cor_vec[, c('TF', 'cor')]
colnames(dabtram_cor_vec) <- c('TF', 'cor.DABTRAM')

cocl2_cor_vec <- cocl2_cor_vec[, c('TF', 'cor')]
colnames(cocl2_cor_vec) <- c('TF', 'cor.COCL2')

cis_cor_vec <- cis_cor_vec[, c('TF', 'cor')]
colnames(cis_cor_vec) <- c('TF', 'cor.CIS')

df <- merge(dabtram_cor_vec, cocl2_cor_vec, by = 'TF')
df <- merge(df, cis_cor_vec, by = 'TF')
df$keyTFs <- ifelse(df$TF %in% keyTFs, 'keyTFs', 'other')

# =============================================================================
# Plot
# =============================================================================
# df <- df[!grepl('\\.', df$TF), ]

# use ggpairs and highligh specific data points and put the keyTFs above others
df <- df %>% arrange(desc(keyTFs))
ggpairs(df, 
        columns = c(2, 3, 4), 
        aes(color = keyTFs),
        lower = list(continuous = 'points'),
        diag = list(continuous = "densityDiag", discrete = "barDiag", na = "naDiag")) +
  theme_bw()

ggpairs(df, 
        columns = c(2, 3, 4), 
        aes(color = keygene),
        lower = list(continuous = 'points'),
        diag = list(continuous = "densityDiag", discrete = "barDiag", na = "naDiag")) +
  theme_bw()

p1 <- ggplot(df, aes(x = cor.DABTRAM, y = cor.COCL2, color = keyTFs)) +
  geom_point() +
  stat_cor() +
  theme_bw() +
  theme(panel.background = element_rect(colour = NA),
        plot.background = element_rect(colour = NA),
        panel.grid.major = element_line(colour="white"),
        panel.grid.minor = element_blank(),
        legend.position = 'none') +
  labs(x = 'corr.DABTRAM', y = 'corr.COCL2')
p2 <- ggplot(df, aes(x = cor.DABTRAM, y = cor.CIS, color = keyTFs)) +
  geom_point() +
  stat_cor() +
  theme_bw() +
  theme(panel.background = element_rect(colour = NA),
        plot.background = element_rect(colour = NA),
        legend.position = 'none',
        panel.grid.major = element_line(colour="white"),
        panel.grid.minor = element_blank()) +
  labs(x = 'corr.DABTRAM', y = 'corr.CIS')
p3 <- ggplot(df, aes(x = cor.CIS, y = cor.COCL2, color = keyTFs)) +
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

ggsave('~/Downloads/cor_d0_w5_TF.pdf', p4, width = 9, height = 3)

# =============================================================================
# Find shared positively regulated
# =============================================================================

df.pos <- df %>% 
  filter(cor.DABTRAM > 0,  cor.CIS > 0)

# =============================================================================
# Compare COCL2 with original
# =============================================================================

out_dir2 <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/'
load(paste0(out_dir2, 'chromVAR_cor_vec.RData'))

treatment <- 'dabtram'
tm_d0_d10 <- chromVAR_cor_vec[[paste0(treatment, "_d0_chromVAR_cor_vec")]]
tm_d10_w5 <- chromVAR_cor_vec[[paste0(treatment, "_d10_chromVAR_cor_vec")]]

tm_d0_d10$TF <- rownames(tm_d0_d10)
tm_d10_w5$TF <- rownames(tm_d10_w5)

colnames(tm_d0_d10) <- c('cor.d0.d10', 'pval.d0.d10', 'TF')
colnames(tm_d10_w5) <- c('cor.d10.w5', 'pval.d10.w5', 'TF')

df2 <- merge(tm_d0_d10, tm_d10_w5, by = 'TF')
df2 <- merge(df2, dabtram_cor_vec, by = 'TF')

df2$keyTFs <- ifelse(df2$TF %in% keyTFs, 'keyTFs', 'other')

df2 <- df2[, c('TF', 'cor.d0.d10', 'cor.d10.w5', paste0('cor.', toupper(treatment)), 
               'keyTFs', 'pval.d0.d10', 'pval.d10.w5')]

# df2 <- df2[!grepl('\\.', df2$TF), ]
df2 <- df2 %>% arrange(desc(keyTFs))

p4 <- ggpairs(df2, 
        columns = c(2, 3, 4), 
        aes(color = keyTFs),
        lower = list(continuous = 'points')) +
  theme_bw() +
  theme(panel.background = element_rect(colour = NA),
        plot.background = element_rect(colour = NA),
        panel.grid.major = element_line(colour="white"),
        panel.grid.minor = element_blank())
ggsave('~/Downloads/cor_d0_w5_TF_dabtram.pdf', p4, width = 5, height = 5)

ggplot(df2) +
  geom_bar(aes(x = TF, y = cor.DABTRAM), stat = 'identity')

p4 <- ggplot(df2, aes(x = cor.d0.d10, y = cor.CIS)) +
  geom_point(aes(color = keyTFs)) +
  stat_cor() +
  theme_bw() +
  theme(panel.background = element_rect(colour = NA),
        plot.background = element_rect(colour = NA),
        panel.grid.major = element_line(colour="white"),
        panel.grid.minor = element_blank(),
        legend.position = 'none')
p4
