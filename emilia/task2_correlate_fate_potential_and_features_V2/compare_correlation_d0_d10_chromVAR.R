library(tidyverse)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

# ==============================================================================
# Read data
# ==============================================================================
load('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/chromVAR_cor_vec.RData')

dabtram_d10_chromVAR_cor_vec <- chromVAR_cor_vec[['dabtram_d10_chromVAR_cor_vec']] 
cocl2_d10_chromVAR_cor_vec <- chromVAR_cor_vec[['cocl2_d10_chromVAR_cor_vec']] 
cis_d10_chromVAR_cor_vec <- chromVAR_cor_vec[['cis_d10_chromVAR_cor_vec']] 

dabtram_d0_chromVAR_cor_vec <- chromVAR_cor_vec[['dabtram_d0_chromVAR_cor_vec']] 
cocl2_d0_chromVAR_cor_vec <- chromVAR_cor_vec[['cocl2_d0_chromVAR_cor_vec']] 
cis_d0_chromVAR_cor_vec <- chromVAR_cor_vec[['cis_d0_chromVAR_cor_vec']] 

# ==============================================================================
# Wrangle data
# ==============================================================================
dabtram_d10_chromVAR_cor_vec <- dabtram_d10_chromVAR_cor_vec %>% 
  rename(correlation.day10 = correlation,
         p_value.day10 = p.value)

cocl2_d10_chromVAR_cor_vec <- cocl2_d10_chromVAR_cor_vec %>% 
  rename(correlation.day10 = correlation,
         p_value.day10 = p.value)

cis_d10_chromVAR_cor_vec <- cis_d10_chromVAR_cor_vec %>% 
  rename(correlation.day10 = correlation,
         p_value.day10 = p.value)

dabtram_d0_chromVAR_cor_vec <- dabtram_d0_chromVAR_cor_vec %>% 
  rename(correlation.day0 = correlation,
         p_value.day0 = p.value)

cocl2_d0_chromVAR_cor_vec <- cocl2_d0_chromVAR_cor_vec %>% 
  rename(correlation.day0 = correlation,
         p_value.day0 = p.value)

cis_d0_chromVAR_cor_vec <- cis_d0_chromVAR_cor_vec %>% 
  rename(correlation.day0 = correlation,
         p_value.day0 = p.value)

dabtram_comp <- merge(dabtram_d10_chromVAR_cor_vec, dabtram_d0_chromVAR_cor_vec, by = 'row.names')
cocl2_comp <- merge(cocl2_d10_chromVAR_cor_vec, cocl2_d0_chromVAR_cor_vec, by = 'row.names')
cis_comp <- merge(cis_d10_chromVAR_cor_vec, cis_d0_chromVAR_cor_vec, by = 'row.names')

# ==============================================================================
# Plot
# ==============================================================================

dabtram_comp$correlation.day0 <- as.numeric(dabtram_comp$correlation.day0)
dabtram_comp$correlation.day10 <- as.numeric(dabtram_comp$correlation.day10)
p1 <- ggplot(dabtram_comp, aes(x = correlation.day0, y = correlation.day10)) +
  geom_point(alpha = 0.2) +
  stat_cor() +
  xlim(-1, 1) +
  ylim(-1, 1) +
  labs(title = 'DABTRAM (chromVAR)',
       x = 'Day 0',
       y = 'Day 10') +
  theme_bw()

p2 <- ggplot(cocl2_comp, aes(x = correlation.day0, y = correlation.day10)) +
  geom_point(alpha = 0.2) +
  stat_cor() +
  xlim(-0.5, 0.5) +
  ylim(-0.5, 0.5) +
  labs(title = 'COCL2 (chromVAR)',
       x = 'Day 0',
       y = 'Day 10') +
  theme_bw()

p3 <- ggplot(cis_comp, aes(x = correlation.day0, y = correlation.day10)) +
  geom_point(alpha = 0.2) +
  stat_cor() +
  xlim(-0.5, 0.5) +
  ylim(-0.5, 0.5) +
  labs(title = 'CIS (chromVAR)',
       x = 'Day 0',
       y = 'Day 10') +
  theme_bw()

grid.arrange(p1, p2, p3, ncol = 3)

