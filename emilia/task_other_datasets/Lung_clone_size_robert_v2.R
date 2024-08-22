library(dplyr)
library(tidyr)

# ==============================================================================
# Read data (Dylan)
# ==============================================================================
metadata <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/Writeup6m_all-data_metadat.csv', row.names = 1)
metadata <- metadata[, c('dataset', 'assigned_lineage')]
metadata$cell_id <- rownames(metadata)

# ==============================================================================
# Calculate frequency (Dylan)
# ==============================================================================
sample.size <- metadata %>% 
  group_by(dataset) %>% 
  summarise(n.total = n())

lin.size <- metadata %>% 
  group_by(dataset, assigned_lineage) %>% 
  summarise(n.cells = n())

lin.size <- merge(lin.size, sample.size, by = 'dataset')
lin.size$freq <- lin.size$n.cells / lin.size$n.total

min.freq <- min(lin.size$freq) # 0.0001089562

# ==============================================================================
# Read data (Robert)
# ==============================================================================
data.lung <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/Lung_robert_data.csv')

sample.size.lung <- data.lung %>% 
  group_by(variable) %>% 
  summarise(n.reads.total = sum(value))

lin.size <- data.lung %>% 
  group_by(variable, closest_match_in_barcodes, category, Num_Of_Cells) %>% 
  summarise(n.reads = sum(value))

lin.size <- merge(lin.size, sample.size.lung, by = 'variable')

lin.size$freq <- lin.size$n.reads / lin.size$n.reads.total

lin.size.filt <- lin.size[lin.size$freq > 0.0001089562, ]
lin.size.filt$n.reads.log10 <- log10(lin.size.filt$n.reads)

lin.size.filt.spike <- lin.size.filt[lin.size.filt$category == 'Spike', ]
ggplot(lin.size.filt) +
  geom_violin(aes(x = variable, y =n.reads.log10)) +
  geom_jitter(aes(x = variable, y =n.reads.log10), alpha = 0.1) +
  geom_jitter(data = lin.size.filt.spike, aes(x = variable, y = n.reads.log10, color = Num_Of_Cells)) +
  scale_color_gradient(low=  'blue', high = 'red') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

lin.size.filt.wide <- lin.size.filt[, c('variable', 'closest_match_in_barcodes', 'n.reads')]
lin.size.filt.wide <- spread(lin.size.filt.wide, key = variable, value = n.reads)

osi <- lin.size.filt.wide[, c('Sample_11_Osi_day10', 'Sample_12_Osi_day10', 'Sample_40_Osi_END', 'Sample_41_Osi_END')]
osi[is.na(osi)] <- 0
osi$all <- rowSums(osi)
osi <- osi[osi$all > 0, ]
osi <- log10(osi + 1)
osi <- subset(osi, select = -c(all))

ggpairs(osi, 
        mapping = aes(alpha = 0.1),
        lower = list(continuous=wrap("points", position=position_jitter(height=0.5, width=0.5)))) +
  theme_bw()

cis <- lin.size.filt.wide[, c('Sample_15_Cisplatin_day10', 'Sample_38_Cisplatin_END')]
cis[is.na(cis)] <- 0
cis$all <- rowSums(cis)
cis <- cis[cis$all > 0, ]
cis <- log10(cis+1)
cis <- subset(cis, select = -c(all))
ggpairs(cis, 
        mapping = aes(alpha = 0.1)) +
  theme_bw()

cocl2 <- lin.size.filt.wide[, c('Sample_5_CoCl2_day10', 'Sample_6_CoCl2_day10', 'Sample_43_CoCl2_END')]
cocl2[is.na(cocl2)] <- 0
cocl2$all <- rowSums(cocl2)
cocl2 <- cocl2[cocl2$all > 0, ]
cocl2 <- log10(cocl2+1)
cocl2 <- subset(cocl2, select = -c(all))
ggpairs(cocl2, 
        mapping = aes(alpha = 0.1),
        lower = list(continuous=wrap("points", position=position_jitter(height=0.5, width=0.5)))) +
  theme_bw()

day10 <- lin.size.filt.wide[, c('Sample_11_Osi_day10', 'Sample_12_Osi_day10', 'Sample_15_Cisplatin_day10',
                                'Sample_16_Cisplatin_day10', 'Sample_5_CoCl2_day10', 'Sample_6_CoCl2_day10')]
day10[is.na(day10)] <- 0
day10$all <- rowSums(day10)
day10 <- day10[day10$all > 0, ]
day10 <- log10(day10+1)
day10 <- subset(day10, select = -c(all))
ggpairs(day10, 
        mapping = aes(alpha = 0.1),
        lower = list(continuous=wrap("points", position=position_jitter(height=0.5, width=0.5)))) +
  theme_bw()

end <- lin.size.filt.wide[, c('Sample_40_Osi_END', 'Sample_41_Osi_END', 'Sample_43_CoCl2_END',
                              'Sample_38_Cisplatin_END')]
end[is.na(end)] <- 0
end$all <- rowSums(end)
end <- end[end$all > 0, ]
end <- log10(end+1)
end <- subset(end, select = -c(all))
ggpairs(end, 
        mapping = aes(alpha = 0.1),
        lower = list(continuous=wrap("points", position=position_jitter(height=0.5, width=0.5)))) +
  theme_bw()
