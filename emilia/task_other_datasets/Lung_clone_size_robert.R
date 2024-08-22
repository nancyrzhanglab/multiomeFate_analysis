library(data.table)
library(ggplot2)
library(GGally)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

in_dir <- '/Users/emiliac/Dropbox/MultiomeFate/data/ShafferLab/robert_2024_05_03_HCC827_gDNA_Sequencing_SMS008/R_Scripts_and_Plots/Organized_Data/'
out_dir <- '~/Dropbox/Thesis/Lineage_trace/outputs/task9_lung_clone_size_robert/'


# ==============================================================================
# Read data
# ==============================================================================
data.without.spike <- read.csv(paste0(in_dir, 'Data_Without_Spike_Ins.csv'), row.names = 1)
data.spike_in <- read.csv(paste0(in_dir, 'Spike_Ins_With_Data.csv'), row.names = 1)
meta.data <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/Lung_robert_metadata.csv')

# ==============================================================================
# Wrangle data
# ==============================================================================
meta.data$Sample <- paste0('Sample_', meta.data$Sample)
shared_samples <- intersect(meta.data$Sample, colnames(data.without.spike))

rownames(meta.data) <- meta.data$Sample
meta.data$Sample_Name <- paste0(meta.data$Sample, '_', meta.data$Treatment, '_', meta.data$Time)

meta.data <- meta.data[colnames(data.without.spike), ]
colnames(data.without.spike) <- meta.data$Sample_Name

meta.data1 <- meta.data[colnames(data.spike_in)[5:16], ]
colnames(data.spike_in) <- c(colnames(data.spike_in)[1:4], meta.data1$Sample_Name) 


# ==============================================================================
# Plotting (scatter plots)
# ==============================================================================
osi <- data.without.spike[, c('Sample_11_Osi_day10', 'Sample_12_Osi_day10', 'Sample_40_Osi_END', 'Sample_41_Osi_END')]
osi$all <- rowSums(osi)
osi <- osi[osi$all > 0, ]
osi <- log10(osi+1)
osi <- subset(osi, select = -c(all))

# osi$count <- rowSums(osi!=0)
# osi <- osi[osi$count > 1, ]
# osi <- subset(osi, select = -c(count))

osi.spike <- data.spike_in[, c('Sample_11_Osi_day10', 'Sample_12_Osi_day10', 'Sample_40_Osi_END', 'Sample_41_Osi_END')]
osi.spike <- log10(osi.spike + 1)

osi$is.spike <- 'non-spike'
osi.spike$is.spike <- 'spike'

osi <- rbind(osi, osi.spike)

ggpairs(osi, 
        mapping = aes(alpha = 0.1),
        upper = list(continuous="density"),
        lower = list(continuous=wrap("points", position=position_jitter(height=0.5, width=0.5)))) +
  theme_bw()

ggpairs(osi, 
        mapping = aes(alpha = 0.1),
        lower = list(continuous=wrap("points", position=position_jitter(height=0.5, width=0.5)))) +
  theme_bw()

ggpairs(osi, columns = 1:4,
        ggplot2::aes(color = is.spike, alpha = 0.3),
        lower = list(continuous=wrap("points", position=position_jitter(height=0.2, width=0.2)))) +
  scale_color_manual(values = c('gray', 'red')) +
  theme_bw()
# ggsave(paste0(out_dir, 'Osi_clone_size_scatter.png'), width = 8, height = 8)


osi <- data.without.spike[, c('Sample_17_DMSO_day10', 'Sample_18_DMSO_day10', 'Sample_11_Osi_day10', 'Sample_12_Osi_day10')]
osi$all <- rowSums(osi)
osi <- osi[osi$all > 0, ]
osi <- log10(osi+1)
osi <- subset(osi, select = -c(all))

osi.spike <- data.spike_in[, c('Sample_17_DMSO_day10', 'Sample_18_DMSO_day10', 'Sample_11_Osi_day10', 'Sample_12_Osi_day10')]
osi.spike <- log10(osi.spike + 1)

osi$is.spike <- 'non-spike'
osi.spike$is.spike <- 'spike'

osi <- rbind(osi, osi.spike)
ggpairs(osi, columns = 1:4,
        ggplot2::aes(color = is.spike, alpha = 0.3),
        lower = list(continuous=wrap("points", position=position_jitter(height=0.2, width=0.2)))) +
  scale_color_manual(values = c('gray', 'red')) +
  theme_bw()
# ggsave(paste0(out_dir, 'Osi_log10_reads_scatter_dmso.png'), width = 8, height = 8)





data.spike_in.melt <- melt(data.spike_in, id.vars= c('Sequence', 'Num_Of_Cells', 'Trimmed_with_StartSeq', 'closest_match_in_barcodes'))
data.spike_in.melt$reads_log10 <- log10(data.spike_in.melt$value + 1)
data.spike_in.melt$Num_Of_Cells_log10 <- log10(data.spike_in.melt$Num_Of_Cells + 1)
ggplot(data.spike_in.melt, aes(x = Num_Of_Cells_log10, y = reads_log10)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  facet_wrap(. ~ variable) +
  theme_bw()

cis <- data.without.spike[, c('Sample_15_Cisplatin_day10', 'Sample_38_Cisplatin_END')]
cis$all <- rowSums(cis)
cis <- cis[cis$all > 0, ]
cis <- log10(cis+1)
cis <- subset(cis, select = -c(all))

cis.spike <- data.spike_in[, c('Sample_15_Cisplatin_day10', 'Sample_38_Cisplatin_END')]
cis.spike <- log10(cis.spike + 1)

cis$is.spike <- 'non-spike'
cis.spike$is.spike <- 'spike'

cis <- rbind(cis, cis.spike)
ggpairs(cis, 
        mapping = aes(alpha = 0.1)) +
  theme_bw()

ggpairs(cis, columns = 1:2,
        ggplot2::aes(color = is.spike, alpha = 0.3),
        lower = list(continuous=wrap("points", position=position_jitter(height=0.2, width=0.2)))) +
  scale_color_manual(values = c('gray', 'red')) +
  theme_bw()
# ggsave(paste0(out_dir, 'Cis_clone_size_scatter.png'), width = 5, height = 5)

cis <- data.without.spike[, c('Sample_17_DMSO_day10', 'Sample_18_DMSO_day10', 'Sample_15_Cisplatin_day10')]
cis$all <- rowSums(cis)
cis <- cis[cis$all > 0, ]
cis <- log10(cis+1)
cis <- subset(cis, select = -c(all))

cis.spike <- data.spike_in[, c('Sample_17_DMSO_day10', 'Sample_18_DMSO_day10', 'Sample_15_Cisplatin_day10')]
cis.spike <- log10(cis.spike + 1)

cis$is.spike <- 'non-spike'
cis.spike$is.spike <- 'spike'

cis <- rbind(cis, cis.spike)
ggpairs(cis, columns = 1:3,
        ggplot2::aes(color = is.spike, alpha = 0.3),
        lower = list(continuous=wrap("points", position=position_jitter(height=0.2, width=0.2)))) +
  scale_color_manual(values = c('gray', 'red')) +
  theme_bw()
ggsave(paste0(out_dir, 'Cis_clone_size_scatter_dmso.png'), width = 5, height = 5)


cocl2 <- data.without.spike[, c('Sample_5_CoCl2_day10', 'Sample_6_CoCl2_day10', 'Sample_43_CoCl2_END')]
cocl2$all <- rowSums(cocl2)
cocl2 <- cocl2[cocl2$all > 0, ]
cocl2 <- log10(cocl2+1)
cocl2 <- subset(cocl2, select = -c(all))

cocl2.spike <- data.spike_in[, c('Sample_5_CoCl2_day10', 'Sample_6_CoCl2_day10', 'Sample_43_CoCl2_END')]
cocl2.spike <- log10(cocl2.spike + 1)

cocl2$is.spike <- 'non-spike'
cocl2.spike$is.spike <- 'spike'
cocl2 <- rbind(cocl2, cocl2.spike)

ggpairs(cocl2, 
        mapping = aes(alpha = 0.1),
        lower = list(continuous=wrap("points", position=position_jitter(height=0.5, width=0.5)))) +
  theme_bw()
ggpairs(cocl2, columns = 1:3,
        ggplot2::aes(color = is.spike, alpha = 0.3),
        lower = list(continuous=wrap("points", position=position_jitter(height=0.2, width=0.2)))) +
  scale_color_manual(values = c('gray', 'red')) +
  theme_bw()
ggsave(paste0(out_dir, 'Cocl2_clone_size_scatter.png'), width = 8, height = 8)


cocl2 <- data.without.spike[, c('Sample_17_DMSO_day10', 'Sample_18_DMSO_day10', 'Sample_5_CoCl2_day10', 'Sample_6_CoCl2_day10')]
cocl2$all <- rowSums(cocl2)
cocl2 <- cocl2[cocl2$all > 0, ]
cocl2 <- log10(cocl2+1)
cocl2 <- subset(cocl2, select = -c(all))

cocl2.spike <- data.spike_in[, c('Sample_17_DMSO_day10', 'Sample_18_DMSO_day10', 'Sample_5_CoCl2_day10', 'Sample_6_CoCl2_day10')]
cocl2.spike <- log10(cocl2.spike + 1)

cocl2$is.spike <- 'non-spike'
cocl2.spike$is.spike <- 'spike'
cocl2 <- rbind(cocl2, cocl2.spike)

ggpairs(cocl2, columns = 1:4,
        ggplot2::aes(color = is.spike, alpha = 0.3),
        lower = list(continuous=wrap("points", position=position_jitter(height=0.2, width=0.2)))) +
  scale_color_manual(values = c('gray', 'red')) +
  theme_bw()
ggsave(paste0(out_dir, 'Cocl2_clone_size_scatter_dmso.png'), width = 8, height = 8)



day10 <- data.without.spike[, c('Sample_11_Osi_day10', 'Sample_12_Osi_day10', 'Sample_15_Cisplatin_day10',
                                'Sample_16_Cisplatin_day10', 'Sample_5_CoCl2_day10', 'Sample_6_CoCl2_day10')]
day10$all <- rowSums(day10)
day10 <- day10[day10$all > 0, ]
day10 <- log10(day10+1)
day10 <- subset(day10, select = -c(all))

day10.spike <- data.spike_in[, c('Sample_11_Osi_day10', 'Sample_12_Osi_day10', 'Sample_15_Cisplatin_day10',
                                 'Sample_16_Cisplatin_day10', 'Sample_5_CoCl2_day10', 'Sample_6_CoCl2_day10')]
day10.spike <- log10(day10.spike + 1)

day10$is.spike <- 'non-spike'
day10.spike$is.spike <- 'spike'

day10 <- rbind(day10, day10.spike)

ggpairs(day10, 
        mapping = aes(alpha = 0.1),
        lower = list(continuous=wrap("points", position=position_jitter(height=0.5, width=0.5)))) +
  theme_bw()
ggpairs(day10, columns = 1:6,
        ggplot2::aes(color = is.spike, alpha = 0.3),
        lower = list(continuous=wrap("points"))) +
  scale_color_manual(values = c('gray', 'red')) +
  theme_bw()
ggsave(paste0(out_dir, 'Day10_clone_size_scatter1.png'), width = 12, height = 12)


end <- data.without.spike[, c('Sample_40_Osi_END', 'Sample_41_Osi_END', 'Sample_43_CoCl2_END',
                              'Sample_38_Cisplatin_END')]
end$all <- rowSums(end)
end <- end[end$all > 0, ]
end <- log10(end+1)
end <- subset(end, select = -c(all))

end.spike <- data.spike_in[, c('Sample_40_Osi_END', 'Sample_41_Osi_END', 'Sample_43_CoCl2_END',
                               'Sample_38_Cisplatin_END')]
end.spike <- log10(end.spike + 1)

end$is.spike <- 'non-spike'
end.spike$is.spike <- 'spike'

end <- rbind(end, end.spike)


ggpairs(end, 
        mapping = aes(alpha = 0.1),
        lower = list(continuous=wrap("points", position=position_jitter(height=0.5, width=0.5)))) +
  theme_bw()

ggpairs(end, columns = 1:4,
        ggplot2::aes(color = is.spike, alpha = 0.3),
        lower = list(continuous=wrap("points", position=position_jitter(height=0.2, width=0.2)))) +
  scale_color_manual(values = c('gray', 'red')) +
  theme_bw()
ggsave(paste0(out_dir, 'End_clone_size_scatter.png'), width = 12, height = 12)

ggpairs(end, columns = 1:4,
        ggplot2::aes(color = is.spike, alpha = 0.3),
        lower = list(continuous=wrap("points"))) +
  scale_color_manual(values = c('gray', 'red')) +
  theme_bw()
ggsave(paste0(out_dir, 'End_clone_size_scatter_nojitter.png'), width = 12, height = 12)

# ==============================================================================
# Plotting (per sample)
# ==============================================================================
overlap <- intersect(colnames(data.spike_in), colnames(data.without.spike))
data.spike_in.melt <- melt(data.spike_in[, c('Num_Of_Cells','closest_match_in_barcodes',  overlap)], id.vars = c('Num_Of_Cells', 'closest_match_in_barcodes'))
data.spike_in.melt$category <- 'Spike'


osi.day10.spike <- data.spike_in.melt[data.spike_in.melt$variable %in% c('Sample_11_Osi_day10', 'Sample_12_Osi_day10'), ]
osi.day10.spike$value_log10 <- log10(osi.day10.spike$value)
ggplot(osi.day10) +
  geom_point(aes(x = Num_Of_Cells, y = value, color = variable)) +
  theme_bw()

osi.day10 <- data.without.spike[, c('Sample_11_Osi_day10', 'Sample_12_Osi_day10')]
osi.day10 <- melt(osi.day10)
osi.day10 <- osi.day10[osi.day10$value > 0, ]
osi.day10$value_log10 <- log10(osi.day10$value)


ggplot(osi.day10) +
  geom_histogram(aes(x = value_log10), bins = 100) +
  theme_bw()


data.without.spike$closest_match_in_barcodes <- rownames(data.without.spike)
data.without.spike.melt <- melt(data.without.spike, id.vars = 'closest_match_in_barcodes')
data.without.spike.melt$category <- 'Non-Spike'
data.without.spike.melt <- data.without.spike.melt[data.without.spike.melt$variable %in% overlap, ]
data.without.spike.melt$Num_Of_Cells <- NA
data.melt <- rbind(data.without.spike.melt, data.spike_in.melt)
write.csv(data.melt, '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/Lung_robert_data.csv', row.names = FALSE)

data.without.spike.melt$value_log10 <- log10(data.without.spike.melt$value + 1)
ggplot(data.without.spike.melt) +
  geom_violin(aes(x = variable, y =value_log10)) +
  geom_jitter(aes(x = variable, y =value_log10), alpha = 0.1) +
  geom_jitter(data = data.spike_in.melt, aes(x = variable, y = value_log10, color = Num_Of_Cells)) +
  scale_color_gradient(low=  'blue', high = 'red') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

data.spike_in.melt$value_log10 <- log10(data.spike_in.melt$value + 1)
data.spike_in.melt$Num_Of_Cells_log10 <- log10(data.spike_in.melt$Num_Of_Cells + 1)
ggplot(data.spike_in.melt, aes(x = Num_Of_Cells, y =value)) +
  geom_smooth(method = 'lm', alpha = 0.3) +
  geom_point() +
  facet_wrap(.~variable, scale = 'free_y') +
  ylab('# seq reads from spike-in cells') +
  xlab('# spike-in cells') +
  theme_bw()
ggsave(paste0(out_dir, 'Spike_in_scatter.png'), width = 9, height = 6)
