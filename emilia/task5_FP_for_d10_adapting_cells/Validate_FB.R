rm(list=ls())
library(Seurat)
library(multiomeFate)
library(tidyverse)
library(ggpubr)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'
figure_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/figures/MultiomeFate_V3/Fig6/'


remove_unassigned_cells <- TRUE


# =============================================================================
# Read data
# =============================================================================
load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))

all_data@misc <- all_data_fatepotential

# remove cells with no lineage
if(remove_unassigned_cells) {
  print("Removing cells with no assigned lineage")
  all_data$keep <- !is.na(all_data$assigned_lineage)
  if(any(!all_data$keep)){
    print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
    all_data <- subset(all_data, keep == TRUE)
  }
}

df.bias.DABTRAM <- read.csv(paste0(out_dir, 'adapting_bias_thres_0_DABTRAM.csv'))
df.bias.COCL2 <- read.csv(paste0(out_dir, 'adapting_bias_thres_0_COCL2.csv'))
df.bias.CIS <- read.csv(paste0(out_dir, 'adapting_bias_thres_0_CIS.csv'))

# =============================================================================
# Wrangle data
# =============================================================================

metadata <- all_data@meta.data
metadata$cell_id <- rownames(metadata)
metadata.day0 <- subset(metadata, dataset == 'day0')
metadata.day0 <- metadata.day0[, c('cell_id', 'assigned_lineage')]
metadata.day10.COCL2 <- subset(metadata, dataset == 'day10_COCL2')
metadata.week5.COCL2 <- subset(metadata, dataset == 'week5_COCL2')
metadata.day10.CIS <- subset(metadata, dataset == 'day10_CIS')
metadata.week5.CIS <- subset(metadata, dataset == 'week5_CIS')
metadata.day10.DABTRAM <- subset(metadata, dataset == 'day10_DABTRAM')
metadata.week5.DABTRAM <- subset(metadata, dataset == 'week5_DABTRAM')

fp.d0.d10.CIS <- all_data@misc[["fatepotential_CIS_d0_d10"]][["lineage_imputed_count"]]
fp.d0.d10.COCL2 <- all_data@misc[["fatepotential_COCL2_d0_d10"]][["lineage_imputed_count"]]
fp.d10.w5.COCL2 <- all_data@misc[["fatepotential_COCL2_d10_w5"]][["lineage_imputed_count"]]
fp.d10.w5.DABTRAM <- all_data@misc[["fatepotential_DABTRAM_d10_w5"]][["lineage_imputed_count"]]
fp.d10.w5.CIS <- all_data@misc[["fatepotential_CIS_d10_w5"]][["lineage_imputed_count"]]

fp.d0.d10.COCL2 <- as.data.frame(fp.d0.d10.COCL2)
fp.d0.d10.COCL2$assigned_lineage <- rownames(fp.d0.d10.COCL2)
fp.d0.d10.COCL2.fast.polif <- fp.d0.d10.COCL2 %>% 
  filter(fp.d0.d10.COCL2 > 3.398170) %>% 
  pull(assigned_lineage)

fp.d0.d10.CIS <- as.data.frame(fp.d0.d10.CIS)
fp.d0.d10.CIS$assigned_lineage <- rownames(fp.d0.d10.CIS)
fp.d0.d10.CIS.fast.polif <- fp.d0.d10.CIS %>% 
  filter(fp.d0.d10.CIS > 3.729425) %>% 
  pull(assigned_lineage)

fp.d10.w5.COCL2 <- as.data.frame(fp.d10.w5.COCL2)
fp.d10.w5.COCL2$assigned_lineage <- rownames(fp.d10.w5.COCL2)
fp.d10.w5.COCL2.winner <- fp.d10.w5.COCL2 %>% 
  filter(fp.d10.w5.COCL2 > 1)

fp.d10.w5.CIS <- as.data.frame(fp.d10.w5.CIS)
fp.d10.w5.CIS$assigned_lineage <- rownames(fp.d10.w5.CIS)
fp.d10.w5.CIS.winner <- fp.d10.w5.CIS %>% 
  filter(fp.d10.w5.CIS > 1)

fp.d10.w5.DABTRAM <- as.data.frame(fp.d10.w5.DABTRAM)
fp.d10.w5.DABTRAM$assigned_lineage <- rownames(fp.d10.w5.DABTRAM)
fp.d10.w5.DABTRAM.winner <- fp.d10.w5.DABTRAM %>% 
  filter(fp.d10.w5.DABTRAM > 1)

df.bias.COCL2 <- merge(df.bias.COCL2, metadata.day0[, c('cell_id', 'assigned_lineage')], by = 'cell_id')
fp.bias.high.COCL2 <- df.bias.COCL2[df.bias.COCL2$bias > 0.031872810, ] %>% pull(assigned_lineage) %>% unique()
fp.adaptingFP.high.COCL2 <- df.bias.COCL2[df.bias.COCL2$adaptingFP > -1.349308, ] %>% pull(assigned_lineage) %>% unique()

intersect(union(union(fp.bias.high.COCL2, fp.adaptingFP.high.COCL2), 
                fp.d0.d10.COCL2.fast.polif),
          metadata.week5.COCL2$assigned_lineage) %>% length()
length(union(union(fp.bias.high.COCL2, fp.adaptingFP.high.COCL2),
             fp.d0.d10.COCL2.fast.polif))


length(intersect(metadata.day0$assigned_lineage, metadata.week5.COCL2$assigned_lineage)) /
length(unique(metadata.day0$assigned_lineage))

intersect(union(fp.bias.high.COCL2, fp.adaptingFP.high.COCL2), fp.d10.w5.COCL2.winner$assigned_lineage) %>% length()
length(union(fp.bias.high.COCL2, fp.adaptingFP.high.COCL2))

intersect(union(union(fp.bias.high.COCL2, fp.adaptingFP.high.COCL2), 
                fp.d0.d10.COCL2.fast.polif),
          fp.d10.w5.COCL2.winner$assigned_lineage) %>% length() /
length(union(union(fp.bias.high.COCL2, fp.adaptingFP.high.COCL2),
             fp.d0.d10.COCL2.fast.polif))

length(intersect(metadata.day0$assigned_lineage, fp.d10.w5.COCL2.winner$assigned_lineage)) /
length(unique(metadata.day0$assigned_lineage))




df.bias.CIS <- merge(df.bias.CIS, metadata.day0[, c('cell_id', 'assigned_lineage')], by = 'cell_id')
fp.bias.high.CIS <- df.bias.CIS[df.bias.CIS$bias > 0.20749704, ] %>% pull(assigned_lineage) %>% unique()
fp.adaptingFP.high.CIS <- df.bias.CIS[df.bias.CIS$adaptingFP > -0.6522971, ] %>% pull(assigned_lineage) %>% unique()

intersect(union(fp.bias.high.CIS, fp.adaptingFP.high.CIS), metadata.week5.CIS$assigned_lineage) %>% length()
length(union(fp.bias.high.CIS, fp.adaptingFP.high.CIS))

intersect(union(union(fp.bias.high.CIS, fp.adaptingFP.high.CIS), 
                fp.d0.d10.CIS.fast.polif),
          metadata.week5.CIS$assigned_lineage) %>% length() /
length(union(union(fp.bias.high.CIS, fp.adaptingFP.high.CIS),
             fp.d0.d10.CIS.fast.polif))

length(intersect(metadata.day0$assigned_lineage, metadata.week5.CIS$assigned_lineage)) /
length(unique(metadata.day0$assigned_lineage))




intersect(union(fp.bias.high.CIS, fp.adaptingFP.high.CIS), fp.d10.w5.CIS.winner$assigned_lineage) %>% length()
length(union(fp.bias.high.CIS, fp.adaptingFP.high.CIS))

intersect(union(union(fp.bias.high.CIS, fp.adaptingFP.high.CIS), 
                fp.d0.d10.CIS.fast.polif),
          fp.d10.w5.CIS.winner$assigned_lineage) %>% length()
length(union(union(fp.bias.high.CIS, fp.adaptingFP.high.CIS),
             fp.d0.d10.CIS.fast.polif))


length(intersect(metadata.day0$assigned_lineage, fp.d10.w5.CIS.winner$assigned_lineage)) /
length(unique(metadata.day0$assigned_lineage))






fp.d0.d10.COCL2.c <- all_data@misc[["fatepotential_COCL2_d0_d10"]][["cell_imputed_score"]]
fp.d0.d10.COCL2.c <- as.data.frame(fp.d0.d10.COCL2.c)
fp.d0.d10.COCL2.c$cell_id <- rownames(fp.d0.d10.COCL2.c)
df.COCL2 <- merge(fp.d0.d10.COCL2.c, df.bias.COCL2, by = 'cell_id') 
df.COCL2$is.lin.in.w5 <- df.COCL2$assigned_lineage %in% metadata.week5.COCL2$assigned_lineage
df.COCL2$is.lin.in.w5.imputed <- df.COCL2$assigned_lineage %in% fp.d10.w5.COCL2.winner$assigned_lineage

ggplot(df.COCL2, aes(x = fp.d0.d10.COCL2.c, y = bias)) +
  geom_point(aes(color = adaptingFP)) +
  facet_wrap(~is.lin.in.w5) +
  geom_smooth(method = 'lm') +
  theme_minimal()


fp.d0.d10.CIS.c <- all_data@misc[["fatepotential_CIS_d0_d10"]][["cell_imputed_score"]]
fp.d0.d10.CIS.c <- as.data.frame(fp.d0.d10.CIS.c)
fp.d0.d10.CIS.c$cell_id <- rownames(fp.d0.d10.CIS.c)
df.CIS <- merge(fp.d0.d10.CIS.c, df.bias.CIS, by = 'cell_id') 
df.CIS$is.lin.in.d10 <- df.CIS$assigned_lineage %in% metadata.day10.CIS$assigned_lineage
df.CIS$is.lin.in.w5 <- df.CIS$assigned_lineage %in% metadata.week5.CIS$assigned_lineage
df.CIS$is.lin.in.w5.imputed <- df.CIS$assigned_lineage %in% fp.d10.w5.CIS.winner$assigned_lineage

ggplot(df.CIS, aes(x = fp.d0.d10.CIS.c, y = bias)) +
  geom_point(aes(color = adaptingFP)) +
  facet_wrap(~is.lin.in.d10) +
  geom_smooth(method = 'lm') +
  theme_minimal()





df.COCL2$adaptingCount <- 10**df.COCL2$adaptingFP
df.COCL2.summary <- df.COCL2 %>% 
  group_by(assigned_lineage, is.lin.in.w5, is.lin.in.w5.imputed) %>% 
  summarize(adaptingCount = sum(adaptingCount))

ggplot(df.COCL2.summary, aes(x = is.lin.in.w5.imputed, y = adaptingCount)) +
  geom_violin(scale = 'width') +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  stat_compare_means() +
  theme_minimal() 


df.CIS$adaptingCount <- 10**df.CIS$adaptingFP
df.CIS.summary <- df.CIS %>% 
  group_by(assigned_lineage, is.lin.in.w5, is.lin.in.w5.imputed) %>% 
  summarize(adaptingCount = sum(adaptingCount))

ggplot(df.CIS.summary, aes(x = is.lin.in.w5.imputed, y = adaptingCount)) +
  geom_violin(scale = 'width') +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  stat_compare_means() +
  theme_minimal() 

df.bias.DABTRAM <- merge(df.bias.DABTRAM, metadata.day0[, c('cell_id', 'assigned_lineage')], by = 'cell_id')
df.bias.DABTRAM$is.lin.in.w5 <- df.bias.DABTRAM$assigned_lineage %in% metadata.week5.DABTRAM$assigned_lineage
df.bias.DABTRAM$is.lin.in.w5.imputed <- df.bias.DABTRAM$assigned_lineage %in% fp.d10.w5.DABTRAM.winner$assigned_lineage
df.bias.DABTRAM$adaptingCount <- 10**df.bias.DABTRAM$adaptingFP
df.DABTRAM.summary <- df.bias.DABTRAM %>% 
  group_by(assigned_lineage, is.lin.in.w5, is.lin.in.w5.imputed) %>% 
  summarize(adaptingCount = sum(adaptingCount))

ggplot(df.DABTRAM.summary, aes(x = is.lin.in.w5.imputed, y = adaptingCount)) +
  geom_violin(scale = 'width') +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  stat_compare_means() +
  theme_minimal() 

df.DABTRAM.summary2 <- df.bias.DABTRAM %>% 
  group_by(assigned_lineage, is.lin.in.w5, is.lin.in.w5.imputed) %>% 
  summarize(bias = sum(bias))

ggplot(df.DABTRAM.summary2, aes(x = is.lin.in.w5, y = bias)) +
  geom_violin(scale = 'width') +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  stat_compare_means() +
  theme_minimal() 



