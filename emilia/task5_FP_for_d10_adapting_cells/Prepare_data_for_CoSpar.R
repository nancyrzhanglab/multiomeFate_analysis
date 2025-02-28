rm(list=ls())
library(Seurat)
library(tidyverse)

set.seed(123)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/day0_for_CoSpar/'

remove_unassigned_cells <- TRUE
adapting_thres <- 0

date_of_run <- Sys.time()
session_info <- devtools::session_info()
# =============================================================================
# Read data
# =============================================================================
load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))
load(paste0(data_dir, 'Writeup10a_data_saver.RData'))
load(paste0(data_dir, 'Writeup10a_data_fasttopic_DABTRAM.RData'))
load(paste0(data_dir, 'Writeup10a_data_fasttopic_COCL2.RData'))
load(paste0(data_dir, 'Writeup10a_data_fasttopic_CIS.RData'))

all_data@misc <- all_data_fatepotential
all_data[['fasttopic_DABTRAM']] <- all_data_fasttopic_DABTRAM
all_data[["ft.DABTRAM.umap"]] <- all_data_ft_DABTRAM_umap

all_data[['fasttopic_COCL2']] <- all_data_fasttopic_COCL2
all_data[["ft.COCL2.umap"]] <- all_data_ft_COCL2_umap

all_data[['fasttopic_CIS']] <- all_data_fasttopic_CIS
all_data[["ft.CIS.umap"]] <- all_data_ft_CIS_umap

all_data[['Saver']] <- all_data_saver
all_data[["Saver.pca"]] <- all_data_saver_pca
all_data[["Saver.umap"]] <- all_data_saver_umap

# remove cells with no lineage
if(remove_unassigned_cells) {
  print("Removing cells with no assigned lineage")
  all_data$keep <- !is.na(all_data$assigned_lineage)
  if(any(!all_data$keep)){
    print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
    all_data <- subset(all_data, keep == TRUE)
  }
}
# =============================================================================
# Classify day10 adapting vs non-adapting cells
# =============================================================================
metadat <- all_data@meta.data

metadat.day0 <- metadat[metadat$dataset == 'day0',]
metadat.day10_DABTRAM <- metadat[metadat$dataset == 'day10_DABTRAM', ]
metadat.day10_COCL2 <- metadat[metadat$dataset == 'day10_COCL2', ]
metadat.day10_CIS <- metadat[metadat$dataset == 'day10_CIS', ]
metadat.week5_DABTRAM <- metadat[metadat$dataset == 'week5_DABTRAM', ]
metadat.week5_COCL2 <- metadat[metadat$dataset == 'week5_COCL2', ]
metadat.week5_CIS <- metadat[metadat$dataset == 'week5_CIS', ]

nrow(metadat) == nrow(metadat.day0) +
                 nrow(metadat.day10_DABTRAM) + nrow(metadat.day10_COCL2) + nrow(metadat.day10_CIS) + 
                 nrow(metadat.week5_DABTRAM) + nrow(metadat.week5_COCL2) + nrow(metadat.week5_CIS)

metadat.day0 <- metadat.day0 %>% 
  mutate(day10_adapting.DABTRAM = 'day0',
         day10_adapting.COCL2 = 'day0',
         day10_adapting.CIS = 'day0')

metadat.day10_DABTRAM <- metadat.day10_DABTRAM %>% 
  mutate(day10_adapting.DABTRAM = ifelse(fatepotential_DABTRAM_d10_w5 > adapting_thres, 'adapting', 'non-adapting'),
         day10_adapting.COCL2 = NA,
         day10_adapting.CIS = NA)

metadat.day10_COCL2 <- metadat.day10_COCL2 %>%
  mutate(day10_adapting.DABTRAM = NA,
         day10_adapting.COCL2 = ifelse(fatepotential_COCL2_d10_w5 > adapting_thres, 'adapting', 'non-adapting'),
         day10_adapting.CIS = NA)

metadat.day10_CIS <- metadat.day10_CIS %>%
  mutate(day10_adapting.DABTRAM = NA,
         day10_adapting.COCL2 = NA,
         day10_adapting.CIS = ifelse(fatepotential_CIS_d10_w5 > adapting_thres, 'adapting', 'non-adapting'))

metadat.week5_DABTRAM <- metadat.week5_DABTRAM %>% 
  mutate(day10_adapting.DABTRAM = NA,
         day10_adapting.COCL2 = NA,
         day10_adapting.CIS = NA)

metadat.week5_COCL2 <- metadat.week5_COCL2 %>%
  mutate(day10_adapting.DABTRAM = NA,
         day10_adapting.COCL2 = NA,
         day10_adapting.CIS = NA)

metadat.week5_CIS <- metadat.week5_CIS %>%
  mutate(day10_adapting.DABTRAM = NA,
         day10_adapting.COCL2 = NA,
         day10_adapting.CIS = NA)

metadat.new <- rbind(metadat.day0, metadat.day10_DABTRAM, metadat.day10_COCL2, metadat.day10_CIS, metadat.week5_DABTRAM, metadat.week5_COCL2, metadat.week5_CIS)

all_data <- AddMetaData(all_data, metadat.new)

# =============================================================================
# Subset to day0 and day10 DABTRAM cells
# =============================================================================
all_data.day0_day10DABTRAM <- subset(all_data, (dataset == 'day0') | (dataset == 'day10_DABTRAM'))
all_data.day0_day10DABTRAM@active.assay <- 'Saver'

dim(all_data.day0_day10DABTRAM)[2] == nrow(metadat.day0) + nrow(metadat.day10_DABTRAM)

# process data
all_data.day0_day10DABTRAM <- NormalizeData(all_data.day0_day10DABTRAM)
all_data.day0_day10DABTRAM <- ScaleData(all_data.day0_day10DABTRAM)
all_data.day0_day10DABTRAM <- RunPCA(all_data.day0_day10DABTRAM, features = VariableFeatures(object = all_data.day0_day10DABTRAM))
DimPlot(all_data.day0_day10DABTRAM, reduction = "pca", group.by = "day10_adapting.DABTRAM")

all_data.day0_day10DABTRAM <- FindNeighbors(all_data.day0_day10DABTRAM, dims = 1:10)
all_data.day0_day10DABTRAM <- FindClusters(all_data.day0_day10DABTRAM, resolution = 0.5)
all_data.day0_day10DABTRAM <- RunUMAP(all_data.day0_day10DABTRAM, dims = 1:10)

DimPlot(all_data.day0_day10DABTRAM, reduction = "umap", group.by = "day10_adapting.DABTRAM")

# clean up
all_data.day0_day10DABTRAM[['fasttopic_COCL2']] <- NULL
all_data.day0_day10DABTRAM[['fasttopic_CIS']] <- NULL
all_data.day0_day10DABTRAM[['ft.COCL2.umap']] <- NULL
all_data.day0_day10DABTRAM[['ft.CIS.umap']] <- NULL

# save data
save(date_of_run, session_info,
     all_data.day0_day10DABTRAM,
     file = paste0(out_dir, 'all_data.day0_day10DABTRAM.RData'))

metadat.all_data.day0_day10DABTRAM <- all_data.day0_day10DABTRAM@meta.data
table(metadat.all_data.day0_day10DABTRAM$day10_adapting.DABTRAM)

# =============================================================================
# Subset to day0 and day10 COCL2 cells
# =============================================================================
all_data.day0_day10COCL2 <- subset(all_data, (dataset == 'day0') | (dataset == 'day10_COCL2'))
all_data.day0_day10COCL2@active.assay <- 'Saver'

dim(all_data.day0_day10COCL2)[2] == nrow(metadat.day0) + nrow(metadat.day10_COCL2)

# process data
all_data.day0_day10COCL2 <- NormalizeData(all_data.day0_day10COCL2)
all_data.day0_day10COCL2 <- ScaleData(all_data.day0_day10COCL2)
all_data.day0_day10COCL2 <- RunPCA(all_data.day0_day10COCL2, features = VariableFeatures(object = all_data.day0_day10COCL2))
DimPlot(all_data.day0_day10COCL2, reduction = "pca", group.by = "day10_adapting.COCL2")

all_data.day0_day10COCL2 <- FindNeighbors(all_data.day0_day10COCL2, dims = 1:10)
all_data.day0_day10COCL2 <- FindClusters(all_data.day0_day10COCL2, resolution = 0.5)
all_data.day0_day10COCL2 <- RunUMAP(all_data.day0_day10COCL2, dims = 1:10)

DimPlot(all_data.day0_day10COCL2, reduction = "umap", group.by = "day10_adapting.COCL2")

# clean up
all_data.day0_day10COCL2[['fasttopic_DABTRAM']] <- NULL
all_data.day0_day10COCL2[['fasttopic_CIS']] <- NULL
all_data.day0_day10COCL2[['ft.DABTRAM.umap']] <- NULL
all_data.day0_day10COCL2[['ft.CIS.umap']] <- NULL

# save data
save(date_of_run, session_info,
     all_data.day0_day10COCL2,
     file = paste0(out_dir, 'all_data.day0_day10COCL2.RData'))

metadat.all_data.day0_day10COCL2 <- all_data.day0_day10COCL2@meta.data
table(all_data.day0_day10COCL2$day10_adapting.COCL2)

# =============================================================================
# Subset to day0 and day10 CIS cells
# =============================================================================
all_data.day0_day10CIS <- subset(all_data, (dataset == 'day0') | (dataset == 'day10_CIS'))
all_data.day0_day10CIS@active.assay <- 'Saver'

dim(all_data.day0_day10CIS)[2] == nrow(metadat.day0) + nrow(metadat.day10_CIS)

# process data
all_data.day0_day10CIS <- NormalizeData(all_data.day0_day10CIS)
all_data.day0_day10CIS <- ScaleData(all_data.day0_day10CIS)
all_data.day0_day10CIS <- RunPCA(all_data.day0_day10CIS, features = VariableFeatures(object = all_data.day0_day10CIS))
DimPlot(all_data.day0_day10CIS, reduction = "pca", group.by = "day10_adapting.CIS")

all_data.day0_day10CIS <- FindNeighbors(all_data.day0_day10CIS, dims = 1:10)
all_data.day0_day10CIS <- FindClusters(all_data.day0_day10CIS, resolution = 0.5)
all_data.day0_day10CIS <- RunUMAP(all_data.day0_day10CIS, dims = 1:10)

# DimPlot(all_data.day0_day10CIS, reduction = "ft.CIS.umap", split.by = 'day10_adapting.CIS')
# FeaturePlot(all_data.day0_day10CIS, 'fatepotential_CIS_d10_w5', reduction = "umap")

# clean up
all_data.day0_day10CIS[['fasttopic_DABTRAM']] <- NULL
all_data.day0_day10CIS[['fasttopic_COCL2']] <- NULL
all_data.day0_day10CIS[['ft.DABTRAM.umap']] <- NULL
all_data.day0_day10CIS[['ft.COCL2.umap']] <- NULL

# save data
save(date_of_run, session_info,
     all_data.day0_day10CIS,
     file = paste0(out_dir, 'all_data.day0_day10CIS.RData'))

metadat.all_data.day0_day10CIS <- all_data.day0_day10CIS@meta.data
table(metadat.all_data.day0_day10CIS$day10_adapting.CIS)

