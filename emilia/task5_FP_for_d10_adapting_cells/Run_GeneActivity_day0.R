rm(list=ls())
library(multiomeFate)
library(tidyverse)
library(ggplot2)
library(GGally)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'

remove_unassigned_cells <- TRUE

# =============================================================================
# Read data
# =============================================================================
load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_atac.RData'))

all_data[['atac']] <- all_data_atac

# remove cells with no lineage
# if(remove_unassigned_cells) {
#   print("Removing cells with no assigned lineage")
#   all_data$keep <- !is.na(all_data$assigned_lineage)
#   if(any(!all_data$keep)){
#     print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
#     all_data <- subset(all_data, keep == TRUE)
#   }
# }
# 
# metadat <- all_data@meta.data

# =============================================================================
# Subset to day0 and run GA
# =============================================================================
all_data_day0 <- subset(all_data, dataset == 'day0')

all_data_day0@assays[["atac"]]@fragments[[1]]@path <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/Raw_and_Processed/d0/atac_fragments.tsv.gz'

all_data_day0@active.assay <- 'atac'

gene.activities <- GeneActivity(all_data_day0)
saveRDS(gene.activities, '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/Raw_and_Processed/d0/gene_activities.rds')
