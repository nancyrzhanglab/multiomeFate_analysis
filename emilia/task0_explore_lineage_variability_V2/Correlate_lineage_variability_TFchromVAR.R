rm(list = ls())

set.seed(123)

library(tidyverse)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
results_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task0_explore_lineage_variability_V2/Feature_correlations/'
figure_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/figures/MultiomeFate_V3/Fig4/'

dataset_colors <- c(day0 = "gray",
                    day10_CIS = "#FBD08C",
                    day10_COCL2 = "#6DC49C",
                    day10_DABTRAM = "#9D85BE",
                    week5_CIS = "#C96D29",
                    week5_COCL2 = "#0F8241",
                    week5_DABTRAM = "#623594")

remove_unassigned_cells <- TRUE

date_of_run <- Sys.time()
session_info <- devtools::session_info()

# =============================================================================
# reading data
# =============================================================================

load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))
load(paste0(data_dir, 'Writeup10a_data_chromVar_day10_DABTRAM.RData'))
load(paste0(data_dir, 'Writeup10a_data_chromVar_day10_COCL2.RData'))
load(paste0(data_dir, 'Writeup10a_data_chromVar_day10_CIS.RData'))

all_data@misc <- all_data_fatepotential
all_data[["chromVar.day10_DABTRAM"]] <- all_data_chromVar_day10_DABTRAM
all_data[["chromVar.day10_COCL2"]] <- all_data_chromVar_day10_COCL2
all_data[["chromVar.day10_CIS"]] <- all_data_chromVar_day10_CIS

# remove cells with no lineage
if(remove_unassigned_cells) {
  print("Removing cells with no assigned lineage")
  all_data$keep <- !is.na(all_data$assigned_lineage)
  if(any(!all_data$keep)){
    print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
    all_data <- subset(all_data, keep == TRUE)
  }
}

lin_var.day10_DABTRAM <- read.csv(paste0(results_dir, 'day10_DABTRAM/lineage_variability_day10_DABTRAM_saver_sample.csv'))
lin_var.day10_COCL2 <- read.csv(paste0(results_dir, 'day10_COCL2/lineage_variability_day10_COCL2_saver_sample.csv'))
lin_var.day10_CIS <- read.csv(paste0(results_dir, 'day10_CIS/lineage_variability_day10_CIS_saver_sample.csv'))

# =============================================================================
# calcualte mean chromvar activity
# =============================================================================
metadat <- all_data@meta.data
metadat$cell_id <- rownames(metadat)

mat.cv.dabtram <- all_data@assays[["chromVar.day10_DABTRAM"]]@data
mat.cv.cocl2 <- all_data@assays[["chromVar.day10_COCL2"]]@data
mat.cv.cis <- all_data@assays[["chromVar.day10_CIS"]]@data

# dabtram
metadat.day10_DABTRAM <- metadat[metadat$dataset == 'day10_DABTRAM', ]
mat.cv.dabtram <- as.data.frame(t(mat.cv.dabtram))
mat.cv.dabtram$cell_id <- rownames(mat.cv.dabtram)
mat.cv.dabtram <- merge(mat.cv.dabtram, metadat.day10_DABTRAM[, c('cell_id', 'assigned_lineage')], by = 'cell_id')
mean.cv.dabtram <- mat.cv.dabtram %>%
  filter(cell_id %in% metadat.day10_DABTRAM$cell_id) %>%
  select(-cell_id) %>%
  group_by(assigned_lineage) %>%
  summarise_all(mean, na.rm = TRUE) %>%
  ungroup()
df.day10_DABTRAM <- merge(mean.cv.dabtram, lin_var.day10_DABTRAM, by = 'assigned_lineage') %>% drop_na()

# COCL2
metadat.day10_COCL2 <- metadat[metadat$dataset == 'day10_COCL2', ]
mat.cv.cocl2 <- as.data.frame(t(mat.cv.cocl2))
mat.cv.cocl2$cell_id <- rownames(mat.cv.cocl2)
mat.cv.cocl2 <- merge(mat.cv.cocl2, metadat.day10_COCL2[, c('cell_id', 'assigned_lineage')], by = 'cell_id')

mean.cv.cocl2 <- mat.cv.cocl2 %>% 
  filter(cell_id %in% metadat.day10_COCL2$cell_id) %>%
  group_by(assigned_lineage) %>%
  select(-cell_id) %>%
  summarise_all(mean, na.rm = TRUE)
df.day10_COCL2 <- merge(mean.cv.cocl2, lin_var.day10_COCL2, by = 'assigned_lineage') %>% drop_na()

# CIS
metadat.day10_CIS <- metadat[metadat$dataset == 'day10_CIS', ]
mat.cv.cis <- as.data.frame(t(mat.cv.cis))
mat.cv.cis$cell_id <- rownames(mat.cv.cis)
mat.cv.cis <- merge(mat.cv.cis, metadat.day10_CIS[, c('cell_id', 'assigned_lineage')], by = 'cell_id')

mean.cv.cis <- mat.cv.cis %>% 
  filter(cell_id %in% metadat.day10_CIS$cell_id) %>%
  group_by(assigned_lineage) %>%
  select(-cell_id) %>%
  summarise_all(mean, na.rm = TRUE)
df.day10_CIS <- merge(mean.cv.cis, lin_var.day10_CIS, by = 'assigned_lineage') %>% drop_na()

# ==============================================================================
# correlate variance and gene expression means
# ==============================================================================
tfs <- rownames(all_data@assays[["chromVar.day10_DABTRAM"]]@data)

# DABTRAM
linvar_cor_vec.day10_DABTRAM <- sapply(tfs, function(j){
  if(quantile(df.day10_DABTRAM[[j]], na.rm = T)[4] == 10) {
    return(c(NA, NA))
  }
  res <- stats::cor.test(df.day10_DABTRAM[['normalized_avg_eud_dist_by_shuffle']], df.day10_DABTRAM[,j],
                         method = "spearman")
  c(res$estimate, res$p.value)
})

linvar_cor_vec.day10_DABTRAM <- as.data.frame(t(linvar_cor_vec.day10_DABTRAM))
colnames(linvar_cor_vec.day10_DABTRAM) <- c("correlation", "p.value")
rownames(linvar_cor_vec.day10_DABTRAM) <- tfs

# COCL2
linvar_cor_vec.day10_COCL2 <- sapply(tfs, function(j){
  if(quantile(df.day10_COCL2[[j]], na.rm = T)[4] == 10) {
    return(c(NA, NA))
  }
  res <- stats::cor.test(df.day10_COCL2[['normalized_avg_eud_dist_by_shuffle']], df.day10_COCL2[,j],
                         method = "spearman")
  c(res$estimate, res$p.value)
})
linvar_cor_vec.day10_COCL2 <- as.data.frame(t(linvar_cor_vec.day10_COCL2))
colnames(linvar_cor_vec.day10_COCL2) <- c("correlation", "p.value")
rownames(linvar_cor_vec.day10_COCL2) <- tfs

# CIS
linvar_cor_vec.day10_CIS <- sapply(tfs, function(j){
  if(quantile(df.day10_CIS[[j]], na.rm = T)[4] == 10) {
    return(c(NA, NA))
  }
  res <- stats::cor.test(df.day10_CIS[['normalized_avg_eud_dist_by_shuffle']], df.day10_CIS[,j],
                         method = "spearman")
  c(res$estimate, res$p.value)
})
linvar_cor_vec.day10_CIS <- as.data.frame(t(linvar_cor_vec.day10_CIS))
colnames(linvar_cor_vec.day10_CIS) <- c("correlation", "p.value")
rownames(linvar_cor_vec.day10_CIS) <- tfs


cor.obj <- list(
  day10_DABTRAM_linvar_cor_df = linvar_cor_vec.day10_DABTRAM,
  day10_COCL2_linvar_cor_df = linvar_cor_vec.day10_COCL2,
  day10_CIS_linvar_cor_df = linvar_cor_vec.day10_CIS
)

# save the correlation objects
saveRDS(cor.obj, file = paste0(results_dir, "day10_linvar_TFchromVAR_cor_df.rds"))
