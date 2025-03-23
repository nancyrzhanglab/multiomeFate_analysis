rm(list = ls())

set.seed(123)

library(tidyverse)
library(ggpubr)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'

treatment <- 'DABTRAM'

remove_unassigned_cells <- TRUE

date_of_run <- Sys.time()
session_info <- devtools::session_info()

# =============================================================================
# reading data
# =============================================================================

load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))
load(paste0(data_dir, 'Writeup10a_data_chromVar_day10_DABTRAM.RData'))

all_data@misc <- all_data_fatepotential
all_data[["day10_DABTRAM"]] <- all_data_chromVar_day10_DABTRAM

# remove cells with no lineage
if(remove_unassigned_cells) {
  print("Removing cells with no assigned lineage")
  all_data$keep <- !is.na(all_data$assigned_lineage)
  if(any(!all_data$keep)){
    print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
    all_data <- subset(all_data, keep == TRUE)
  }
}

# ==============================================================================
# Calculate fate potential variance
# ==============================================================================
metadat <- all_data@meta.data
metadat$cell_id <- rownames(metadat)

# DABTRAM
fp.d10_w5.DABTRAM <- as.data.frame(all_data@misc[["fatepotential_DABTRAM_d10_w5"]][['cell_imputed_score']])
colnames(fp.d10_w5.DABTRAM) <- 'fatepotential_DABTRAM_d10_w5'
fp.d10_w5.DABTRAM$cell_id <- rownames(fp.d10_w5.DABTRAM)
fp.d10_w5.DABTRAM <- merge(fp.d10_w5.DABTRAM, metadat[, c('cell_id', 'assigned_lineage')], by = 'cell_id')

var.fp.d10_w5.DABTRAM <- fp.d10_w5.DABTRAM %>% 
  group_by(assigned_lineage) %>%
  summarise(variance = var(fatepotential_DABTRAM_d10_w5, na.rm = TRUE),
            n_cells = n())
# var.fp.d10_w5.DABTRAM <- var.fp.d10_w5.DABTRAM[var.fp.d10_w5.DABTRAM$n_cells > lin.size.thres, ]

# ==============================================================================
# Calculate mean gene expression
# ==============================================================================

mat.cv <- as.data.frame(t(all_data@assays[["day10_DABTRAM"]]@data))
mat.cv$cell_id <- rownames(mat.cv)
mat.cv <- merge(mat.cv, metadat[, c('cell_id', 'assigned_lineage')], by = 'cell_id')

# DABTRAM
metadat.day10_DABTRAM <- metadat[metadat$dataset == 'day10_DABTRAM', ]

mean.cv.day10_DABTRAM <- mat.cv %>% 
  filter(cell_id %in% metadat.day10_DABTRAM$cell_id) %>%
  group_by(assigned_lineage) %>%
  select(-cell_id) %>%
  summarise_all(mean, na.rm = TRUE)


df.day10_DABTRAM <- merge(mean.cv.day10_DABTRAM, var.fp.d10_w5.DABTRAM, by = 'assigned_lineage') %>% drop_na()

# ==============================================================================
# correlate variance and chromVar means
# ==============================================================================
tfs <- rownames(all_data@assays[["day10_DABTRAM"]]@data)

# DABTRAM
chromVar_cor_vec.day10_DABTRAM <- sapply(tfs, function(j){
  res <- stats::cor.test(df.day10_DABTRAM[['variance']], df.day10_DABTRAM[,j],
                         method = "spearman")
  c(res$estimate, res$p.value)
})

chromVar_cor_vec.day10_DABTRAM <- as.data.frame(t(chromVar_cor_vec.day10_DABTRAM))
colnames(chromVar_cor_vec.day10_DABTRAM) <- c("correlation", "p.value")
rownames(chromVar_cor_vec.day10_DABTRAM) <- tfs

getAndSortTable <- function(df) {
  cor_vec <- as.data.frame(df) %>% drop_na()
  cor_vec$p.value <- as.numeric(cor_vec$p.value)
  cor_vec$p_adj <- p.adjust(cor_vec$p.value, method = 'BH')
  cor_vec$tf <- rownames(cor_vec)
  
  cor_vec$correlation <- as.numeric(cor_vec$correlation)
  cor_vec <- cor_vec[order(cor_vec$correlation, decreasing = FALSE),]
  
  return(cor_vec)
}

chromVar_cor_vec.day10_DABTRAM <- getAndSortTable(chromVar_cor_vec.day10_DABTRAM)

# =============================================================================
# Plot
# =============================================================================
tfs_toplot <- tfs[grepl('JUN', tfs) | grepl('FOS', tfs) | grepl('SOX10', tfs) | grepl('MITF', tfs) | grepl('TEAD', tfs) | grepl('STAT', tfs)]
tfs_toplot <- tfs_toplot[!grepl('(var.2)', tfs_toplot)]

keyTFs.dabtram <- chromVar_cor_vec.day10_DABTRAM[chromVar_cor_vec.day10_DABTRAM$tf %in% tfs_toplot, 'correlation']


chromVar_cor_vec.day10_DABTRAM$order.DABTRAM <- 1:nrow(chromVar_cor_vec.day10_DABTRAM)
p1.TF <- ggplot(chromVar_cor_vec.day10_DABTRAM, aes(x = order.DABTRAM, y = correlation)) +
  geom_point(size = 1) +
  geom_point(data = chromVar_cor_vec.day10_DABTRAM[chromVar_cor_vec.day10_DABTRAM$tf %in% tfs_toplot, ], color = 'red') +
  ggrepel::geom_text_repel(data = subset(chromVar_cor_vec.day10_DABTRAM, tf %in% tfs_toplot), 
                           aes(label = tf), nudge_x = -15, nudge_y = 0.05, max.overlaps = 20) +
  xlab('TF') +
  ylab('Correlation') +
  ggtitle('Day10 DABTRAM') +
  ylim(-1, 1) +
  theme_Publication()
p1.TF
