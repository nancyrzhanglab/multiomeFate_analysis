rm(list=ls())
library(Seurat)
library(clusterProfiler)
library(UCell)
library(tidyverse)
library(ggpubr)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
ref_dir <- '/Users/emiliac/Dropbox/Epi Evolution Paper/Data/Gavish Clinical Analysis/Resources/'
out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'

remove_unassigned_cells <- TRUE
treatment <- 'DABTRAM'


# ==============================================================================
# Read data
# ==============================================================================
load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_chromVar_week5_DABTRAM.RData'))

all_data[['chromVar.week5_DABTRAM']] <- all_data_chromVar_week5_DABTRAM



# remove cells with no lineage
if(remove_unassigned_cells) {
  print("Removing cells with no assigned lineage")
  all_data$keep <- !is.na(all_data$assigned_lineage)
  if(any(!all_data$keep)){
    print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
    all_data <- subset(all_data, keep == TRUE)
  }
}

metadat <- all_data@meta.data
metadat$cell_id <- rownames(metadat)

# ==============================================================================
# Get high fp cells from d10 and cells in w5
# ==============================================================================

# d0
metadat.d0 <- metadat %>% filter(dataset == 'day0')

# w5
metadat.DABTRAM <- read.csv(paste0(data_dir, 'metadat.all_data_dabtram_recluster.csv'), row.names = 1)
metadat.week5_DABTRAM <- metadat.DABTRAM[metadat.DABTRAM$dataset == 'week5_DABTRAM', ]
metadat.week5_DABTRAM.clust0 <- metadat.week5_DABTRAM %>% filter(seurat_clusters == 0)
metadat.week5_DABTRAM.clust3 <- metadat.week5_DABTRAM %>% filter(seurat_clusters == 3)

# ==============================================================================
# DE
# ==============================================================================

all_data.week5_DABTRAM <- subset(all_data, dataset == 'week5_DABTRAM')

pvalue_list <- lapply(all_data.week5_DABTRAM[["chromVar.week5_DABTRAM"]]@var.features, function(tf){
  tf_mat <- as.matrix(all_data.week5_DABTRAM[["chromVar.week5_DABTRAM"]]@data)
  tf_name <- unname(unlist(tf))
  x_vec <- tf_mat[tf_name, rownames(metadat.week5_DABTRAM.clust0)]
  y_vec <- tf_mat[tf_name, rownames(metadat.week5_DABTRAM.clust3)]
  
  # remove NAs
  x_vec <- x_vec[!is.na(x_vec)]
  y_vec <- y_vec[!is.na(y_vec)]
  
  if(diff(range(x_vec, na.rm = T)) <= 1e-6 || diff(range(y_vec, na.rm = T)) <= 1e-6 ) return(NA)
  
  test_res <- stats::wilcox.test(x = x_vec, y = y_vec)
  list(stat = mean(x_vec, na.rm = T) - mean(y_vec, na.rm = T), pvalue = test_res$p.value)
})
names(pvalue_list) <- all_data.week5_DABTRAM[["chromVar.week5_DABTRAM"]]@var.features
pvalue_list <- pvalue_list[sapply(pvalue_list, function(x){!all(is.na(x))})]
pvalue_vec <- sapply(pvalue_list, function(x){x$pvalue})
names(pvalue_vec) <- names(pvalue_list)
pvalue_vec <- stats::p.adjust(pvalue_vec, method = "BH")

diff_vec <- sapply(pvalue_list, function(x){x$stat})
logpval_vec <- sapply(pvalue_list, function(x){-log10(x$pvalue)})


df <- data.frame(difference = diff_vec,
                 pvalue_adj =  pvalue_vec,
                 log10pval = logpval_vec,
                 name = names(logpval_vec))
# ==============================================================================
# Plot
# ==============================================================================
# df$log10pval <- ifelse(df$log10pval > 100, 100, df$log10pval)
ggplot(df, aes(x = difference, y = log10pval)) +
  geom_point() +
  geom_hline(yintercept=min(logpval_vec[which(pvalue_vec <= 0.05)]), linetype="dashed", 
             color = "gray", linewidth=2) +
  geom_vline(xintercept=0, linetype="dashed",
             color = "gray", linewidth=2) 


tf_mat <- as.matrix(all_data.week5_DABTRAM[["chromVar.week5_DABTRAM"]]@data)
x_vec <- tf_mat['FOSL2', rownames(metadat.week5_DABTRAM.clust0)]
y_vec <- tf_mat['FOSL2', rownames(metadat.week5_DABTRAM.clust3)]

df2 <- data.frame(x = x_vec)
df2$group <- 'clust0'
df3 <- data.frame(x = y_vec)
df3$group <- 'clust3'
df4 <- rbind(df2, df3)

ggplot(df4, aes(x = x, fill = group)) +
  geom_density(alpha = 0.5) +
  theme_minimal()
