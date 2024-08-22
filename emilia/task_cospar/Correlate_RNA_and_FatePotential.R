library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

# ==============================================================================
# Read data
# ==============================================================================

in_dir <- '~/Dropbox/Thesis/Lineage_trace/data/Watermelon/'
out_dir <- '~/Dropbox/Thesis/Lineage_trace/outputs/task5_cospar/Watermelon/'

## Read fate potentials
tp_early <- "0"
tp_late <- "14"

load(paste0(out_dir, 'PC9_time_course_d', tp_early, "_d", tp_late, ".RData"))
cell_imputed_score <- final_fit$cell_imputed_score

## Read RNA
load(paste0(in_dir, 'PC9_time_course_fasttopics.RData'))
seurat_subset <- subset(seurat_object, subset = time_point == tp_early)
metadat <- seurat_subset@meta.data
metadat$cell_barcode <- rownames(metadat)

mat.RNA <- seurat_subset@assays[["RNA"]]@scale.data

## Read raw correlations
supp2g <- read.csv('~/Downloads/watermelon_fig2g.csv', row.names =1)
supp2g$gene <- rownames(supp2g)

# ==============================================================================
# Calculate correlation data
# ==============================================================================

mat.RNA.use <- mat.RNA[, names(cell_imputed_score)]
mat.RNA.use <- t(mat.RNA.use)
print(all(rownames(mat.RNA.use) == names(cell_imputed_score)))
rna_cor_vec <- sapply(1:ncol(mat.RNA.use), function(j){
  res <- stats::cor.test(cell_imputed_score, mat.RNA.use[,j],
                         alternative = "two.sided",
                         method = "pearson")
  c(res$estimate, res$p.value)
})
rna_cor_vec <- as.data.frame(t(rna_cor_vec))
colnames(rna_cor_vec) <- c("correlation", "p.value")
rownames(rna_cor_vec) <- colnames(mat.RNA.use)

hist(rna_cor_vec$correlation, breaks = 100)

rna_cor_vec <- rna_cor_vec %>% drop_na()
rna_cor_vec$gene <- rownames(rna_cor_vec)
rna_cor_vec <- rna_cor_vec[order(rna_cor_vec$correlation), ]
rna_cor_vec$order <- seq(1: nrow(rna_cor_vec))

bottom10 <- head(rna_cor_vec, 10)
top10 <- tail(rna_cor_vec, 10)
inPaper <- rna_cor_vec[c('GPX2', 'AKR1B10', 'AKR1C2', 'UCHL1', 'GSTM3'), ]
top_shared <- rna_cor_vec[c('GHR', 'PRKG2', 'BNC2', 'ADCY2', 'IGSF11', 'B4GALT5', 'IL16', 'PDE3B','SLCO5A1','NEDD4L', 'TMEM163',
                            'FRMD4B', 'PMEL', 'NRG3', 'MYO1D','SEMA5A', 'BACE2', 'PLCB4', 'KCNQ5', 'ST6GALNAC3', 'UBE3A', 'WDR33', 'DTD1', 'ACBD6', 'LRMDA', 'ILKAP', 'EEFSEC'), ]
bottom_shared <- rna_cor_vec[c('SPOCD1', 'ARL4C', 'MAMLD1', 'CAV1', 'COL8A1', 'F2R', 'PHLDB2',
                               'CSF1', 'ITGA3', 'INHBA', 'COL6A1', 'LYPD6B' ,'PDGFRB', 'BST2', 'DKK1'), ]

ggplot(rna_cor_vec, aes(x = order, y = correlation)) +
  geom_point() +
  geom_text_repel(data = bottom10, aes(label = gene), color = 'blue') +
  geom_text_repel(data = top10, aes(label = gene), color = 'red') +
  geom_text_repel(data = inPaper, aes(label = gene),nudge_y = 0.03, color = 'purple') +
  geom_text_repel(data = top_shared, aes(label = gene),nudge_y = 0.03, color = '#DA7297') +
  geom_text_repel(data = bottom_shared, aes(label = gene),nudge_y = -0.03, color = '#80C4E9') +
  theme_bw()

# ==============================================================================
# Compare with their data
# ==============================================================================

comp <- merge(supp2g, rna_cor_vec, by = 'gene')

ggplot(comp, aes(x = day_7, y = day_14)) +
  geom_point(alpha=0.5) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01)+
  theme_bw()

# ==============================================================================
# Fate potential distribution
# ==============================================================================
cell_imputed_score_df <- as.data.frame(cell_imputed_score)
cell_imputed_score_df$cell_barcode <- rownames(cell_imputed_score_df)
cell_imputed_score_df <- merge(cell_imputed_score_df, metadat, by = 'cell_barcode')

# by lineage size
lineage_size <- cell_imputed_score_df %>% 
  group_by(lineage_barcode) %>% 
  summarise(n_cells = n())
lineage_size <- lineage_size[order(lineage_size$n_cells), ]
lineage_size <- lineage_size[lineage_size$n_cells > 2, ]
largest_lbs <- tail(lineage_size, 10)
smallest_lbs <- head(lineage_size, 10)

lin_to_plot <- cell_imputed_score_df[cell_imputed_score_df$lineage_barcode %in% c(largest_lbs$lineage_barcode, smallest_lbs$lineage_barcode), ]
lin_to_plot$category <- ifelse(lin_to_plot$lineage_barcode %in% largest_lbs$lineage_barcode, 'Large', 'Small')

ggplot(lin_to_plot, aes(x = lineage_barcode, y = cell_imputed_score)) +
  geom_violin(aes(fill = category), scale = 'width') +
  geom_jitter(width = 0.1) +
  theme_bw()

# by mean
mean_fp <- cell_imputed_score_df %>% 
  group_by(lineage_barcode) %>% 
  summarise(mean_fp = mean(cell_imputed_score),
            n_cells = n())
mean_fp <- mean_fp[order(mean_fp$mean_fp), ]
mean_fp <- mean_fp[mean_fp$n_cells > 2, ]
top_lbs <- tail(mean_fp, 15)
bottom_lbs <- head(mean_fp, 15)

lin_to_plot2 <- cell_imputed_score_df[cell_imputed_score_df$lineage_barcode %in% c(top_lbs$lineage_barcode, bottom_lbs$lineage_barcode), ]
lin_to_plot2$category <- ifelse(lin_to_plot2$lineage_barcode %in% top_lbs$lineage_barcode, 'High mean', 'Low mean')
lin_to_plot2 <- lin_to_plot2 %>%
                  mutate(lineage_barcode = fct_reorder(lineage_barcode, cell_imputed_score))

ggplot(lin_to_plot2, aes(x = lineage_barcode, y = cell_imputed_score)) +
  geom_violin(aes(fill = category), scale = 'width') +
  stat_summary(fun = "median",
               geom = "crossbar",
               color = 'red') +
  geom_jitter(width = 0.1) +
  theme_bw()
