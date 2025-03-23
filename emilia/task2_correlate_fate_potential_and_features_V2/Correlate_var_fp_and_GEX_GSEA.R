rm(list = ls())
library(Seurat)
library(tidyverse)
library(data.table)
library(GSEABase)
library(Biobase) #base functions for bioconductor; required by GSEABase
library(clusterProfiler) # provides a suite of tools for functional enrichment analysis
library(enrichplot) # great for making the standard GSEA enrichment plots
library(ggplot2)
library(ggpubr)

theme_Publication<- function(base_size=12, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            plot.subtitle = element_text(face = "bold", hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold"),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line.x = element_line(colour="black"),
            axis.line.y = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="white"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.box.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "mm"),
            legend.direction = "vertical",
            legend.key.size= unit(0.5, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"),
            strip.background=element_rect(colour="#F0F0F0",fill="#F0F0F0"),
            strip.text = element_text(face="bold")
    ))
}


data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
result_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/'

remove_unassigned_cells <- TRUE
lin.size.thres <- 3

# ==============================================================================
# Read data general
# ==============================================================================

load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))
load(paste0(data_dir, 'Writeup10a_data_saver.RData'))

all_data@misc <- all_data_fatepotential
all_data[["Saver"]] <- all_data_saver

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
var.fp.d10_w5.DABTRAM <- var.fp.d10_w5.DABTRAM[var.fp.d10_w5.DABTRAM$n_cells > lin.size.thres, ]

hist(var.fp.d10_w5.DABTRAM$variance, breaks = 50, main = 'DABTRAM', xlab = 'Variance of fate potential')

# COCL2
fp.d10_w5.COCL2 <- as.data.frame(all_data@misc[["fatepotential_COCL2_d10_w5"]][['cell_imputed_score']])
colnames(fp.d10_w5.COCL2) <- 'fatepotential_COCL2_d10_w5'
fp.d10_w5.COCL2$cell_id <- rownames(fp.d10_w5.COCL2)
fp.d10_w5.COCL2 <- merge(fp.d10_w5.COCL2, metadat[, c('cell_id', 'assigned_lineage')], by = 'cell_id')

var.fp.d10_w5.COCL2 <- fp.d10_w5.COCL2 %>% 
  group_by(assigned_lineage) %>%
  summarise(variance = var(fatepotential_COCL2_d10_w5, na.rm = TRUE),
            n_cells = n())
var.fp.d10_w5.COCL2 <- var.fp.d10_w5.COCL2[var.fp.d10_w5.COCL2$n_cells > lin.size.thres, ]

# CIS
fp.d10_w5.CIS <- as.data.frame(all_data@misc[["fatepotential_CIS_d10_w5"]][['cell_imputed_score']])
colnames(fp.d10_w5.CIS) <- 'fatepotential_CIS_d10_w5'
fp.d10_w5.CIS$cell_id <- rownames(fp.d10_w5.CIS)
fp.d10_w5.CIS <- merge(fp.d10_w5.CIS, metadat[, c('cell_id', 'assigned_lineage')], by = 'cell_id')

var.fp.d10_w5.CIS <- fp.d10_w5.CIS %>% 
  group_by(assigned_lineage) %>%
  summarise(variance = var(fatepotential_CIS_d10_w5, na.rm = TRUE),
            n_cells = n())
var.fp.d10_w5.CIS <- var.fp.d10_w5.CIS[var.fp.d10_w5.CIS$n_cells > lin.size.thres, ]

quantile(var.fp.d10_w5.CIS$n_cells)
# ==============================================================================
# Calculate mean gene expression
# ==============================================================================

mat.rna <- as.data.frame(t(all_data@assays[["Saver"]]@data))
mat.rna$cell_id <- rownames(mat.rna)
mat.rna <- merge(mat.rna, metadat[, c('cell_id', 'assigned_lineage')], by = 'cell_id')

# DABTRAM
metadat.day10_DABTRAM <- metadat[metadat$dataset == 'day10_DABTRAM', ]

mean.rna.day10_DABTRAM <- mat.rna %>% 
  filter(cell_id %in% metadat.day10_DABTRAM$cell_id) %>%
  group_by(assigned_lineage) %>%
  select(-cell_id) %>%
  summarise_all(mean, na.rm = TRUE)


df.day10_DABTRAM <- merge(mean.rna.day10_DABTRAM, var.fp.d10_w5.DABTRAM, by = 'assigned_lineage') %>% drop_na()

# COCL2
metadat.day10_COCL2 <- metadat[metadat$dataset == 'day10_COCL2', ]
mean.rna.day10_COCL2 <- mat.rna %>% 
  filter(cell_id %in% metadat.day10_COCL2$cell_id) %>%
  group_by(assigned_lineage) %>%
  select(-cell_id) %>%
  summarise_all(mean, na.rm = TRUE)
df.day10_COCL2 <- merge(mean.rna.day10_COCL2, var.fp.d10_w5.COCL2, by = 'assigned_lineage') %>% drop_na()

# CIS
metadat.day10_CIS <- metadat[metadat$dataset == 'day10_CIS', ]
mean.rna.day10_CIS <- mat.rna %>% 
  filter(cell_id %in% metadat.day10_CIS$cell_id) %>%
  group_by(assigned_lineage) %>%
  select(-cell_id) %>%
  summarise_all(mean, na.rm = TRUE)
df.day10_CIS <- merge(mean.rna.day10_CIS, var.fp.d10_w5.CIS, by = 'assigned_lineage') %>% drop_na()

## check mean scales with variance

var.rna.day10_DABTRAM <- mat.rna %>%
  filter(cell_id %in% metadat.day10_DABTRAM$cell_id) %>%
  group_by(assigned_lineage) %>%
  select(-cell_id) %>%
  summarise_all(var, na.rm = TRUE)

# lin size
lin.size.day10_DABTRAM <- mat.rna %>%
  filter(cell_id %in% metadat.day10_DABTRAM$cell_id) %>%
  group_by(assigned_lineage) %>%
  summarize(n_cells = n())

mean.rna.day10_DABTRAM.melt <- melt(mean.rna.day10_DABTRAM, id.vars = 'assigned_lineage')
colnames(mean.rna.day10_DABTRAM.melt) <- c('assigned_lineage', 'gene', 'mean_expression')

var.rna.day10_DABTRAM.melt <- melt(var.rna.day10_DABTRAM, id.vars = 'assigned_lineage')
colnames(var.rna.day10_DABTRAM.melt) <- c('assigned_lineage', 'gene', 'variance')

df.comp <- merge(mean.rna.day10_DABTRAM.melt, var.rna.day10_DABTRAM.melt, by = c('assigned_lineage', 'gene'))

df.comp.one_gene <- df.comp[df.comp$gene == 'FN1', ]
df.comp.one_gene <- merge(df.comp.one_gene, lin.size.day10_DABTRAM, by = 'assigned_lineage')
ggplot(df.comp.one_gene, aes(x = mean_expression, y = sqrt(variance))) +
  geom_point(aes(size = n_cells), shape = 21, fill = 'gray') +
  stat_cor() +
  ggtitle('day10 DABTRAM: FN1') +
  theme_bw()

# COCL2
var.rna.day10_COCL2 <- mat.rna %>%
  filter(cell_id %in% metadat.day10_COCL2$cell_id) %>%
  group_by(assigned_lineage) %>%
  select(-cell_id) %>%
  summarise_all(var, na.rm = TRUE)

mean.rna.day10_COCL2.melt <- melt(mean.rna.day10_COCL2, id.vars = 'assigned_lineage')
colnames(mean.rna.day10_COCL2.melt) <- c('assigned_lineage', 'gene', 'mean_expression')

var.rna.day10_COCL2.melt <- melt(var.rna.day10_COCL2, id.vars = 'assigned_lineage')
colnames(var.rna.day10_COCL2.melt) <- c('assigned_lineage', 'gene', 'variance')

df.comp <- merge(mean.rna.day10_COCL2.melt, var.rna.day10_COCL2.melt, by = c('assigned_lineage', 'gene'))

df.comp.one_gene <- df.comp[df.comp$gene == 'DKK1', ]

# lin size
lin.size.day10_COCL2 <- mat.rna %>%
  filter(cell_id %in% metadat.day10_COCL2$cell_id) %>%
  group_by(assigned_lineage) %>%
  summarize(n_cells = n())

df.comp.one_gene <- merge(df.comp.one_gene, lin.size.day10_COCL2, by = 'assigned_lineage')
ggplot(df.comp.one_gene, aes(x = mean_expression, y = sqrt(variance))) +
  geom_point(aes(size = n_cells), shape = 21, fill = 'gray') +
  stat_cor() +
  ggtitle('day10 COCL2: FN1') +
  theme_bw()

# CIS
# var.rna.day10_CIS <- mat.rna %>%
#   filter(cell_id %in% metadat.day10_CIS$cell_id) %>%
#   group_by(assigned_lineage) %>%
#   select(-cell_id) %>%
#   summarise_all(var, na.rm = TRUE)
# 
# mean.rna.day10_CIS.melt <- melt(mean.rna.day10_CIS, id.vars = 'assigned_lineage')
# colnames(mean.rna.day10_CIS.melt) <- c('assigned_lineage', 'gene', 'mean_expression')
# 
# var.rna.day10_CIS.melt <- melt(var.rna.day10_CIS, id.vars = 'assigned_lineage')
# colnames(var.rna.day10_CIS.melt) <- c('assigned_lineage', 'gene', 'variance')
# 
# df.comp <- merge(mean.rna.day10_CIS.melt, var.rna.day10_CIS.melt, by = c('assigned_lineage', 'gene'))
# 
# df.comp.one_gene <- df.comp[df.comp$gene == 'FN1', ]
# colnames(df.comp.one_gene)[4] <- 'variance.gene'
# df.comp.one_gene <- merge(df.comp.one_gene, var.fp.d10_w5.CIS, by = c('assigned_lineage'))
# 
# 
# ggplot(df.comp.one_gene, aes(x = mean_expression, y = sqrt(variance))) +
#   geom_point() +
#   stat_cor()
# 
# ggplot(df.comp.one_gene, aes(x = n_cells, y = mean_expression)) +
#   geom_point() +
#   stat_cor()

# ==============================================================================
# correlate variance and gene expression means
# ==============================================================================
genes <- rownames(all_data@assays[["Saver"]]@data)

# DABTRAM
saver_cor_vec.day10_DABTRAM <- sapply(genes, function(j){
  if(quantile(df.day10_DABTRAM[[j]], na.rm = T)[4] == 10) {
    return(c(NA, NA))
  }
  res <- stats::cor.test(df.day10_DABTRAM[['variance']], df.day10_DABTRAM[,j],
                         method = "spearman")
  c(res$estimate, res$p.value)
})

saver_cor_vec.day10_DABTRAM <- as.data.frame(t(saver_cor_vec.day10_DABTRAM))
colnames(saver_cor_vec.day10_DABTRAM) <- c("correlation", "p.value")
rownames(saver_cor_vec.day10_DABTRAM) <- genes

# COCL2
saver_cor_vec.day10_COCL2 <- sapply(genes, function(j){
  if(quantile(df.day10_COCL2[[j]], na.rm = T)[4] == 10) {
    return(c(NA, NA))
  }
  res <- stats::cor.test(df.day10_COCL2[['variance']], df.day10_COCL2[,j],
                         method = "spearman")
  c(res$estimate, res$p.value)
})
saver_cor_vec.day10_COCL2 <- as.data.frame(t(saver_cor_vec.day10_COCL2))
colnames(saver_cor_vec.day10_COCL2) <- c("correlation", "p.value")
rownames(saver_cor_vec.day10_COCL2) <- genes

# CIS
saver_cor_vec.day10_CIS <- sapply(genes, function(j){
  if(quantile(df.day10_CIS[[j]], na.rm = T)[4] == 10) {
    return(c(NA, NA))
  }
  res <- stats::cor.test(df.day10_CIS[['variance']], df.day10_CIS[,j],
                         method = "spearman")
  c(res$estimate, res$p.value)
})
saver_cor_vec.day10_CIS <- as.data.frame(t(saver_cor_vec.day10_CIS))
colnames(saver_cor_vec.day10_CIS) <- c("correlation", "p.value")
rownames(saver_cor_vec.day10_CIS) <- genes


# ==============================================================================
# Plot gene and fp var correlation
# ==============================================================================

df.day10_COCL2.one_gene <- df.day10_COCL2[, c('assigned_lineage', 'FAP', 'variance')]
ggplot(df.day10_COCL2.one_gene, aes(x = FAP, y = sqrt(variance))) +
  geom_point() +
  stat_cor()
quantile(df.day10_COCL2.one_gene$CD44, 0.90)

df.day10_CIS.one_gene <- df.day10_CIS[, c('assigned_lineage', 'FN1', 'variance')]
ggplot(df.day10_CIS.one_gene, aes(x = FN1, y = sqrt(variance))) +
  geom_point() +
  stat_cor()
quantile(df.day10_CIS.one_gene$FN1)

# ==============================================================================
# Read signatures
# ==============================================================================
ref_dir <- '/Users/emiliac/Dropbox/Thesis/resources/GSEA_pathways/'
hallmark <- read.gmt(paste0(ref_dir, "h.all.v2024.1.Hs.symbols.gmt"))
reactome <- read.gmt(paste0(ref_dir, "c2.cp.reactome.v2024.1.Hs.symbols.gmt"))

# ==============================================================================
# GSEA
# ==============================================================================

getAndSortTable <- function(df) {
  cor_vec <- as.data.frame(df) %>% drop_na()
  cor_vec$p.value <- as.numeric(cor_vec$p.value)
  cor_vec$p_adj <- p.adjust(cor_vec$p.value, method = 'BH')
  cor_vec <- cor_vec[cor_vec$p_adj < 0.05, ]
  cor_vec$gene <- rownames(cor_vec)
  
  cor_vec$correlation <- as.numeric(cor_vec$correlation)
  cor_vec <- cor_vec[order(cor_vec$correlation, decreasing = TRUE),]
  
  return(cor_vec)
}

saver_cor_vec.day10_DABTRAM <- getAndSortTable(saver_cor_vec.day10_DABTRAM)
saver_cor_vec.day10_COCL2 <- getAndSortTable(saver_cor_vec.day10_COCL2)
saver_cor_vec.day10_CIS <- getAndSortTable(saver_cor_vec.day10_CIS)

# GSEA DABTRAM
gsea_input.day10.DABTRAM <- saver_cor_vec.day10_DABTRAM$correlation
names(gsea_input.day10.DABTRAM) <- saver_cor_vec.day10_DABTRAM$gene

set.seed(123)
GSEA_res.day10.DABTRAM <- GSEA(geneList = gsea_input.day10.DABTRAM, 
                               TERM2GENE = hallmark, 
                               pvalueCutoff = 0.2,
                               seed = T,
                               verbose = F)
GSEA_res.day10.DABTRAM.df <- as_tibble(GSEA_res.day10.DABTRAM@result)
GSEA_res.day10.DABTRAM.df <- GSEA_res.day10.DABTRAM.df[, c('ID', 'NES', 'p.adjust', 'qvalue', 'setSize')]
colnames(GSEA_res.day10.DABTRAM.df) <- c('ID', 'NES.DABTRAM', 'p.adjust.DABTRAM', 'qvalue.DABTRAM', 'setSize')

# GSEA COCL2
gsea_input.day10.COCL2 <- saver_cor_vec.day10_COCL2$correlation
names(gsea_input.day10.COCL2) <- saver_cor_vec.day10_COCL2$gene

GSEA_res.day10.COCL2 <- GSEA(geneList = gsea_input.day10.COCL2, 
                               TERM2GENE = hallmark, 
                               pvalueCutoff = 0.2,
                               seed = T,
                               verbose = F)
GSEA_res.day10.COCL2.df <- as_tibble(GSEA_res.day10.COCL2@result)
GSEA_res.day10.COCL2.df <- GSEA_res.day10.COCL2.df[, c('ID', 'NES', 'p.adjust', 'qvalue', 'setSize')]
colnames(GSEA_res.day10.COCL2.df) <- c('ID', 'NES.COCL2', 'p.adjust.COCL2', 'qvalue.COCL2', 'setSize')

# GSEA CIS
gsea_input.day10.CIS <- saver_cor_vec.day10_CIS$correlation
names(gsea_input.day10.CIS) <- saver_cor_vec.day10_CIS$gene

GSEA_res.day10.CIS <- GSEA(geneList = gsea_input.day10.CIS, 
                               TERM2GENE = hallmark, 
                               pvalueCutoff = 0.2,
                               seed = T,
                               verbose = F)
GSEA_res.day10.CIS.df <- as_tibble(GSEA_res.day10.CIS@result)
GSEA_res.day10.CIS.df <- GSEA_res.day10.CIS.df[, c('ID', 'NES', 'p.adjust', 'qvalue', 'setSize')]
colnames(GSEA_res.day10.CIS.df) <- c('ID', 'NES.CIS', 'p.adjust.CIS', 'qvalue.CIS', 'setSize')


saver_cor_vec.day10_CIS.1 <- merge(saver_cor_vec.day10_CIS, hallmark, by = 'gene')
# ==============================================================================
# Plot
# ==============================================================================
GSEA_res.day10.DABTRAM.df <- GSEA_res.day10.DABTRAM.df[order(GSEA_res.day10.DABTRAM.df$NES.DABTRAM, decreasing = TRUE),]
GSEA_res.day10.DABTRAM.df <- GSEA_res.day10.DABTRAM.df[GSEA_res.day10.DABTRAM.df$qvalue.DABTRAM < 0.05,]
GSEA_res.day10.DABTRAM.df$neg_log10_pval.DABTRAM <- -log10(GSEA_res.day10.DABTRAM.df$p.adjust.DABTRAM)
GSEA_res.day10.DABTRAM.df$NES.DABTRAM.abs <- abs(GSEA_res.day10.DABTRAM.df$NES.DABTRAM)

GSEA_res.day10.COCL2.df <- GSEA_res.day10.COCL2.df[order(GSEA_res.day10.COCL2.df$NES.COCL2, decreasing = TRUE),]
GSEA_res.day10.COCL2.df <- GSEA_res.day10.COCL2.df[GSEA_res.day10.COCL2.df$qvalue.COCL2 < 0.05,]
GSEA_res.day10.COCL2.df$neg_log10_pval.COCL2 <- -log10(GSEA_res.day10.COCL2.df$p.adjust.COCL2)
GSEA_res.day10.COCL2.df$NES.COCL2.abs <- abs(GSEA_res.day10.COCL2.df$NES.COCL2)

GSEA_res.day10.CIS.df <- GSEA_res.day10.CIS.df[order(GSEA_res.day10.CIS.df$NES.CIS, decreasing = TRUE),]
GSEA_res.day10.CIS.df <- GSEA_res.day10.CIS.df[GSEA_res.day10.CIS.df$qvalue.CIS < 0.05,]
GSEA_res.day10.CIS.df$neg_log10_pval.CIS <- -log10(GSEA_res.day10.CIS.df$p.adjust.CIS)
GSEA_res.day10.CIS.df$NES.CIS.abs <- abs(GSEA_res.day10.CIS.df$NES.CIS)

ggplot(GSEA_res.day10.DABTRAM.df, aes(x = NES.DABTRAM, y = reorder(ID, NES.DABTRAM), color = neg_log10_pval.DABTRAM, size = setSize)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'black') +
  scale_size_continuous(range = c(1, 8), breaks = c(20, 40, 60)) +
  labs(title = 'GSEA DABTRAM Day 10', x = '', y = '') +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggplot(GSEA_res.day10.DABTRAM.df, aes(x = 'NES.DABTRAM', y = reorder(ID, NES.DABTRAM))) +
  geom_point(aes(size = neg_log10_pval.DABTRAM, color = NES.DABTRAM)) +
  scale_color_gradient2(low = "#604CC3", mid = 'bisque', high = "#FFA500", midpoint = 0) +
  # geom_vline(xintercept = 0, linetype = 'dashed', color = 'black') +
  scale_size_continuous(range = c(1, 8), breaks = c(2, 4, 6)) +
  labs(title = 'GSEA DABTRAM Day 10', x = '', y = '') +
  theme_Publication(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


GSEA_res.day10.DABTRAM.df.use <- GSEA_res.day10.DABTRAM.df[, c('ID', 'NES.DABTRAM', 'NES.DABTRAM.abs', 'neg_log10_pval.DABTRAM')]
GSEA_res.day10.COCL2.df.use <- GSEA_res.day10.COCL2.df[, c('ID', 'NES.COCL2', 'NES.COCL2.abs', 'neg_log10_pval.COCL2')]
GSEA_res.day10.CIS.df.use <- GSEA_res.day10.CIS.df[, c('ID', 'NES.CIS', 'NES.CIS.abs', 'neg_log10_pval.CIS')]

colnames(GSEA_res.day10.DABTRAM.df.use) <- c('ID', 'NES', 'NES.abs', 'neg_log10_pval')
colnames(GSEA_res.day10.COCL2.df.use) <- c('ID', 'NES', 'NES.abs', 'neg_log10_pval')
colnames(GSEA_res.day10.CIS.df.use) <- c('ID', 'NES', 'NES.abs', 'neg_log10_pval')

GSEA_res.day10.DABTRAM.df.use$treatment <- 'DABTRAM'
GSEA_res.day10.COCL2.df.use$treatment <- 'COCL2'
GSEA_res.day10.CIS.df.use$treatment <- 'CIS'

GSEA_res <- rbind(GSEA_res.day10.DABTRAM.df.use, GSEA_res.day10.COCL2.df.use)
GSEA_res <- rbind(GSEA_res, GSEA_res.day10.CIS.df.use)
GSEA_res$NES <- ifelse(GSEA_res$NES > 2, 2, GSEA_res$NES)
ggplot(GSEA_res, aes(x = treatment, y = reorder(ID, NES), color = NES, size = neg_log10_pval)) +
  geom_point() +
  scale_color_gradient2(low = "#604CC3", mid = 'bisque', high = "#FFA500", midpoint = 0) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'black') +
  scale_size_continuous(range = c(1, 8), breaks = c(2, 4, 6)) +
  labs(title = 'GSEA Day 10', x = '', y = '') +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
# ==============================================================================
# Compare 
# ==============================================================================
saver_cor_vec.all <- merge(saver_cor_vec.day10_DABTRAM, saver_cor_vec.day10_COCL2, by = 'gene')
saver_cor_vec.all <- merge(saver_cor_vec.all, saver_cor_vec.day10_CIS, by = 'gene')

ggpairs(saver_cor_vec.all[, c('correlation', 'correlation.x', 'correlation.y')])
