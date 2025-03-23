rm(list = ls())

set.seed(123)

library(GSEABase) #functions and methods for Gene Set Enrichment Analysis
library(Biobase) #base functions for bioconductor; required by GSEABase
library(clusterProfiler) # provides a suite of tools for functional enrichment analysis
library(tidyverse)
library(ggpubr)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
results_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task0_explore_lineage_variability_V2/'
figure_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/figures/MultiomeFate_V3/Fig4/'

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

lin_var.day10_DABTRAM <- read.csv(paste0(results_dir, 'day10_DABTRAM/lineage_variability_day10_DABTRAM_saver_sample.csv'))
lin_var.day10_COCL2 <- read.csv(paste0(results_dir, 'day10_COCL2/lineage_variability_day10_COCL2_saver_sample.csv'))
lin_var.day10_CIS <- read.csv(paste0(results_dir, 'day10_CIS/lineage_variability_day10_CIS_saver_sample.csv'))

# =============================================================================
# calcualte mean gene expression
# =============================================================================
metadat <- all_data@meta.data
metadat$cell_id <- rownames(metadat)

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


df.day10_DABTRAM <- merge(mean.rna.day10_DABTRAM, lin_var.day10_DABTRAM, by = 'assigned_lineage') %>% drop_na()

# COCL2
metadat.day10_COCL2 <- metadat[metadat$dataset == 'day10_COCL2', ]
mean.rna.day10_COCL2 <- mat.rna %>% 
  filter(cell_id %in% metadat.day10_COCL2$cell_id) %>%
  group_by(assigned_lineage) %>%
  select(-cell_id) %>%
  summarise_all(mean, na.rm = TRUE)
df.day10_COCL2 <- merge(mean.rna.day10_COCL2, lin_var.day10_COCL2, by = 'assigned_lineage') %>% drop_na()

# CIS
metadat.day10_CIS <- metadat[metadat$dataset == 'day10_CIS', ]
mean.rna.day10_CIS <- mat.rna %>% 
  filter(cell_id %in% metadat.day10_CIS$cell_id) %>%
  group_by(assigned_lineage) %>%
  select(-cell_id) %>%
  summarise_all(mean, na.rm = TRUE)
df.day10_CIS <- merge(mean.rna.day10_CIS, lin_var.day10_CIS, by = 'assigned_lineage') %>% drop_na()

# ==============================================================================
# correlate variance and gene expression means
# ==============================================================================
genes <- rownames(all_data@assays[["Saver"]]@data)

# DABTRAM
linvar_cor_vec.day10_DABTRAM <- sapply(genes, function(j){
  if(quantile(df.day10_DABTRAM[[j]], na.rm = T)[4] == 10) {
    return(c(NA, NA))
  }
  res <- stats::cor.test(df.day10_DABTRAM[['normalized_avg_eud_dist_by_shuffle']], df.day10_DABTRAM[,j],
                         method = "spearman")
  c(res$estimate, res$p.value)
})

linvar_cor_vec.day10_DABTRAM <- as.data.frame(t(linvar_cor_vec.day10_DABTRAM))
colnames(linvar_cor_vec.day10_DABTRAM) <- c("correlation", "p.value")
rownames(linvar_cor_vec.day10_DABTRAM) <- genes

# COCL2
linvar_cor_vec.day10_COCL2 <- sapply(genes, function(j){
  if(quantile(df.day10_COCL2[[j]], na.rm = T)[4] == 10) {
    return(c(NA, NA))
  }
  res <- stats::cor.test(df.day10_COCL2[['normalized_avg_eud_dist_by_shuffle']], df.day10_COCL2[,j],
                         method = "spearman")
  c(res$estimate, res$p.value)
})
linvar_cor_vec.day10_COCL2 <- as.data.frame(t(linvar_cor_vec.day10_COCL2))
colnames(linvar_cor_vec.day10_COCL2) <- c("correlation", "p.value")
rownames(linvar_cor_vec.day10_COCL2) <- genes

# CIS
linvar_cor_vec.day10_CIS <- sapply(genes, function(j){
  if(quantile(df.day10_CIS[[j]], na.rm = T)[4] == 10) {
    return(c(NA, NA))
  }
  res <- stats::cor.test(df.day10_CIS[['normalized_avg_eud_dist_by_shuffle']], df.day10_CIS[,j],
                         method = "spearman")
  c(res$estimate, res$p.value)
})
linvar_cor_vec.day10_CIS <- as.data.frame(t(linvar_cor_vec.day10_CIS))
colnames(linvar_cor_vec.day10_CIS) <- c("correlation", "p.value")
rownames(linvar_cor_vec.day10_CIS) <- genes

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

linvar_cor_vec.day10_DABTRAM <- getAndSortTable(linvar_cor_vec.day10_DABTRAM)
linvar_cor_vec.day10_COCL2 <- getAndSortTable(linvar_cor_vec.day10_COCL2)
linvar_cor_vec.day10_CIS <- getAndSortTable(linvar_cor_vec.day10_CIS)

# GSEA DABTRAM
gsea_input.day10.DABTRAM <- linvar_cor_vec.day10_DABTRAM$correlation
names(gsea_input.day10.DABTRAM) <- linvar_cor_vec.day10_DABTRAM$gene

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
gsea_input.day10.COCL2 <- linvar_cor_vec.day10_COCL2$correlation
names(gsea_input.day10.COCL2) <- linvar_cor_vec.day10_COCL2$gene

GSEA_res.day10.COCL2 <- GSEA(geneList = gsea_input.day10.COCL2, 
                             TERM2GENE = hallmark, 
                             pvalueCutoff = 0.2,
                             seed = T,
                             verbose = F)
GSEA_res.day10.COCL2.df <- as_tibble(GSEA_res.day10.COCL2@result)
GSEA_res.day10.COCL2.df <- GSEA_res.day10.COCL2.df[, c('ID', 'NES', 'p.adjust', 'qvalue', 'setSize')]
colnames(GSEA_res.day10.COCL2.df) <- c('ID', 'NES.COCL2', 'p.adjust.COCL2', 'qvalue.COCL2', 'setSize')

# GSEA CIS
gsea_input.day10.CIS <- linvar_cor_vec.day10_CIS$correlation
names(gsea_input.day10.CIS) <- linvar_cor_vec.day10_CIS$gene

GSEA_res.day10.CIS <- GSEA(geneList = gsea_input.day10.CIS, 
                           TERM2GENE = hallmark, 
                           pvalueCutoff = 0.2,
                           seed = T,
                           verbose = F)
GSEA_res.day10.CIS.df <- as_tibble(GSEA_res.day10.CIS@result)
GSEA_res.day10.CIS.df <- GSEA_res.day10.CIS.df[, c('ID', 'NES', 'p.adjust', 'qvalue', 'setSize')]
colnames(GSEA_res.day10.CIS.df) <- c('ID', 'NES.CIS', 'p.adjust.CIS', 'qvalue.CIS', 'setSize')

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
GSEA_res$NES <- ifelse(GSEA_res$NES < -2, -2, GSEA_res$NES)

GSEA_res$ID.2 <- gsub('HALLMARK_', '', GSEA_res$ID)
GSEA_res$ID.2 <- gsub('_', ' ', GSEA_res$ID.2)

# only capitalize first letter
GSEA_res$ID.2 <- tolower(GSEA_res$ID.2)
GSEA_res$ID.2 <- gsub('(^|\\W)(\\w)', '\\1\\U\\2', GSEA_res$ID.2, perl = TRUE) 


p1 <- ggplot(GSEA_res, aes(x = treatment, y = reorder(ID.2, NES), color = NES, size = neg_log10_pval)) +
  geom_point() +
  scale_color_gradient2(low = "#606676", mid = '#C0C0C0', high = "#FFA725", midpoint = 0) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'black') +
  scale_size_continuous(range = c(1, 8), breaks = c(2, 4, 6)) +
  labs( x = '', y = '') +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

p1
ggsave(paste0(figure_dir, 'GSEA_res_day10.pdf'), p1, width = 6, height = 6)



