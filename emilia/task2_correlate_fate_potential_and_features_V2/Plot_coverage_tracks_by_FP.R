library(Seurat)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'

############################################################################################################
# Read data general
############################################################################################################

load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))
load(paste0(data_dir, 'Writeup10a_data_atac.RData'))

all_data@misc <- all_data_fatepotential
all_data[["ATAC"]] <- all_data_atac
all_data@active.assay <- 'ATAC'

############################################################################################################
# DAY 10 DABTRAM
############################################################################################################

raw_data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/day10_DABTRAM/Raw/'
# ==============================================================================
# Read data
# ==============================================================================
all_data@assays[["ATAC"]]@fragments[[1]]@path <- paste0(raw_data_dir, 'atac_fragments.tsv.gz')
all_data@assays[["ATAC"]]@fragments[[2]]@path <- paste0(raw_data_dir, 'atac_fragments.tsv.gz')
all_data@assays[["ATAC"]]@fragments[[3]]@path <- paste0(raw_data_dir, 'atac_fragments.tsv.gz')
all_data@assays[["ATAC"]]@fragments[[4]]@path <- paste0(raw_data_dir, 'atac_fragments.tsv.gz')
all_data@assays[["ATAC"]]@fragments[[5]]@path <- paste0(raw_data_dir, 'atac_fragments.tsv.gz')
all_data@assays[["ATAC"]]@fragments[[6]]@path <- paste0(raw_data_dir, 'atac_fragments.tsv.gz')
all_data@assays[["ATAC"]]@fragments[[7]]@path <- paste0(raw_data_dir, 'atac_fragments.tsv.gz')

# ==============================================================================
# Plot coverage
# ==============================================================================

# subset data to day 10
all_data_d10_dabtram <- subset(all_data, dataset == 'day10_DABTRAM')
all_data_d10_dabtram_meta <- all_data_d10_dabtram@meta.data
all_data_d10_dabtram_meta$category_0 <- ifelse(all_data_d10_dabtram[['fatepotential_DABTRAM_d10_w5']] > 0, 'Winner', 'Other')
all_data_d10_dabtram_meta$category_0 <- factor(all_data_d10_dabtram_meta[['category_0']], levels = c('Winner', 'Other'))

all_data_d10_dabtram <- AddMetaData(all_data_d10_dabtram, all_data_d10_dabtram_meta)
all_data_d10_dabtram[['category_0']] <- factor(all_data_d10_dabtram[['category_0']], levels = c('Winner', 'Other'))


cells_to_kepp <- all_data_d10_dabtram_meta[!is.na(all_data_d10_dabtram_meta$fatepotential_DABTRAM_d10_w5), ]

all_data_d10_dabtram <- all_data_d10_dabtram[, rownames(cells_to_kepp)]

# set identity by category_0
all_data_d10_dabtram <- SetIdent(all_data_d10_dabtram, value = 'category_0')


CoveragePlot(
  object = all_data_d10_dabtram,
  region = "ACTB",
  extend.upstream = 1000,
  extend.downstream = 1000,
  window = 500,
  peaks = F,
  annotation.order = c("Winner", 'Other')
)

############################################################################################################
# DAY 10 COCL2
############################################################################################################
raw_data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/day10_COCL2/Raw/'
# ==============================================================================
# Read data
# ==============================================================================
all_data@assays[["ATAC"]]@fragments[[1]]@path <- paste0(raw_data_dir, 'atac_fragments.tsv.gz')
all_data@assays[["ATAC"]]@fragments[[2]]@path <- paste0(raw_data_dir, 'atac_fragments.tsv.gz')
all_data@assays[["ATAC"]]@fragments[[3]]@path <- paste0(raw_data_dir, 'atac_fragments.tsv.gz')
all_data@assays[["ATAC"]]@fragments[[4]]@path <- paste0(raw_data_dir, 'atac_fragments.tsv.gz')
all_data@assays[["ATAC"]]@fragments[[5]]@path <- paste0(raw_data_dir, 'atac_fragments.tsv.gz')
all_data@assays[["ATAC"]]@fragments[[6]]@path <- paste0(raw_data_dir, 'atac_fragments.tsv.gz')
all_data@assays[["ATAC"]]@fragments[[7]]@path <- paste0(raw_data_dir, 'atac_fragments.tsv.gz')

# ==============================================================================
# Plot coverage
# ==============================================================================

# subset data to day 10
all_data_d10_cocl2 <- subset(all_data, dataset == 'day10_COCL2')
all_data_d10_cocl2_meta <- all_data_d10_cocl2@meta.data
all_data_d10_cocl2_meta$category_0 <- ifelse(all_data_d10_cocl2$`fatepotential_COCL2_d10_w5` > 0, 'Winner', 'Other')
all_data_d10_cocl2_meta$category_0 <- factor(all_data_d10_cocl2_meta$`category_0`, levels = c('Winner', 'Other'))

all_data_d10_cocl2 <- AddMetaData(all_data_d10_cocl2, all_data_d10_cocl2_meta)

cells_to_keep <- all_data_d10_cocl2_meta[!is.na(all_data_d10_cocl2_meta$fatepotential_COCL2_d10_w5), ]

all_data_d10_cocl2 <- all_data_d10_cocl2[, rownames(cells_to_keep)]

# set identity by category_0
all_data_d10_cocl2 <- SetIdent(all_data_d10_cocl2, value = 'category_0')


CoveragePlot(
  object = all_data_d10_cocl2,
  region = "TIMP3",
  extend.upstream = 1000,
  extend.downstream = 1000,
  window = 500,
  peaks = F,
  annotation.order = c("Winner", 'Other')
)

############################################################################################################
# DAY 0
############################################################################################################

raw_data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/day0/Raw/'

# ==============================================================================
# Read data
# ==============================================================================
all_data@assays[["ATAC"]]@fragments[[1]]@path <- paste0(raw_data_dir, 'atac_fragments.tsv.gz')
all_data@assays[["ATAC"]]@fragments[[2]]@path <- paste0(raw_data_dir, 'atac_fragments.tsv.gz')
all_data@assays[["ATAC"]]@fragments[[3]]@path <- paste0(raw_data_dir, 'atac_fragments.tsv.gz')
all_data@assays[["ATAC"]]@fragments[[4]]@path <- paste0(raw_data_dir, 'atac_fragments.tsv.gz')
all_data@assays[["ATAC"]]@fragments[[5]]@path <- paste0(raw_data_dir, 'atac_fragments.tsv.gz')
all_data@assays[["ATAC"]]@fragments[[6]]@path <- paste0(raw_data_dir, 'atac_fragments.tsv.gz')
all_data@assays[["ATAC"]]@fragments[[7]]@path <- paste0(raw_data_dir, 'atac_fragments.tsv.gz')

# ==============================================================================
# Plot coverage
# ==============================================================================

# subset data to day 10
all_data_d0 <- subset(all_data, dataset == 'day0')
all_data_d0_meta <- all_data_d0@meta.data
all_data_d0_meta$category_0 <- ifelse(all_data_d0_meta$`fatepotential_DABTRAM_d0_d10` > 0, 'Winner', 'Other')
all_data_d0_meta$category_0 <- factor(all_data_d0_meta$`category_0`, levels = c('Winner', 'Other'))

all_data_d0 <- AddMetaData(all_data_d0, all_data_d0_meta)

cells_to_keep <- all_data_d0_meta[!is.na(all_data_d0_meta$fatepotential_DABTRAM_d0_d10), ]

all_data_d0 <- all_data_d0[, rownames(cells_to_keep)]

# set identity by category_0
all_data_d0 <- SetIdent(all_data_d0, value = 'category_0')


CoveragePlot(
  object = all_data_d0,
  region = "FOSL1",
  extend.upstream = 1000,
  extend.downstream = 1000,
  window = 500,
  peaks = F,
  annotation.order = c("Winner", 'Other')
)
