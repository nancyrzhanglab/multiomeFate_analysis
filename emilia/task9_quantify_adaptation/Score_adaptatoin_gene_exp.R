rm(list = ls())

set.seed(123)

library(Seurat)
library(UCell)
library(tidyverse)
library(ggpubr)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'
# figure_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/figures/MultiomeFate_V3/Fig6/'

treatment <- 'DABTRAM'

remove_unassigned_cells <- TRUE

date_of_run <- Sys.time()
session_info <- devtools::session_info()

# =============================================================================
# reading data
# =============================================================================

load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_saver.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))

all_data[['saver']] <- all_data_saver
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

# ==============================================================================
# Sigantures
# ==============================================================================

custom.dabtram <- c('ACTB', 'TMEM43', 'TPM4', 'CALM2', 'FN1', 'PALLD', 'LMO7',
                    'ACTN1', 'HSPG2', 'MYOF', 'TNFRSF12A', 'TUBB', 'RCN1', 'CRIM1',
                    'COL5A2', 'SAMD5', 'TPM1', 'OXSR1', 'CBX5')
custom.cis <- c('NUCKS1','TUBB','HMGB2', 'ICMT', 'CBX5', 'TUBA1B', 'ANP32B','TYMS',
                'GMNN', 'USP1', 'NASP', 'TMPO', 'NCAPH', 'TK1', 'TUBG1', 'PRC1',
                'PBK', 'SMC3', 'RRM2', 'RAD51AP1')
custom.cocl2 <- c('GXYLT2', 'ANTXR1', 'CADM1', 'ITGB3', 'BICC1', 'SLC1A4',
                  'CADPS', 'HMGA2', 'TIMP3', 'PTPRG', 'SERPINE2', 'IMMP2L', 'LRMDA',
                  'MFSD12', 'SOX5', 'EPHA3', 'PRKG2', 'IL1RAP', 'SLC44A1', 'KCNQ5')
jackpot = c('SOX10', 'MITF', 'FN1', 'AXL', 'EGFR', 'NT5E',
            'C1S', 'FRZB', 'SERPINB2', 'SERPINE1', 'NGFR',
            'NDRG1', "FEZF1", 'EGR3', 'VGF',
            'WNT5A', 'POSTN', 'PDGFRB', 'NRG1', 'VEGFC', 'FOSL1',
            'RUNX2', 'LOXL2', 'JUN', 'PDGFRC', 'CD44', 'ID3')

DABTRAM = c('AXL', 'EGFP', 'NGFR', 'IGFBP5', 'ANXA1',
            'IGFBP7', 'JUNB', 'BASP1', 'IER2', 'JUN',
            'CXCL12', 'ANXA2', 'FOS', 'MMP2', 'GLRX',
            'IL6ST', 'PRNP', 'FOSB', 'CTSL', 'SLC12A8',
            'TFPI2', 'MYL6', 'IFITM3', 'CAV1', 'CD44')

isg.rs <- c('IFI27', 'IRF7', 'USP18', 'BST2', 'CXCL10',
            'DDX60', 'HERC6', 'HLA-B', 'HLA-G',
            'IFI35', 'IFI44', 'IFI44L', 'IFIT1', 'IFIT3',
            'ISG15', 'LGALS3BP', 'LY6E', 'MX1', 'MX2',
            'OAS3', 'OASL', 'PLSCR1', 'STAT1', 'TRIM14',
            'HSD17B1', 'OAS1', 'CA2', 'CCNA1', 'CXCL1',
            'GALC', 'IFI6', 'IFITM1', 'LAMP3', 'MCL1',
            'ROBO1', 'SLC6A15', 'THBS1', 'TIMP3')

jackpot <- jackpot[jackpot %in% all_data@assays[["saver"]]@var.features]
DABTRAM <- DABTRAM[DABTRAM %in% all_data@assays[["saver"]]@var.features]
isg.rs <- isg.rs[isg.rs %in% all_data@assays[["saver"]]@var.features]

scores <- ScoreSignatures_UCell(all_data@assays[["saver"]]@data, features=list('dabtram.adaptation.genes' = custom.dabtram,
                                                                               'cis.adaptation.genes' = custom.cis,
                                                                               'cocl2.adaptation.genes' = custom.cocl2,
                                                                               'jackpot' = jackpot,
                                                                               'dabtram' = DABTRAM,
                                                                               'isg.rs' = isg.rs))

write.csv(scores, paste0(out_dir, 'UCell_scores_adaptation_gene_sets.csv'), row.names = TRUE)
