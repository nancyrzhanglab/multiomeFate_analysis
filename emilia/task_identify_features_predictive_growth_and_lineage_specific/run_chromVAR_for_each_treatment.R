library(Seurat)
library(Signac)
library(BSgenome.Hsapiens.UCSC.hg38)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

cis <- readRDS('/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/emiliac/lineage_trace/task_identify_features_predictive_growth_and_lineage_specific/data/week5_CIS.rds')
cocl2 <- readRDS('/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/emiliac/lineage_trace/task_identify_features_predictive_growth_and_lineage_specific/data/week5_COCL2.rds')
dabtram <- readRDS('/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/emiliac/lineage_trace/task_identify_features_predictive_growth_and_lineage_specific/data/week5_DABTRAM.rds')

Seurat::DefaultAssay(cis) <- "ATAC"
Seurat::DefaultAssay(cocl2) <- "ATAC"
Seurat::DefaultAssay(dabtram) <- "ATAC"

cis <- Signac::RunChromVAR(
  object = cis,
  genome = "hg38"
)

cocl2 <- Signac::RunChromVAR(
  object = cocl2,
  genome = "hg38"
)

dabtram <- Signac::RunChromVAR(
  object = dabtram,
  genome = "hg38"
)

chromvar_results_cis <- cis@assays[["chromvar"]]@data
chromvar_results_cocl2 <- cocl2@assays[["chromvar"]]@data
chromvar_results_dabtram <- dabtram@assays[["chromvar"]]@data

save(date_of_run, session_info, 
     cis, cocl2, dabtram,
     file = "/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/emiliac/lineage_trace/task_identify_features_predictive_growth_and_lineage_specific/data/chromVar_week5_data.RData")