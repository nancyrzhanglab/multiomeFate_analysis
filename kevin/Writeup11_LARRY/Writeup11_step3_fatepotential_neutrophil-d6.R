rm(list=ls())
library(Seurat)
library(multiomeFate)

out_folder <- "~/project/Multiome_fate/out/kevin/Writeup11/"
load(paste0(out_folder, "Writeup11_larry_seurat_fasttopics.RData"))

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

###############
# we'll focus on predict all the day6 cells at multiple timepoints from day4

treatment_vec <- as.character(sort(unique(seurat_obj$time_celltype)))
day_early_vec <- treatment_vec[grep("^.*-4", treatment_vec)]
treatment_vec <- treatment_vec[grep("^.*-6", treatment_vec)]

treatment <- treatment_vec[2]
day_early <- "4"
day_later <- treatment
seurat_obj2 <- seurat_obj

# keep only the relevant cells
keep_vec <- rep(FALSE, ncol(seurat_obj2))
idx <- which(seurat_obj2$time_celltype %in% c(day_early_vec,treatment))
keep_vec[idx] <- TRUE
seurat_obj2$keep <- keep_vec
seurat_obj2 <- subset(seurat_obj2, keep == TRUE)

tab_mat <- table(seurat_obj2$assigned_lineage, droplevels(seurat_obj2$time_celltype))

# keep only the relevant cells for this analysis
keep_vec <- rep(FALSE, ncol(seurat_obj2))
keep_vec[which(seurat_obj2$time_celltype %in% day_early_vec)] <- TRUE
seurat_obj2$keep <- keep_vec
seurat_obj2 <- subset(seurat_obj2, keep == TRUE)

# construct cell_features matrix
cell_features <- seurat_obj2[["fasttopic"]]@cell.embeddings
cell_features <- scale(cell_features)
p <- ncol(cell_features)

cell_lineage <- seurat_obj2$assigned_lineage
uniq_lineage <- sort(unique(cell_lineage))
lineage_future_count <- tab_mat[uniq_lineage, day_later]
lineage_current_count <- rowSums(tab_mat[uniq_lineage, day_early_vec])
tab_mat <- tab_mat[uniq_lineage,]

lambda_initial <- 3

fit_res <- multiomeFate:::lineage_cv(
  cell_features = cell_features,
  cell_lineage = cell_lineage,
  future_timepoint = treatment,
  lineage_future_count = lineage_future_count,
  lambda_initial = lambda_initial,
  lambda_sequence_length = 20,
  tab_mat = tab_mat,
  num_folds = 20,
  savefile_tmp = paste0(out_folder, "Writeup11_fatepotential_for_", day_later, "_tmp.RData"),
  verbose = 4
)

save(date_of_run, session_info,
     fit_res,
     file = paste0(out_folder, "Writeup11_fatepotential_for_", day_later, ".RData"))

final_fit <- multiomeFate:::lineage_cv_finalize(
  cell_features = cell_features,
  cell_lineage = cell_lineage,
  fit_res = fit_res,
  lineage_future_count = lineage_future_count
)
lineage_imputed_count <- final_fit$lineage_imputed_count
cell_imputed_score <- final_fit$cell_imputed_score

save(date_of_run, session_info,
     cell_features,
     cell_imputed_score,
     cell_lineage,
     final_fit,
     fit_res,
     lineage_current_count,
     lineage_imputed_count,
     lineage_future_count,
     tab_mat,
     treatment,
     file = paste0(out_folder, "Writeup11_fatepotential_for_", day_later, ".RData"))

print("Done! :)")
