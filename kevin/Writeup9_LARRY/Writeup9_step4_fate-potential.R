rm(list=ls())
library(Seurat)
library(multiomeFate)
load("~/project/Multiome_fate/out/kevin/Writeup8/Writeup8_larry-dataset_step3_fasttopics.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

# we'll focus on predict all the day6 cells at multiple timepoints from day4

treatment_vec <- as.character(sort(unique(seurat_object$time_celltype)))
day_early_vec <- treatment_vec[grep("^.*-4", treatment_vec)]
day_early <- "4"
treatment_vec <- treatment_vec[grep("^.*-6", treatment_vec)]

treatment <- treatment_vec[1]
day_later <- treatment
seurat_object2 <- seurat_object

# keep only the relevant cells
keep_vec <- rep(FALSE, ncol(seurat_object2))
idx <- which(seurat_object2$time_celltype %in% c(day_early_vec,treatment))
keep_vec[idx] <- TRUE
seurat_object2$keep <- keep_vec
seurat_object2 <- subset(seurat_object2, keep == TRUE)

tab_mat <- table(seurat_object2$assigned_lineage, droplevels(seurat_object2$time_celltype))

# keep only the relevant cells for this analysis
keep_vec <- rep(FALSE, ncol(seurat_object2))
keep_vec[which(seurat_object2$time_celltype %in% day_early_vec)] <- TRUE
seurat_object2$keep <- keep_vec
seurat_object2 <- subset(seurat_object2, keep == TRUE)

# construct cell_features matrix
rna_mat <- SeuratObject::LayerData(seurat_object2, 
                                   layer = "data", 
                                   assay = "RNA",
                                   features = Seurat::VariableFeatures(seurat_object2))
rna_mat <- as.matrix(t(rna_mat))
mean_vec <- colMeans(rna_mat)
sd_vec <- apply(rna_mat, 2, stats::sd)
rm_idx <- which(sd_vec <= 1e-3)
if(length(rm_idx) > 0){
  rna_mat <- rna_mat[,-rm_idx,drop = FALSE]
}
cell_features <- scale(rna_mat)
p <- ncol(cell_features)

cell_lineage <- seurat_object2$assigned_lineage
uniq_lineage <- sort(unique(cell_lineage))
lineage_future_count <- tab_mat[uniq_lineage, day_later]
lineage_current_count <- rowSums(tab_mat[uniq_lineage, day_early_vec])
tab_mat <- tab_mat[uniq_lineage,]

lambda_initial <- 10

fit_res <- multiomeFate:::lineage_cv(
  cell_features = cell_features,
  cell_lineage = cell_lineage,
  future_timepoint = treatment,
  lineage_future_count = lineage_future_count,
  lambda_initial = 30,
  lambda_sequence_length = 4,
  tab_mat = tab_mat,
  num_folds = 5,
  savefile_tmp = paste0("~/project/Multiome_fate/out/kevin/Writeup9/Writeup9_", treatment, "_from_day", day_early, "_lineage-imputation_tmp.RData"),
  verbose = 4
)


save(date_of_run, session_info,
     fit_res,
     file = paste0("~/project/Multiome_fate/out/kevin/Writeup9/Writeup9_", treatment, "_from_day", day_early, "_lineage-imputation.RData"))

final_fit <- lineage_cv_finalize(
  cell_features = cell_features,
  cell_lineage = cell_lineage,
  fit_res = fit_res,
  lineage_future_count = lineage_future_count
)

save(date_of_run, session_info,
     cell_features,
     cell_lineage,
     lineage_current_count,
     lineage_future_count,
     tab_mat,
     treatment,
     file = paste0("~/project/Multiome_fate/out/kevin/Writeup9/Writeup9_", treatment, "_from_day", day_early, "_lineage-imputation.RData"))

print("Done! :)")
