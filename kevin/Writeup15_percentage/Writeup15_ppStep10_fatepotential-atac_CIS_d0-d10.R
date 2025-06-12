rm(list=ls())
library(Seurat)
library(multiomeFate)

out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup15/"
all_data <- multiomeFate:::data_loader(which_files = c("peakvi"))

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

# we'll focus on predict all the day6 cells at multiple timepoints from day4

treatment <- "CIS"
day_early_vec <- "day0"
day_early <- "d0"
day_later_vec <- paste0("day10_", treatment)
day_later <- "d10"

# keep only the relevant cells
keep_vec <- rep(FALSE, ncol(all_data))
idx <- which(all_data$dataset %in% c("day0", day_early_vec, day_later_vec))
keep_vec[idx] <- TRUE
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == TRUE)

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)

# keep only the relevant cells for this analysis
keep_vec <- rep(FALSE, ncol(all_data))
keep_vec[which(all_data$dataset %in% day_early_vec)] <- TRUE
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == TRUE)

# construct cell_features matrix
atac_features <- all_data[[paste0("peakVI.", treatment)]]@cell.embeddings
atac_features <- scale(atac_features)
cell_features <- atac_features
p <- ncol(cell_features)

cell_lineage <- all_data$assigned_lineage
uniq_lineage <- sort(unique(cell_lineage))
lineage_future_count <- tab_mat[uniq_lineage, day_later_vec]
lineage_current_count <- tab_mat[uniq_lineage, day_early_vec]
tab_mat <- tab_mat[uniq_lineage,,drop = FALSE]

lambda_initial <- 3

fit_res <- multiomeFate:::lineage_cv(
  cell_features = cell_features,
  cell_lineage = cell_lineage,
  future_timepoint = day_later_vec,
  lineage_future_count = lineage_future_count,
  lambda_initial = lambda_initial,
  lambda_sequence_length = 20,
  tab_mat = tab_mat,
  num_folds = 20,
  savefile_tmp = paste0(out_folder, "Writeup15_", treatment, "-from-", day_early, "_fatepotential-atac_tmp.RData"),
  verbose = 4
)

save(date_of_run, session_info,
     fit_res,
     file = paste0(out_folder, "Writeup15_", treatment, "-from-", day_early, "_fatepotential-atac.RData"))

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
     file = paste0(out_folder, "Writeup15_", treatment, "-from-", day_early, "_fatepotential-atac.RData"))

print("Done! :)")
