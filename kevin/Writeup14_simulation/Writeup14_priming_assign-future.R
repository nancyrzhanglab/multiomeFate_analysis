rm(list=ls())
library(Seurat)
library(multiomeFate)

out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup14/Writeup14_priming-setting_"
plot_folder <- "~/project/Multiome_fate/git/multiomeFate_analysis_kevin/fig/kevin/Writeup14/Writeup14_priming-setting_"
func_folder <- "~/project/Multiome_fate/git/multiomeFate_analysis_kevin/kevin/Writeup14_simulation/"

source(paste0(func_folder, "func.R"))
source(paste0(func_folder, "func_priming.R"))
source(paste0(func_folder, "func_assign-future-cells.R"))
source(paste0(func_folder, "func_seurat.R"))

#######

all_data <- multiomeFate:::data_loader(which_files = c( "rna", "saver", "fasttopics"))
keep_vec <- rep(FALSE, length(Seurat::Cells(all_data)))
keep_vec[all_data$dataset %in% c("day10_COCL2", "week5_COCL2")] <- TRUE
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == TRUE)
all_data[["fasttopic.CIS"]] <- NULL
all_data[["ft.CIS.umap"]] <- NULL
all_data[["fasttopic.DABTRAM"]] <- NULL
all_data[["ft.DABTRAM.umap"]] <- NULL

load(paste0(out_folder, "simulation.RData"))

# count how many future cells there's supposed to be
total_future_cells <- sum(simulation_res$lineage_future_size)
current_future_cells <- length(which(all_data$dataset == "week5_COCL2"))

total_future_cells; current_future_cells

# here, since total_future_cells < current_future_cells, we simply throw out cells
names(simulation_res$lineage_assignment) <- rownames(simulation_res$embedding_mat)
simulation_data <- .form_simulation_seurat_fate(final_fit = final_fit,
                                                simulation_res = simulation_res)

# subset only "day10_COCL2" and "week5_COCL2"
lineage_vec <- rep(NA, length(Seurat::Cells(all_data)))
names(lineage_vec) <- Seurat::Cells(all_data)
lineage_vec[Seurat::Cells(simulation_data)] <- as.character(simulation_data$assigned_lineage)
all_data$assigned_lineage <- lineage_vec

set.seed(10)
all_data <- .assign_future_cells(
  seurat_obj = all_data,
  future_ident = "week5_COCL2",
  ident = "dataset",
  lineage_future_size = simulation_res$lineage_future_size,
  lineage_variable = "assigned_lineage"
)
fatepotential_true <- rep(NA, length(Seurat::Cells(all_data)))
names(fatepotential_true) <- Seurat::Cells(all_data)
fatepotential_true[Seurat::Cells(simulation_data)] <- simulation_data$fatepotential_true
all_data$fatepotential_true <- fatepotential_true

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
sum(abs(tab_mat[names(simulation_res$lineage_future_size),"week5_COCL2"])-sort(simulation_res$lineage_future_size))

all_data
table(all_data$dataset)
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
head(tab_mat)
quantile(tab_mat[,"day10_COCL2"])
quantile(tab_mat[,"week5_COCL2"])
length(unique(all_data$assigned_lineage))
head(all_data$fatepotential_true)

date_of_run <- Sys.time()
session_info <- devtools::session_info()

save(
  all_data,
  date_of_run, session_info,
  file = paste0(out_folder, "seurat_CoSpar_prepared.RData")
)

