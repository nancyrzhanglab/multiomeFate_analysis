rm(list=ls())

# priming
load("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/out/Writeup14/Writeup14_priming-setting_simulation_v2.RData")

# names(simulation_res)
# dim(simulation_res$embedding_mat)

set.seed(10)
tmp <- simulation_res$summary_mat["future_size",]
lineage_ordering <- names(sort(tmp, decreasing = TRUE))

# have only 15 lineages
lineage_keep <- lineage_ordering[round(seq(1, length(lineage_ordering), length.out = 15))]
cell_idx <- which(simulation_res$lineage_assignment %in% lineage_keep)

simulation_res$summary_mat <- NULL
simulation_res$embedding_mat <- simulation_res$embedding_mat[cell_idx,]
simulation_res$lineage_assignment <- simulation_res$lineage_assignment[cell_idx]
simulation_res$lineage_assignment <- droplevels(simulation_res$lineage_assignment)
simulation_res$lineage_future_size <- simulation_res$lineage_future_size[lineage_keep]
rownames(simulation_res$embedding_mat) <- paste0("cell:", 1:nrow(simulation_res$embedding_mat))

cell_features <- simulation_res$embedding_mat
cell_lineage <- as.character(simulation_res$lineage_assignment)
uniq_lineage <- sort(unique(cell_lineage))
lineage_future_count <- simulation_res$lineage_future_size
tmp <- table(simulation_res$lineage_assignment)
lineage_current_count <- as.numeric(tmp); names(lineage_current_count) <- names(tmp)
lineage_current_count <- lineage_current_count[names(lineage_future_count)]
tab_mat <- cbind(lineage_current_count, lineage_future_count)
colnames(tab_mat) <- c("now", "future")

priming_simulation <- list(
  cell_features = cell_features,
  cell_lineage = cell_lineage,
  lineage_future_count = lineage_future_count,
  tab_mat = tab_mat
)

usethis::use_data(priming_simulation, overwrite = TRUE)

##########################
rm(list=ls())

load("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/out/Writeup14/Writeup14_plastic-setting_simulation.RData")

set.seed(10)
tmp <- simulation_res$summary_mat["future_size",]
lineage_ordering <- names(sort(tmp, decreasing = TRUE))

# have only 15 lineages
lineage_keep <- lineage_ordering[round(seq(1, length(lineage_ordering), length.out = 15))]
cell_idx <- which(simulation_res$lineage_assignment %in% lineage_keep)

simulation_res$summary_mat <- NULL
simulation_res$embedding_mat <- simulation_res$embedding_mat[cell_idx,]
simulation_res$lineage_assignment <- simulation_res$lineage_assignment[cell_idx]
simulation_res$lineage_assignment <- droplevels(simulation_res$lineage_assignment)
simulation_res$lineage_future_size <- simulation_res$lineage_future_size[lineage_keep]
rownames(simulation_res$embedding_mat) <- paste0("cell:", 1:nrow(simulation_res$embedding_mat))

cell_features <- simulation_res$embedding_mat
cell_lineage <- as.character(simulation_res$lineage_assignment)
uniq_lineage <- sort(unique(cell_lineage))
lineage_future_count <- simulation_res$lineage_future_size
tmp <- table(simulation_res$lineage_assignment)
lineage_current_count <- as.numeric(tmp); names(lineage_current_count) <- names(tmp)
lineage_current_count <- lineage_current_count[names(lineage_future_count)]
tab_mat <- cbind(lineage_current_count, lineage_future_count)
colnames(tab_mat) <- c("now", "future")

plastic_simulation <- list(
  cell_features = cell_features,
  cell_lineage = cell_lineage,
  lineage_future_count = lineage_future_count,
  tab_mat = tab_mat
)

usethis::use_data(plastic_simulation, overwrite = TRUE)


