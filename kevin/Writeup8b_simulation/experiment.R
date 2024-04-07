rm(list=ls())
library(Seurat)
library(multiomeFate)

load("../../out/Writeup8b/Writeup8b_simulation_day10-COCL2.RData")

par(mfrow = c(1,1))
Seurat::DimPlot(all_data,
                group.by = "dataset",
                reduction = "umap")

embedding_mat <- all_data[["fasttopic_COCL2"]]@cell.embeddings
embedding_mat <- scale(embedding_mat)
embedding_mat <- pmin(embedding_mat, 2)
coefficient_intercept <- -1
coefficient_vec <- rep(1, ncol(embedding_mat))
coefficient_vec <- coefficient_vec/2

# double-check the fate potentials are ok
tmp <- exp((embedding_mat %*% coefficient_vec) + coefficient_intercept)
quantile(tmp)
sum(tmp)

lineage_concentration <- 1
num_lineages <- 50
lineage_prior <- rep(1/num_lineages, length = num_lineages)

############################
# not much to change after this line

early_idx <- which(all_data$dataset == "day10_COCL2")
set.seed(10)
simulation_res <- multiomeFate:::generate_simulation(
  embedding_mat = embedding_mat[early_idx,,drop=FALSE],
  coefficient_intercept = coefficient_intercept,
  coefficient_vec = coefficient_vec,
  lineage_concentration = lineage_concentration,
  lineage_prior = lineage_prior,
  num_lineages = num_lineages,
  verbose = 3
)

# check the simulation to make the sizes look alright
table(simulation_res$lineage_assignment)
hist(simulation_res$cell_fate_potential)
hist(10^(simulation_res$cell_fate_potential))
hist(simulation_res$lineage_future_size)
sum(simulation_res$lineage_future_size)
sum(10^(simulation_res$cell_fate_potential))


par(mfrow = c(1,1))
idx <- which(simulation_res$lineage_assignment == "lineage:1")
plot(all_data[["umap"]]@cell.embeddings,
     pch = 16, col = "gray",
     xlab = "UMAP1", ylab = "UMAP2")
for(i in idx){
  points(all_data[["umap"]]@cell.embeddings[early_idx[idx],,drop = FALSE],
         pch = 16, col = "coral3", cex = 1)
}

##################

# try out our fate potential method
cell_lineage <- as.character(simulation_res$lineage_assignment)
uniq_lineage <- sort(unique(cell_lineage))
lineage_future_count <- simulation_res$lineage_future_size
tmp <- table(simulation_res$lineage_assignment)
lineage_current_count <- as.numeric(tmp); names(lineage_current_count) <- names(tmp)
lineage_current_count <- lineage_current_count[names(lineage_future_count)]
tab_mat <- cbind(lineage_current_count, lineage_future_count)
colnames(tab_mat) <- c("now", "future")

#################
# start cross validation
cell_features = embedding_mat[early_idx,,drop=FALSE]
future_timepoint = "future"
lambda_initial = 3
lambda_sequence_length = 10
num_folds = 10
verbose = 10

####################

tmp <- construct_folds(
  cell_lineage = cell_lineage,
  tab_mat = tab_mat,
  future_timepoint = future_timepoint,
  num_folds = num_folds
)
cv_cell_list <- tmp$cv_cell_list
fold_lineage_list <- tmp$fold_lineage_list

cv_fit_list <- vector("list", length = num_folds)
names(cv_fit_list) <- names(cv_cell_list)

###

i = 1
fold <- names(cv_cell_list)[i]
if(verbose > 0) print(paste0("Dropping fold #", i, " out of ", num_folds))
cell_features_train <- cell_features[-cv_cell_list[[fold]],,drop = F]
cell_lineage_train <- cell_lineage[-cv_cell_list[[fold]]]
lineage_future_count_train <- lineage_future_count[-which(names(lineage_future_count) %in% fold_lineage_list[[fold]])]

set.seed(10)
cell_features = cell_features_train
cell_lineage = cell_lineage_train
lineage_future_count = lineage_future_count_train

###

lambda_initial = NA
lambda_max = 101
lambda_min = 101
multipler = 1e4

res <- .compute_initial_parameters(cell_features = cell_features,
                                   cell_lineage = cell_lineage,
                                   lineage_future_count = lineage_future_count,
                                   lambda_max = lambda_max,
                                   lambda_min = lambda_min,
                                   multipler = multipler)
coefficient_initial <- res$coefficient_initial
if(is.na(lambda_initial)) lambda_initial <- res$lambda_initial

lambda_sequence <- exp(seq(log(lambda_initial+1), 0, 
                           length.out = lambda_sequence_length))-1
fit_list <- vector("list", length = lambda_sequence_length)
i = 1
coefficient_vec <- coefficient_initial
coefficient_initial_list = coefficient_vec
lambda = lambda_sequence[i]
random_initializations = 10
upper_randomness = 5

###

if(!is.list(coefficient_initial_list)) coefficient_initial_list <- list(coefficient_initial_list)
list_len <- length(coefficient_initial_list)

# some cleanup
tmp <- .lineage_cleanup(cell_features = cell_features,
                        cell_lineage = cell_lineage,
                        lineage_future_count = lineage_future_count,
                        verbose = verbose)
cell_features <- tmp$cell_features
cell_lineage <- tmp$cell_lineage
cell_lineage_idx_list <- tmp$cell_lineage_idx_list
lineage_future_count <- tmp$lineage_future_count
uniq_lineages <- tmp$uniq_lineages
coefficient_initial_list <- .append_intercept_term(coefficient_initial_list)
p <- ncol(cell_features)

for(i in 1:list_len){
  if(length(names(coefficient_initial_list[[i]])) != 0){
    stopifnot(all(names(coefficient_initial_list[[i]]) == colnames(cell_features)))
  } else {
    names(coefficient_initial_list[[i]]) <- colnames(cell_features)
  }
}

# rearrange arguments
optim_fn <- function(coefficient_vec,
                     cell_features,
                     cell_lineage,
                     cell_lineage_idx_list,
                     lambda,
                     lineage_future_count){
  .lineage_objective(cell_features = cell_features,
                     cell_lineage = cell_lineage,
                     cell_lineage_idx_list = cell_lineage_idx_list,
                     coefficient_vec = coefficient_vec,
                     lambda = lambda,
                     lineage_future_count = lineage_future_count)
}

optim_gr <- function(coefficient_vec,
                     cell_features,
                     cell_lineage,
                     cell_lineage_idx_list,
                     lambda,
                     lineage_future_count){
  .lineage_gradient(cell_features = cell_features,
                    cell_lineage = cell_lineage,
                    cell_lineage_idx_list = cell_lineage_idx_list,
                    coefficient_vec = coefficient_vec,
                    lambda = lambda,
                    lineage_future_count = lineage_future_count)
}

res_list <- vector("list", length = list_len+random_initializations)

max_feature <- stats::quantile(abs(cell_features), probs = 0.95)
num_cells_per_lineage <- sapply(cell_lineage_idx_list, length)
names(num_cells_per_lineage) <- uniq_lineages
max_count_ratio <- max(lineage_future_count[uniq_lineages]/num_cells_per_lineage[uniq_lineages])
max_limit <- 2*log(max_count_ratio)/(p*max_feature)
min_value <- ifelse(max_limit < 0, 2*max_limit, 0)

coef_vec <- pmin(stats::runif(p, min = min_value, max = max_limit), upper_randomness)
names(coef_vec) <- colnames(cell_features)


