embedding_mat = cell_features
bool_replace = TRUE
num_rounds = 20
num_trials = 20
percentage_randomize = 0.05
tol = 1e-5

##########

n <- nrow(embedding_mat)

progenies_vec <- .compute_mean_total_cells(
  cell_features = cell_features,
  coefficient_intercept = coefficient_intercept,
  coefficient_vec = coefficient_vec,
  return_sum = FALSE
)

# start by assigning lineages in order
lineage_names <- paste0("lineage:", 1:num_lineages)
lineage_assignment <- .initial_assignment(
  lineage_names,
  n,
  num_lineages,
  progenies_vec
)

# try reordering the lineages, which the goal of: correlation of range is non-positive, and maximize correlation of median
for(round in 1:num_rounds){
  idx <- sample(1:n, size = round(percentage_randomize * n), replace = FALSE)
  
  obj_vec <- rep(0, num_trials+1)
  assignment_list <- vector("list", length = num_trials+1)
  
  for(trial in 1:num_trials){
    attempted_lineage_assignment <- lineage_assignment
    attempted_lineage_assignment[idx] <- sample(attempted_lineage_assignment[idx], replace = bool_replace)
    
    obj_vec[trial] <- .evaluate_lineage_assignment_priming(
      lineage_assignment = attempted_lineage_assignment,
      lineage_names = lineage_names,
      progenies_vec = progenies_vec
    )
    assignment_list[[trial]] <- attempted_lineage_assignment[idx]
  }
  # put the current trial in
  obj_vec[num_trials+1] <- .evaluate_lineage_assignment_priming(
    lineage_assignment = lineage_assignment,
    lineage_names = lineage_names,
    progenies_vec = progenies_vec
  )
  assignment_list[[num_trials+1]] <- lineage_assignment[idx]
  
  # pick the best
  if(any(!is.na(obj_vec))){
    trial_select <- which.max(obj_vec)
    lineage_assignment[idx] <- assignment_list[[trial_select]]
  }
}


# finalize
lineage_assignment <- factor(lineage_assignment)
lineage_future_size <- sapply(levels(lineage_assignment), function(lev) {
  idx <- which(lineage_assignment == lev)
  round(sum(progenies_vec[idx]))
})
names(lineage_future_size) <- levels(lineage_assignment)

summary_mat <- multiomeFate:::.compute_summary_lineages(cell_fate_potential_truth = log10(progenies_vec),
                                                        lineage_assignment = lineage_assignment,
                                                        lineage_future_size = lineage_future_size)
summary_mat

#######

tmp <- .compute_future_lineage_statistics(lineage_assignment = lineage_assignment,
                                          lineage_names = lineage_names,
                                          progenies_vec = progenies_vec)

#######

lineage_assignment_list <- lapply(lineage_names, function(lineage){
  which(lineage_assignment == lineage)
})

size_vec <- .compute_future_lineage_size(lineage_assignment_list, progenies_vec)
median_vec <- .compute_future_lineage_median(lineage_assignment_list, progenies_vec)
range_vec <- .compute_future_lineage_range(lineage_assignment_list, progenies_vec)

range_cor <- stats::cor(size_vec, range_vec)
median_cor <- stats::cor(size_vec, median_vec)
