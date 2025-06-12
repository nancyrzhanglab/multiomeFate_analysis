embedding_mat = cell_features
bool_replace = TRUE
num_rounds = 20
num_trials = 20
min_range_percentage = 0.5
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

###### ???

sum_vec <- sapply(lineage_names, function(lineage){
  idx <- which(lineage_assignment == lineage)
  sum(progenies_vec[idx])
})

#############

.evaluate_lineage_assignment_priming(
  lineage_assignment = lineage_assignment,
  lineage_names = lineage_names,
  min_range_percentage = min_range_percentage,
  progenies_vec = progenies_vec
)

######

tmp <- .compute_future_lineage_statistics(lineage_assignment = lineage_assignment,
                                          lineage_names = lineage_names,
                                          progenies_vec = progenies_vec)
median_cor <- tmp$median_cor
range_cor <- tmp$range_cor
size_vec <- tmp$size_vec