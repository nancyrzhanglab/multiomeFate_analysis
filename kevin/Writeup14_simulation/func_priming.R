.assign_lineages_priming <- function(embedding_mat, 
                                     bool_replace = TRUE,
                                     coefficient_intercept = 0, 
                                     coefficient_vec = rep(1, ncol(embedding_mat)),
                                     num_lineages = 50,
                                     min_range_percentage = 0.5,
                                     num_rounds = 20,
                                     num_trials = 20,
                                     percentage_randomize = 0.05,
                                     verbose = 0){
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
    if(verbose > 0) print(paste0("On round: ", round))
    idx <- sample(1:n, size = round(percentage_randomize * n), replace = FALSE)
    
    obj_vec <- rep(0, num_trials+1)
    assignment_list <- vector("list", length = num_trials+1)
    
    for(trial in 1:num_trials){
      attempted_lineage_assignment <- lineage_assignment
      attempted_lineage_assignment[idx] <- sample(attempted_lineage_assignment[idx], replace = bool_replace)
      
      obj_vec[trial] <- .evaluate_lineage_assignment_priming(
        lineage_assignment = attempted_lineage_assignment,
        lineage_names = lineage_names,
        min_range_percentage = min_range_percentage,
        progenies_vec = progenies_vec
      )
      assignment_list[[trial]] <- attempted_lineage_assignment[idx]
    }
    # put the current trial in
    obj_vec[num_trials+1] <- .evaluate_lineage_assignment_priming(
      lineage_assignment = lineage_assignment,
      lineage_names = lineage_names,
      min_range_percentage = min_range_percentage,
      progenies_vec = progenies_vec
    )
    assignment_list[[num_trials+1]] <- lineage_assignment[idx]
    
    # pick the best
    if(any(!is.na(obj_vec))){
      trial_select <- which.max(obj_vec)
      lineage_assignment[idx] <- assignment_list[[trial_select]]
    } else {
      if(verbose) print("Failed to update")
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
  
  # return
  list(coefficient_intercept = coefficient_intercept,
       coefficient_vec = coefficient_vec,
       embedding_mat = embedding_mat,
       lineage_assignment = lineage_assignment,
       lineage_future_size = lineage_future_size,
       summary_mat = summary_mat)
}

##########################

.initial_assignment <- function(lineage_names,
                                n,
                                num_lineages,
                                progenies_vec){
  tmp <- rep(lineage_names, each = ceiling(n/num_lineages))[1:n]
  tmp[rank(progenies_vec)]
}

.evaluate_lineage_assignment_priming <- function(lineage_assignment,
                                                 lineage_names,
                                                 min_range_percentage,
                                                 progenies_vec){
  tmp <- .compute_future_lineage_statistics(lineage_assignment = lineage_assignment,
                                            lineage_names = lineage_names,
                                            progenies_vec = progenies_vec)
  median_cor <- tmp$median_cor
  range_cor <- tmp$range_cor
  size_vec <- tmp$size_vec
  
  if(diff(range(size_vec)) < min_range_percentage*max(size_vec)) return(NA)
  
  if(range_cor > 0) return(NA)
  return(median_cor)
}

########

.compute_future_lineage_statistics <- function(lineage_assignment,
                                               lineage_names,
                                               progenies_vec){
  lineage_assignment_list <- lapply(lineage_names, function(lineage){
    which(lineage_assignment == lineage)
  })
  
  size_vec <- .compute_future_lineage_size(lineage_assignment_list, progenies_vec)
  median_vec <- .compute_future_lineage_median(lineage_assignment_list, progenies_vec)
  range_vec <- .compute_future_lineage_range(lineage_assignment_list, progenies_vec)
  
  range_cor <- stats::cor(size_vec, range_vec)
  median_cor <- stats::cor(size_vec, median_vec)
  
  return(list(median_cor = median_cor,
              range_cor = range_cor,
              size_vec = size_vec))
}

.compute_future_lineage_size <- function(lineage_assignment_list,
                                         progenies_vec){
  sapply(1:length(lineage_assignment_list), function(i){
    idx <- lineage_assignment_list[[i]]
    if(length(idx) == 0) return(0)
    return(sum(progenies_vec[idx]))
  })
}

.compute_future_lineage_median <- function(lineage_assignment_list,
                                           progenies_vec){
  sapply(1:length(lineage_assignment_list), function(i){
    idx <- lineage_assignment_list[[i]]
    if(length(idx) == 0) return(0)
    return(stats::median(progenies_vec[idx]))
  })
}

.compute_future_lineage_range <- function(lineage_assignment_list,
                                          progenies_vec){
  sapply(1:length(lineage_assignment_list), function(i){
    idx <- lineage_assignment_list[[i]]
    if(length(idx) == 0) return(0)
    return(diff(range(progenies_vec[idx])))
  })
}