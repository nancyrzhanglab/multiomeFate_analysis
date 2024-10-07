.preprocess_rna <- function(rna_mat, treatment){
  cell_idx <- which(all_data$dataset == treatment)
  rna_mat <- rna_mat[cell_idx,]
  rna_mat <- pmin(rna_mat, 10)
  
  mean_vec <- colMeans(rna_mat)
  sd_vec <- apply(rna_mat, 2, stats::sd)
  rm_idx <- which(sd_vec <= 1e-3)
  if(length(rm_idx) > 0){
    rna_mat <- rna_mat[,-rm_idx,drop = FALSE]
  }
  cell_features <- scale(rna_mat)
  cell_features
}

.compute_mean_total_cells <- function(cell_features,
                                      coefficient_intercept,
                                      coefficient_vec,
                                      return_sum){
  tmp <- exp((cell_features %*% coefficient_vec) + coefficient_intercept)[,1]
  
  if(return_sum){
    return(sum(tmp))
  } else {
    return(tmp)
  }
}

.search_for_priming_parameters <- function(cell_features,
                                           min_cells = 6000,
                                           min_maximum = 10,
                                           num_random = 20,
                                           num_fixed = 20,
                                           num_trials = 40,
                                           range_vec = c(-2, 2),
                                           verbose = 0){
  p <- ncol(cell_features)
  
  current_intercept <- 0.1
  current_vec <- rep(0.1, p)
  vec_fixed <- seq(range_vec[1], range_vec[2], length.out = num_fixed)
  
  for(trial in 1:num_trials){
    if(verbose > 0) print(paste0("Working on trial: ", trial))
    
    # first the intercept
    vec_try <- .generate_tries(num_random, range_vec, vec_fixed)
    vec_try <- c(vec_try, current_intercept)
    current_intercept <- .evalute_even_spread_vec(current_intercept = current_intercept,
                                                  current_vec = current_vec,
                                                  j = 0,
                                                  min_cells = min_cells,
                                                  min_maximum = min_maximum,
                                                  vec_try = vec_try)
    
    # now for the coefficients
    for(j in 1:p){
      vec_try <- .generate_tries(num_random, range_vec, vec_fixed)
      vec_try <- c(vec_try, current_vec[j])
      current_vec[j] <- .evalute_even_spread_vec(current_intercept = current_intercept,
                                                 current_vec = current_vec,
                                                 j = j,
                                                 min_cells = min_cells,
                                                 min_maximum = min_maximum,
                                                 vec_try = vec_try)
    }
  }
  
  return(list(coefficient_intercept = current_intercept,
              coefficient_vec = current_vec))
}

.generate_tries <- function(num_random,
                            range_vec,
                            vec_fixed){
  c(vec_fixed, stats::runif(num_random, min = range_vec[1], max = range_vec[2]))
}

.evalute_even_spread_vec <- function(current_intercept,
                                     current_vec,
                                     j,
                                     min_cells,
                                     min_maximum,
                                     vec_try){
  obj_vec <- sapply(vec_try, function(x){
    coefficient_vec <- current_vec
    if(j == 0){
      coefficient_intercept <- x
    } else {
      coefficient_intercept <- current_intercept
      coefficient_vec[j] <- x
    }
    
    num_cells <- .compute_mean_total_cells(cell_features,
                                           coefficient_intercept,
                                           coefficient_vec,
                                           return_sum = FALSE)
    if(sum(num_cells) < min_cells) return(NA)
    if(max(num_cells) < min_maximum) return(NA)
    .evaluate_even_spread(num_bins = 4,
                          vec = num_cells)
  })
  
  stopifnot(!all(is.na(obj_vec)))
  
  vec_try[which.min(obj_vec)]
}

.evaluate_even_spread <- function(num_bins,
                                  vec){
  if(diff(range(vec)) <= 1e-5) return(NA)
  
  break_vec <- seq(0, max(vec), length.out = num_bins)
  tab_vec <- table(cut(vec, breaks = break_vec))
  diff(range(tab_vec))
}

##################

.plot_mean_variance <- function(filename,
                                simulation_res){
  
  median_vec <- simulation_res$summary_mat["median",]
  range_vec <- simulation_res$summary_mat["range",]
  future_vec <- simulation_res$summary_mat["future_size",]
  
  grDevices::png(filename = filename,
                 width = 15, 
                 height = 5, 
                 units = "in", 
                 res = 300)
  par(mfrow = c(1,3))
  
  graphics::plot(median_vec, 
                 future_vec, 
                 pch = 16,
                 xlab = "Median (per lineage)",
                 ylab = "Future size (per lineage)",
                 main = paste0("Future size vs. current median potential\n",
                               "Corr: ", round(cor(median_vec, future_vec), 2)))
  graphics::plot(range_vec,
                 future_vec,
                 pch = 16,
                 xlab = "Range (per lineage)",
                 ylab = "Future size (per lineage)",
                 main = paste0("Future size vs. current range potential\n",
                               "Corr: ", round(cor(range_vec, future_vec), 2)))
  
  graphics::plot(x = median_vec,
                 y = range_vec, 
                 main = paste0("Corr: ", round(cor(range_vec, median_vec), 2)),
                 xlab = "Median (per lineage)", 
                 ylab = "Range (per lineage)",
                 pch = 16, 
                 asp = TRUE)
  
  grDevices::graphics.off()
}