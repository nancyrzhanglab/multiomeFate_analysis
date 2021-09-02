chromatin_potential_prepare2 <- function(mat_x, mat_y, 
                                         snn,
                                         diffusion_dist,
                                         vec_start, list_end,
                                         ht_map){
  # [[NOTE: mat_x and mat_y here are the expressions of the metacells here]]
  stopifnot(nrow(mat_x) == nrow(mat_y), ncol(mat_x) == nrow(df_x), 
            ncol(mat_y) == nrow(df_y),
            length(unlist(list_end)) > 0)
  stopifnot(all(mat_x >= 0), all(mat_y >= 0))
  n <- nrow(mat_x); p1 <- ncol(mat_x)
  p2 <- ncol(mat_y); cell_name <- rownames(mat_x)
  
  # initialize
  df_res <- .init_chrom_df2(cell_name, vec_start, list_end)
  
  structure(list(mat_x = mat_x, 
                 mat_y = mat_y,
                 ht_map = ht_map,
                 df_res = df_res, 
                 snn = snn,
                 diffusion_dist = diffusion_dist),
            class = "chromatin_potential_prep")
}

.init_chrom_df2 <- function(cell_name, vec_start, list_end){
  tmp <- c(vec_start, unlist(list_end))
  n <- length(cell_name)
  df_res <- data.frame(idx = 1:n, 
                       init_state = rep(NA, n), 
                       num_cand = rep(0, n),
                       order_rec = rep(NA, n))
  rownames(df_res) <- cell_name
  
  df_res$init_state[which(cell_name %in% vec_start)] <- -1
  for(i in 1:length(list_end)){
    df_res$init_state[which(cell_name %in% list_end[[i]])] <- i
    df_res$order_rec[which(cell_name %in% list_end[[i]])] <- 0
  }
  
  df_res
}
