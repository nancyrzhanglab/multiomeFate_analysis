chromatin_potential_prepare2 <- function(mat_x, mat_y, df_x, df_y,
                                         snn,
                                         diffusion_dist,
                                         vec_start, list_end,
                                         ht_map = NA, options = list(), verbose = T){
  # [[NOTE: mat_x and mat_y here are the expressions of the metacells here]]
  stopifnot(nrow(mat_x) == nrow(mat_y), ncol(mat_x) == nrow(df_x), 
            ncol(mat_y) == nrow(df_y), is.list(options),
            length(unlist(list_end)) > 0)
  stopifnot(all(mat_x >= 0), all(mat_y >= 0))
  n <- nrow(mat_x); p1 <- ncol(mat_x)
  p2 <- ncol(mat_y); cell_name <- rownames(mat_x)
  
  # initialize
  df_res <- multiomeFate:::.init_chrom_df(n, vec_start, list_end, cell_name)
  if(full_options$est_options$enforce_cis){
    if(class(ht_map) != "hash" && all(is.na(ht_map))){
      full_options$est_options <- multiomeFate:::.gene_peak_map(df_x, df_y, full_options$est_options)
    } else {
      full_options$est_options$ht_map <- ht_map
    }
  }
  
  structure(list(mat_x = mat_x, 
                 mat_y = mat_y, 
                 df_x = df_x, 
                 df_y = df_y,
                 df_res = df_res, 
                 snn = snn,
                 diffusion_dist = diffusion_dist),
            class = "chromatin_potential_prep")
}
