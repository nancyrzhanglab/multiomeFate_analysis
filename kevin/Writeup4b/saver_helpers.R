# modified from https://github.com/mohuangx/SAVER/blob/master/R/cor_adjust.R

compute_cor_genes <- function(x, cell_idx = NULL,
                              gene_idx = NULL) {
 
  if(all(is.null(cell_idx))){
    cell_idx <- 1:ncol(x$estimate)
  }
  if(all(is.null(gene_idx))){
    gene_idx <- 1:nrow(x$estimate)
  }
  ngenes <- length(gene_idx)
  adj_vec <- rep(0, ngenes)
  
  cor_mat <- cor(t(x$estimate[gene_idx,cell_idx]))
  for (i in 1:ngenes) {
    adj_vec[i] <- sqrt(var(x$estimate[gene_idx[i],cell_idx], na.rm = TRUE)/
                         (var(x$estimate[gene_idx[i],cell_idx], na.rm = TRUE) +
                            mean(x$se[gene_idx[i],cell_idx]^2, na.rm = TRUE)))
  }
  adj_mat <- outer(adj_vec, adj_vec)
  cor_adj <- adj_mat*cor_mat

  cor_adj
}