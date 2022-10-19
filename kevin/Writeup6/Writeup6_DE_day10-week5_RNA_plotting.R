rm(list=ls())
library(Seurat)
library(Signac)
library(igraph)

load("../../../../out/kevin/Writeup6/Writeup6_DE_day10-week5.RData")
sapply(1:length(de_day10vsweek5), function(i){
  print(names(de_day10vsweek5)[i])
  tmp <- length(which(de_day10vsweek5[[i]][,"p_val"] <= 1e-5))
  print(tmp)
  invisible()
})

############

sapply(1:length(de_ExpandvsShrunk), function(i){
  print(names(de_ExpandvsShrunk)[i])
  tmp <- length(which(de_ExpandvsShrunk[[i]][,"p_val"] <= 1e-5))
  print(tmp)
  invisible()
})

#############
treatment_vec <- names(de_week5_pairwise)

for(kk in 1:3){
  print(kk)
  combn_mat <- de_week5_pairwise[[kk]]$combn_mat
  de_list <- de_week5_pairwise[[kk]]$de_list
  lin_uniq <- sort(unique(as.vector(combn_mat)))
  num_uniq <- length(lin_uniq)
  adj_mat <- matrix(0, num_uniq, num_uniq)
  colnames(adj_mat) <- lin_uniq
  rownames(adj_mat) <- lin_uniq
  thres <- 1e-4
  for(x in 1:ncol(combn_mat)){
    i <- which(lin_uniq == combn_mat[1,x])
    j <- which(lin_uniq == combn_mat[2,x])
    val <- length(which(de_list[[x]][,"p_val"] <= 1e-5))
    adj_mat[i,j] <- val
    adj_mat[j,i] <- val
  }
  adj_mat <- adj_mat/max(adj_mat)
  
  hclust_res <- stats::hclust(stats::dist(adj_mat))
  adj_mat <- adj_mat[hclust_res$order, hclust_res$order]
  
  png(paste0("../../../../out/figures/Writeup6/Writeup6_DE_", treatment_vec[kk],"-expanded_heatmap_all.png"),
      height = 2500, width = 2500, units = "px", res = 300)
  par(mar = rep(0.5,4))
  image(adj_mat, asp = T)
  graphics.off()
  
  apply(adj_mat, 2, function(x){round(quantile(x),2)})
  
  adj_mat2 <- (1-adj_mat)^7
  g <- igraph::graph_from_adjacency_matrix(adj_mat2, mode = "undirected",
                                           weighted = T, diag = F)
  quantile(igraph::E(g)$weight)
  col_vec <- scales::hue_pal()(150)
  if(kk == 1){
    igraph::V(g)$color <- c("white", "black", col_vec[c(1:6)+75], col_vec[3*c(1:17)])
  } else if(kk == 2){
    igraph::V(g)$color <- c("white", "black", col_vec[3*c(1:5)], col_vec[c(1:11)+75], col_vec[c(1:11)+120])
  } else if(kk == 3){
    igraph::V(g)$color <- c(grDevices::colorRampPalette(c("white","black"))(6), 
                            col_vec[c(1:12)+120],
                            col_vec[c(1:31)+75], 
                            col_vec[3*c(1:6)])
  }
  
  png(paste0("../../../../out/figures/Writeup6/Writeup6_DE_", treatment_vec[kk],"-expanded_graph_all.png"),
      height = 2500, width = 2500, units = "px", res = 300)
  par(mar = rep(0.5,4))
  set.seed(10)
  plot(g, 
       vertex.size = 5, vertex.label = "",
       edge.width = 3*igraph::E(g)$weight)
  graphics.off()
}


