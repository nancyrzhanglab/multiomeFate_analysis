rm(list=ls())
library(Seurat)
library(Signac)
library(igraph)

load("../../../../out/kevin/Writeup6/Writeup6_DE_day10-week5.RData")
load("../../../../out/kevin/Writeup6/Writeup6_all-data_lineage-assigned.RData")

for(i in 1:length(de_day10vsweek5)){
  print(names(de_day10vsweek5)[i])
  tmp <- length(which(de_day10vsweek5[[i]][,"p_val"] <= 1e-6))
  print(tmp)
  print("====")
}

############

for(i in 1:length(de_ExpandvsShrunk)){
  print(names(de_ExpandvsShrunk)[i])
  tmp <- length(which(de_ExpandvsShrunk[[i]][,"p_val"] <= 1e-6))
  print(tmp)
  print("====")
}

#############

treatment_vec <- names(de_week5_pairwise)
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
max_size <- max(tab_mat)
size_bin <- seq(log(10), log(max_size), length.out = 100)
size_vec <- seq(3,20,length.out=100)

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
  
  # apply(adj_mat, 2, function(x){round(quantile(x),2)})
  
  adj_mat2 <- (1-adj_mat)^10
  g <- igraph::graph_from_adjacency_matrix(adj_mat2, mode = "undirected",
                                           weighted = T, diag = F)
  # quantile(igraph::E(g)$weight)
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
  
  lineage_names <- rownames(adj_mat2)
  if(kk == 1){
    col_idx <- which(colnames(tab_mat) == "week5_CIS")
  } else if(kk == 2){
    col_idx <- which(colnames(tab_mat) == "week5_COCL2")
  } else {
    col_idx <- which(colnames(tab_mat) == "week5_DABTRAM")
  }
  tmp <- sapply(lineage_names, function(lineage_name){
    size_vec[which.min(abs(size_bin - log(tab_mat[lineage_name,col_idx])))]
  })
  igraph::V(g)$size <- tmp
  
  png(paste0("../../../../out/figures/Writeup6/Writeup6_DE_", treatment_vec[kk],"-expanded_graph_all.png"),
      height = 2500, width = 2500, units = "px", res = 300)
  par(mar = rep(0.5,4))
  set.seed(10)
  plot(g, 
       vertex.size = igraph::V(g)$size, 
       edge.width = 3*igraph::E(g)$weight)
  graphics.off()
}

############################################

for(kk in 1:3){
  tmp <- quantile(sapply(de_day10_pairwise[[kk]]$de_list, function(x){
    length(which(x[,"p_val"] <= 1e-4))
  }))
  print(tmp)
}

max_size_day10 <- max(tab_mat[,c("day10_CIS", "day10_COCL2", "day10_DABTRAM")])
size_bin_day10 <- seq(log(10), log(max_size_day10), length.out = 100)
max_size_week5 <- max(tab_mat)
size_bin_week5 <- seq(log(10), log(max_size_week5), length.out = 100)
size_vec <- seq(3,20,length.out=100)
treatment_vec <- names(de_week5_pairwise)
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
for(kk in 1:3){
  print(kk)
  combn_mat <- de_day10_pairwise[[kk]]$combn_mat
  de_list <- de_day10_pairwise[[kk]]$de_list
  lin_uniq <- sort(unique(as.vector(combn_mat)))
  num_uniq <- length(lin_uniq)
  adj_mat <- matrix(0, num_uniq, num_uniq)
  colnames(adj_mat) <- lin_uniq
  rownames(adj_mat) <- lin_uniq
  thres <- 1e-4
  for(x in 1:ncol(combn_mat)){
    i <- which(lin_uniq == combn_mat[1,x])
    j <- which(lin_uniq == combn_mat[2,x])
    val <- length(which(de_list[[x]][,"p_val"] <= 0.05))
    adj_mat[i,j] <- val
    adj_mat[j,i] <- val
  }
  adj_mat <- adj_mat/pmax(max(adj_mat),1)
  
  hclust_res <- stats::hclust(stats::dist(adj_mat))
  adj_mat <- adj_mat[hclust_res$order, hclust_res$order]
  
  png(paste0("../../../../out/figures/Writeup6/Writeup6_DE_day10_", treatment_vec[kk],"-expanded_heatmap_all.png"),
      height = 2500, width = 2500, units = "px", res = 300)
  par(mar = rep(0.5,4))
  image(adj_mat, asp = T)
  graphics.off()
  
  adj_mat2 <- (1-adj_mat)^10
  g <- igraph::graph_from_adjacency_matrix(adj_mat2, mode = "undirected",
                                           weighted = T, diag = F)
  # quantile(igraph::E(g)$weight)
  col_palette <- grDevices::colorRampPalette(c("white","firebrick1"))(100)
  lineage_names <- rownames(adj_mat2)
  if(kk == 1){
    col_idx <- which(colnames(tab_mat) == "week5_CIS")
  } else if(kk == 2){
    col_idx <- which(colnames(tab_mat) == "week5_COCL2")
  } else {
    col_idx <- which(colnames(tab_mat) == "week5_DABTRAM")
  }
  col_vec <- sapply(lineage_names, function(lineage_name){
    col_palette[which.min(abs(size_bin_week5 - log(tab_mat[lineage_name,col_idx])))]
  })
  igraph::V(g)$color <- col_vec
  
  if(kk == 1){
    col_idx <- which(colnames(tab_mat) == "day10_CIS")
  } else if(kk == 2){
    col_idx <- which(colnames(tab_mat) == "day10_COCL2")
  } else {
    col_idx <- which(colnames(tab_mat) == "day10_DABTRAM")
  }
  tmp <- sapply(lineage_names, function(lineage_name){
    size_vec[which.min(abs(size_bin_day10 - log(tab_mat[lineage_name,col_idx])))]
  })
  igraph::V(g)$size <- tmp
  
  png(paste0("../../../../out/figures/Writeup6/Writeup6_DE_day10_", treatment_vec[kk],"-expanded_graph_all.png"),
      height = 2500, width = 2500, units = "px", res = 300)
  par(mar = rep(0.5,4))
  set.seed(10)
  plot(g, 
       vertex.size = igraph::V(g)$size, 
       edge.width = 3*igraph::E(g)$weight)
  graphics.off()
}



