rm(list=ls())

load("../../../../out/kevin/Writeup3c/20210823_10x_embryo_result.RData")

# for now, a small percentage of cells haven't been matched yet
# let's remove those cells
zz <- which(apply(res$df_res, 1, function(x){is.na(x["init_state"]) & is.na(x["order_rec"])}))
stopifnot(length(zz) == 0)

#####################################

# let's do the naive fate calculations first
list_diagnos <- res$list_diagnos
n <- nrow(mat_x)
celltype <- as.numeric(as.character(celltype))
end_states <- list("16", "6", c("1","2","4"), "9")

# select endstate
prob_mat <- sapply(end_states, function(end_state){
  prob_vec <- rep(NA, n)
  prob_vec[which(celltype %in% end_state)] <- 1
  
  # set all other endstates to prob = 0
  prob_vec[which(celltype %in% end_states[which(!end_states %in% end_state)])] <- 0
  
  # iterate backward
  total_iter <- 1
  while(any(is.na(prob_vec)) && total_iter <= 10){
    print(sum(is.na(prob_vec)))
    
    for(iter in 1:length(list_diagnos)){
      for(i in 1:length(list_diagnos[[iter]]$recruit$postprocess)){
        idx_from <- list_diagnos[[iter]]$recruit$postprocess[[i]]$from
        idx_to <- list_diagnos[[iter]]$recruit$postprocess[[i]]$to
        
        ## ignore all "another nodes" that don't have a value assigned yet
        valid_to <- idx_to[which(!is.na(prob_vec[idx_to]))]
        valid_from <- idx_from[which(is.na(prob_vec[idx_from]))]
        
        if(length(valid_from) == 0 | length(valid_to) == 0) next()
        ## do multiplication: prob of a node going into another node * prob of that node going to end state
        prob_vec[valid_from] <- sum(prob_vec[valid_to])/length(idx_to)
        prob_vec[setdiff(idx_to, valid_to)] <- sum(prob_vec[valid_to])/length(idx_to)
      }
    }
    
    total_iter <- total_iter + 1
  }
  
  prob_vec
})

sum(is.na(prob_mat))
prob_mat[is.na(prob_mat)] <- 0
prob_mat <- t(apply(prob_mat, 1, function(x){
  if(all(x == 0)) return(x)
  x/sum(x)
}))

#####################################

rownames(prob_mat) <- rownames(mbrain3@meta.data)
colnames(prob_mat) <- paste0("prob_", 1:ncol(prob_mat))
mbrain3[["prob"]] <- Seurat::CreateDimReducObject(embedding = prob_mat, key = "prob_")

plot1 <- Seurat::FeaturePlot(mbrain3, features = "prob_1", reduction = "wnn.umap")
plot1 <- plot1 + ggplot2::ggtitle("To Oligodendrocyte")
plot2 <- Seurat::FeaturePlot(mbrain3, features = "prob_2", reduction = "wnn.umap")
plot2 <- plot2 + ggplot2::ggtitle("To Forebrain gabaergic")
plot3 <- Seurat::FeaturePlot(mbrain3, features = "prob_3", reduction = "wnn.umap")
plot3 <- plot3 + ggplot2::ggtitle("To one of the Corticals")
plot4 <- Seurat::FeaturePlot(mbrain3, features = "prob_4", reduction = "wnn.umap")
plot4 <- plot4 + ggplot2::ggtitle("To the other Cortical")

tmp2 <- cowplot::plot_grid(plot1, plot2, plot3, plot4)
cowplot::save_plot(filename =  "../../../../out/figures/Writeup3c/Writeup3c_10x_embryo_de_fateprob_naive.png",
                   tmp2, ncol = 2, nrow = 2, base_height = 3.5, base_asp = 4/3, device = "png")

##################################

Seurat::DefaultAssay(mbrain3) <- "SCT"
for(i in 1:length(de_combined)){
  gene_vec <- rownames(de_combined[[i]])[which(rownames(de_combined[[i]])[1:50] %in% colnames(mat_y))]
  gene_vec <- gene_vec[1:min(length(gene_vec), 9)]
  plot1 <- Seurat::FeaturePlot(mbrain3, features = gene_vec, reduction = "wnn.umap")
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup3c/Writeup3c_10x_embryo_de_", names(de_combined)[i], ".png"), 
                  plot1, device = "png", width = 12, height = 10, units = "in")
}
