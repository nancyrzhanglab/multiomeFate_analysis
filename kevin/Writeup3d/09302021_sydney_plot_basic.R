rm(list=ls())
library(Seurat); library(Signac)
source("funcs.R")

load("../../../../data/Sydney_stressors_2021-09-24/all_data_SCT.RData")

lin_mat <- all_data[["lineage"]]@counts
lin_mat@x <- rep(1, length(lin_mat@x))
dim(lin_mat)
lin_mat[1:5,1:5]
count_vec <- sparseMatrixStats::rowSums2(lin_mat)
lin_mat <- lin_mat[count_vec > 1, ] # since we don't have the barcodes themselves, the best we can do right now is focus on the lineages with many cells
sum_vec <- sparseMatrixStats::colSums2(lin_mat)
lin_mat <- lin_mat[,sum_vec > 0] 
dim(lin_mat)

#####################

while(TRUE){
  print(dim(lin_mat))
  bool <- TRUE
  
  # need lineages to be in the columns
  lin_mat <- Matrix::t(lin_mat)
  lin_idx_list <- lapply(1:ncol(lin_mat), function(j){
    .nonzero_col(lin_mat, j)
  })
  names(lin_idx_list) <- colnames(lin_mat)
  factor_vec <- as.factor(all_data@meta.data$Original_condition)
  tabulate_mat <- t(sapply(lin_idx_list, function(idx){
    table(factor_vec[idx])
  }))
  rownames(tabulate_mat) <- colnames(lin_mat)
  naive_idx <- which(colnames(tabulate_mat) == "naive")
  if(any(which(tabulate_mat[,naive_idx] == 0))) bool <- FALSE
  tabulate_mat <- tabulate_mat[which(tabulate_mat[,naive_idx] != 0),]
  
  # need lineage to be in the rows
  lin_mat <- Matrix::t(lin_mat)
  lin_mat <- lin_mat[rownames(lin_mat) %in% rownames(tabulate_mat),]
  sum_vec <- sparseMatrixStats::rowSums2(lin_mat)
  if(any(sum_vec < 2)) bool <- FALSE
  lin_mat <- lin_mat[sum_vec >= 2,] 
  
  sum_vec <- sparseMatrixStats::colSums2(lin_mat)
  if(any(sum_vec == 0)) bool <- FALSE
  lin_mat <- lin_mat[,sum_vec > 0] 
  
  if(bool) break()
}

tabulate_mat <- tabulate_mat[order(tabulate_mat[,"naive"], decreasing = T),]
tabulate_mat[1:100,]

###############################3

n <- ncol(all_data[["SCT"]])
keep_vec <- rep(0, n)
keep_vec[rownames(all_data@meta.data) %in% colnames(lin_mat)] <- 1
all_data[["keep"]] <- keep_vec
all_data <- subset(all_data, keep == 1)
all(colnames(all_data) == colnames(lin_mat))

###########################

col_palette <- scales::hue_pal()(6)
lineage_vec <- c("Lin26933", "Lin216958", "Lin331208", "Lin611467", "Lin283219", "Lin156893", "Lin206249")
xlim <- range(all_data[["umap"]]@cell.embeddings[,1])
ylim <- range(all_data[["umap"]]@cell.embeddings[,2])

for(lineage in lineage_vec){
  cells <- which(lin_mat[lineage,] != 0)
  col_vec <- col_palette[which(sort(unique(all_data@meta.data$Original_condition)) %in% all_data@meta.data$Original_condition[cells])]
  plot1 <- Seurat::DimPlot(all_data, reduction = "umap", group.by = "Original_condition", 
                           label = TRUE, repel = TRUE, label.size = 2.5,
                           cols = col_vec,
                           cells = cells) + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  plot1 <- plot1 + ggplot2::ggtitle(paste0("Lineage ", lineage, "\nColored by treatment"))
  plot1 <- plot1 + ggplot2::expand_limits(x = xlim, y = ylim)
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup3d/Writeup3d_sydney_", lineage, ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}

####################



