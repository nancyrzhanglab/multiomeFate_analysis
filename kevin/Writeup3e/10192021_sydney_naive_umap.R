# in this script: try finding the UMAP via highly-variable genes among the naive cells
# (as opposed to DE genes)
# apply this to the SCTransform scale.data [[subject to change since I'm not sure if this is the correct assay still]]

rm(list=ls())
# meant to be run from the laptop currently
library(Seurat)
source("../multiomeFate_analysis/kevin/Writeup3d/funcs.R")
source("../multiomeFate_analysis/kevin/Writeup3e/select_cells.R")

load("../../dbox_MultiomeFate/data/ShafferLab/10192021_kevin_preprocess/10192021_sydney_preprocess.RData")
tabulate_mat <- .tabulate_lineages(all_data)
max_val <- sapply(1:nrow(tabulate_mat), function(i){max(tabulate_mat[i,-5])})
tabulate_mat2 <- tabulate_mat[order(max_val, decreasing = T),]
head(tabulate_mat2)

###########

naive_mat <- all_data[["RNA"]]@counts[,which(all_data@meta.data[,"Original_condition"] == "naive")]
naive_data <- Seurat::CreateSeuratObject(counts = naive_mat)
Seurat::DefaultAssay(naive_data) <- "RNA"
naive_data <- Seurat::SCTransform(naive_data, clip.range = c(-30,30))
quantile(naive_data[["SCT"]]@scale.data)
naive_data <- Seurat::RunPCA(naive_data, verbose = F)
naive_data[["Original_condition"]] <- rep("naive", ncol(naive_data))

set.seed(10)
naive_data <- Seurat::RunUMAP(naive_data, dims = 1:50)
Seurat::DimPlot(naive_data, reduction = "umap")

gene_vec <- c("EGFR", "NGFR", "NT5E")
for(gene in gene_vec){
  Seurat::FeaturePlot(naive_data,
                      features = gene,
                      reduction = "umap")
}

lin_mat <-all_data[["lineage"]]@counts
lin_mat <- lin_mat[,which(colnames(lin_mat) %in% colnames(naive_data))]
naive_data[["lineage"]] <- Seurat::CreateAssayObject(counts = lin_mat)
naive_data2 <- .collapse_by_lineage(naive_data, assay = "SCT", slot = "scale.data")
naive_data2[["Original_condition"]] <- rep("naive", ncol(naive_data2))
#############

# for a particular treatment, plot cells based on their
# lineage's ratio of (log(treated/naive), which comes from
# treated = naive*exp(growth))

treatment <- "dab"

growth_factor <- rep(0, ncol(naive_data))
names(growth_factor) <- colnames(naive_mat)
lin_mat <- Matrix::t(all_data[["lineage"]]@counts)
for(lin_idx in 1:nrow(tabulate_mat)){
  if(lin_idx %% floor(nrow(tabulate_mat)/10) == 0) cat('*')
  
  value <- log((tabulate_mat[lin_idx, treatment]+.5) / tabulate_mat[lin_idx,"naive"])
  
  lin_idx2 <- which(colnames(lin_mat) == rownames(tabulate_mat)[lin_idx])
  cell_idx <- .nonzero_col(lin_mat, lin_idx2)
  cell_names <- colnames(all_data)[cell_idx]
  growth_factor[intersect(cell_names, colnames(naive_mat))] <- value
}
quantile(growth_factor)

naive_data[[paste0(treatment, "_growth")]] <- growth_factor
Seurat::FeaturePlot(naive_data, features = paste0(treatment, "_growth"), reduction = "umap")

################################

# when looking at cells
for(treatment in treatment_vec){
  res <- .select_expansion_naives(naive_data, 
                                  tabulate_mat, 
                                  treatment = treatment,
                                  threshold = 3, 
                                  type = 1)
  naive_deathtoall <- .select_naive_deathtoall(naive_data,
                                               tabulate_mat,
                                               threshold = 1, 
                                               type = 1)
  newgroup <- rep("none", ncol(naive_data))
  names(newgroup) <- rownames(naive_data@meta.data)
  newgroup[res$naive_terminal] <- paste0("naive_survive", treatment)
  newgroup[naive_deathtoall] <- paste0("naive_no", treatment)
  table(newgroup)
  naive_data[[paste0(treatment, "Status")]] <- newgroup
  Seurat::Idents(naive_data) <- paste0(treatment, "Status")
  
  marker_genes <- c("EGFR", "NGFR", "NT5E")
  for(i in 1:length(marker_genes)){
    vec <- naive_data[["SCT"]]@scale.data[marker_genes[i],]
    idx1 <- which(Seurat::Idents(naive_data) == paste0("naive_no", treatment))
    percentage1 <- round(length(which(vec[idx1] > 0))/length(idx1),2)
    idx2 <- which(Seurat::Idents(naive_data) == paste0("naive_survive", treatment))
    percentage2 <- round(length(which(vec[idx2] > 0))/length(idx2),2)
    idx3 <- which(Seurat::Idents(naive_data) == treatment)
    percentage3 <- round(length(which(vec[idx3] > 0))/length(idx3),2)
    
    plot1 <- Seurat::VlnPlot(naive_data, features = marker_genes[i],
                             idents = c(paste0("naive_no", treatment), paste0("naive_survive", treatment), treatment),
                             assay = "SCT",
                             slot = "data") + ggplot2::theme_classic()
    plot1 <- plot1 + ggplot2::scale_x_discrete(labels = c(
      paste0("naive - will not survive\n(",percentage1, "% of ", length(idx1), " cells)"),
      paste0("naive - will survive\n(",percentage2, "% of ", length(idx2), " cells)"),
      paste0(treatment, "\n(",percentage3, "% of ", length(idx3), " cells)")))
    plot1 <- plot1 + Seurat::NoLegend()
    plot1
  }
}


# when looking at lineages
for(treatment in treatment_vec){
  res <- .select_expansion_naives(naive_data2, 
                                  tabulate_mat, 
                                  treatment = treatment,
                                  threshold = 3, 
                                  type = 2)
  naive_deathtoall <- .select_naive_deathtoall(naive_data2,
                                               tabulate_mat,
                                               threshold = 1, 
                                               type = 2)
  newgroup <- rep("none", ncol(naive_data2))
  names(newgroup) <- rownames(naive_data2@meta.data)
  newgroup[res$naive_terminal] <- paste0("naive_survive", treatment)
  newgroup[naive_deathtoall] <- paste0("naive_no", treatment)
  table(newgroup)
  naive_data2[[paste0(treatment, "Status")]] <- newgroup
  Seurat::Idents(naive_data2) <- paste0(treatment, "Status")
  
  marker_genes <- c("EGFR", "NGFR", "NT5E")
  for(i in 1:length(marker_genes)){
    vec <- naive_data2[["RNA"]]@counts[marker_genes[i],]
    idx1 <- which(Seurat::Idents(naive_data2) == paste0("naive_no", treatment))
    percentage1 <- round(length(which(vec[idx1] > 0))/length(idx1),2)
    idx2 <- which(Seurat::Idents(naive_data2) == paste0("naive_survive", treatment))
    percentage2 <- round(length(which(vec[idx2] > 0))/length(idx2),2)
    idx3 <- which(Seurat::Idents(naive_data2) == treatment)
    percentage3 <- round(length(which(vec[idx3] > 0))/length(idx3),2)
    
    plot1 <- Seurat::VlnPlot(naive_data2, features = marker_genes[i],
                             idents = c(paste0("naive_no", treatment), paste0("naive_survive", treatment), treatment),
                             assay = "RNA",
                             slot = "counts") + ggplot2::theme_classic()
    plot1 <- plot1 + ggplot2::scale_x_discrete(labels = c(
      paste0("naive - will not survive\n(",percentage1, "% of ", length(idx1), " cells)"),
      paste0("naive - will survive\n(",percentage2, "% of ", length(idx2), " cells)"),
      paste0(treatment, "\n(",percentage3, "% of ", length(idx3), " cells)")))
    plot1 <- plot1 + Seurat::NoLegend()
    plot1
  }
}


