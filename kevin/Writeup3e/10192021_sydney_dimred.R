rm(list=ls())
library(Seurat)
load("../../../../out/kevin/Writeup3e/10192021_sydney_preprocess.RData")
source("../Writeup3d/funcs.R")
source("select_cells.R")

zz <- all_data[["SCT"]]@scale.data
zz <- zz[,which(all_data$Original_condition == "naive")]
zz <- t(zz)
zz <- scale(zz)
zz <- t(zz)

svd_res <- irlba::irlba(zz, nv = 100)
svd_res$d
kk <- 100
smoothed_mat <- tcrossprod(.mult_mat_vec(svd_res$u[,1:kk], svd_res$d[1:kk]), svd_res$v[,1:kk]) 
rownames(smoothed_mat) <- rownames(zz)
colnames(smoothed_mat) <- colnames(zz)

smoothed_mat_full <- matrix(0, nrow(all_data[["SCT"]]@scale.data), ncol(all_data[["SCT"]]@counts))
rownames(smoothed_mat_full) <- rownames(all_data[["SCT"]]@scale.data)
colnames(smoothed_mat_full) <- colnames(all_data[["SCT"]]@counts)
for(i in 1:ncol(smoothed_mat)){
  smoothed_mat_full[,which(colnames(smoothed_mat_full) == colnames(smoothed_mat)[i])] <- smoothed_mat[,i]
}

all_data[["smooth"]] <- Seurat::CreateAssayObject(data = smoothed_mat_full)

tabulate_mat <- .tabulate_lineages(all_data)
apply(tabulate_mat, 2, quantile)
treatment_vec <- c("tram", "dab")

for(treatment in treatment_vec){
  res <- .select_expansion_naives(all_data, 
                                  tabulate_mat, 
                                  treatment = treatment,
                                  threshold = 3)
  newgroup <- all_data@meta.data$Original_condition
  names(newgroup) <- rownames(all_data@meta.data)
  newgroup[res$naive_terminal] <- paste0("naive_survive", treatment)
  newgroup[setdiff(res$naive_all, res$naive_terminal)] <- paste0("naive_no", treatment)
  table(newgroup)
  all_data[[paste0(treatment, "Status")]] <- newgroup
  Seurat::Idents(all_data) <- paste0(treatment, "Status")
  
  marker_genes <- c("EGFR", "NGFR", "NT5E")
  for(i in 1:length(marker_genes)){
    vec <- all_data[["smooth"]]@data[ marker_genes[i],]
    idx1 <- which(Seurat::Idents(all_data) == paste0("naive_no", treatment))
    percentage1 <- round(length(which(vec[idx1] > 0))/length(idx1),2)
    idx2 <- which(Seurat::Idents(all_data) == paste0("naive_survive", treatment))
    percentage2 <- round(length(which(vec[idx2] > 0))/length(idx2),2)

    plot1 <- Seurat::VlnPlot(all_data, 
                             features = marker_genes[i],
                             idents = c(paste0("naive_no", treatment), paste0("naive_survive", treatment)),
                             assay = "smooth",
                             slot = "data") + ggplot2::theme_classic()
    plot1 <- plot1 + ggplot2::scale_x_discrete(labels = c(
      paste0("naive - will not survive\n(",percentage1, "% of ", length(idx1), " cells)"),
      paste0("naive - will survive\n(",percentage2, "% of ", length(idx2), " cells)")))
    plot1 <- plot1 + Seurat::NoLegend()
    ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup3e/Writeup3e_sydney_smooth_naive_within", treatment, "_vln_", marker_genes[i], ".png"),
                    plot1, device = "png", width = 6, height = 5, units = "in")
  }
}