rm(list=ls())
library(Seurat)
load("../../../../out/kevin/Writeup3e/10192021_sydney_preprocess.RData")
source("../Writeup3d/funcs.R")
source("select_cells.R")

###############

Seurat::DefaultAssay(all_data) <- "RNA"
all_data <- Seurat::NormalizeData(all_data)
all_data <- Seurat::ScaleData(all_data) 

###############

Seurat::Idents(all_data) <- "Original_condition"
table(Seurat::Idents(all_data))

tabulate_mat <- .tabulate_lineages(all_data)
head(tabulate_mat)
quantile(tabulate_mat[,"naive"], probs = seq(0,1,length.out=11))
max_val <- sapply(1:nrow(tabulate_mat), function(i){max(tabulate_mat[i,-5])})
tabulate_mat <- tabulate_mat[order(max_val, decreasing = T),]

#################

# first up -- cocl2
treatment_vec <- c("cocl2", "acid", "dab", "tram")
pval_thres <- c(0.05, 0.05, 0.05, 0.05)
de_list <- vector("list", 4)

table(all_data@meta.data$Original_condition)
for(kk in 1:length(treatment_vec)){
  print(treatment_vec[kk])
  treatment <- treatment_vec[kk]
  
  # set up identities
  res <- .select_expansion_naives(all_data, 
                                  tabulate_mat, 
                                  treatment = treatment,
                                  threshold = 3)
  newgroup <- all_data@meta.data$Original_condition
  names(newgroup) <- rownames(all_data@meta.data)
  newgroup[res$naive_terminal] <- paste0("naive_survive", treatment)
  newgroup[setdiff(res$naive_all, res$naive_terminal)] <- paste0("naive_no", treatment)
  all_data[[paste0(treatment, "Status")]] <- newgroup
  Seurat::Idents(all_data) <- paste0(treatment, "Status")
  
  set.seed(10)
  seurat_de <- Seurat::FindMarkers(all_data,
                                   assay = "RNA",
                                   slot = "scale.data",
                                   ident.1 = paste0("naive_survive", treatment),
                                   ident.2 = paste0("naive_no", treatment),
                                   min.diff.pct = 0.2,
                                   min.pct = 0.05,
                                   only.pos = T,
                                   verbose = T)
  print(dim(seurat_de))
  de_list[[kk]] <- seurat_de
  
  tentative_genes <- rownames(seurat_de)[which(seurat_de$p_val <= pval_thres[kk])]
  print(length(tentative_genes))
  # proportion_list <- .order_genes_by_threshold(genes = tentative_genes,
  #                                              naive_terminal = res$naive_terminal,
  #                                              seurat_obj = all_data,
  #                                              assay = "RNA")
  # marker_genes <- unique(c(names(proportion_list)[1:min(length(proportion_list), 15)], "EGFR", "NGFR", "NT5E"))
  
  marker_genes <- unique(c(tentative_genes[1:15], "EGFR", "NGFR", "NT5E"))
  marker_genes <- marker_genes[!is.na(marker_genes)]
  for(i in 1:length(marker_genes)){
    vec <- all_data[["RNA"]]@scale.data[ marker_genes[i],]
    idx1 <- which(Seurat::Idents(all_data) == paste0("naive_no", treatment))
    percentage1 <- round(length(which(vec[idx1] > 0))/length(idx1),2)
    idx2 <- which(Seurat::Idents(all_data) == paste0("naive_survive", treatment))
    percentage2 <- round(length(which(vec[idx2] > 0))/length(idx2),2)
    idx3 <- which(Seurat::Idents(all_data) == treatment)
    percentage3 <- round(length(which(vec[idx3] > 0))/length(idx3),2)
    
    plot1 <- Seurat::VlnPlot(all_data, features = marker_genes[i],
                             idents = c(paste0("naive_no", treatment), paste0("naive_survive", treatment), treatment),
                             assay = "RNA",
                             slot = "scale.data") + ggplot2::theme_classic()
    plot1 <- plot1 + ggplot2::scale_x_discrete(labels = c(
      paste0("naive - will not survive\n(",percentage1, "% of ", length(idx1), " cells)"),
      paste0("naive - will survive\n(",percentage2, "% of ", length(idx2), " cells)"),
      paste0(treatment, "\n(",percentage3, "% of ", length(idx3), " cells)")))
    plot1 <- plot1 + Seurat::NoLegend()
    ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup3e/10192021_sydney_naive_within", treatment, "_vln_", marker_genes[i], ".png"),
                    plot1, device = "png", width = 6, height = 5, units = "in")
  }
}
