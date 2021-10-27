rm(list=ls())
# meant to be run from the laptop currently
library(Seurat)
source("../multiomeFate_analysis/kevin/Writeup3d/funcs.R")
source("../multiomeFate_analysis/kevin/Writeup3e/select_cells.R")

load("../../dbox_MultiomeFate/data/ShafferLab/10192021_kevin_preprocess/10192021_sydney_preprocess.RData")
tabulate_mat <- .tabulate_lineages(all_data)
all_data2 <- .collapse_by_lineage(all_data, assay = "SCT", slot = "scale.data")
all_data3 <- .collapse_by_lineage(all_data, assay = "SCT", slot = "data")

Seurat::DefaultAssay(all_data) <- "RNA"
all_data <- Seurat::NormalizeData(all_data)
all_data <- Seurat::ScaleData(all_data) 
all_data4 <- .collapse_by_lineage(all_data, assay = "RNA", slot = "scale.data")

##################

# experiment 1: try SCTransform, but for naive that survive vs naive that die to everything
tabulate_mat <- .tabulate_lineages(all_data)
apply(tabulate_mat, 2, quantile)
treatment_vec <- c("tram", "dab")

for(treatment in treatment_vec){
  res <- .select_expansion_naives(all_data, 
                                  tabulate_mat, 
                                  treatment = treatment,
                                  threshold = 3)
  naive_deathtoall <- .select_naive_deathtoall(all_data,
                                               tabulate_mat,
                                               threshold = 1)
  newgroup <- all_data@meta.data$Original_condition
  names(newgroup) <- rownames(all_data@meta.data)
  newgroup[res$naive_terminal] <- paste0("naive_survive", treatment)
  newgroup[naive_deathtoall] <- paste0("naive_no", treatment)
  table(newgroup)
  all_data[[paste0(treatment, "Status")]] <- newgroup
  Seurat::Idents(all_data) <- paste0(treatment, "Status")
  
  marker_genes <- c("EGFR", "NGFR", "NT5E")
  for(i in 1:length(marker_genes)){
    vec <- all_data[["SCT"]]@scale.data[ marker_genes[i],]
    idx1 <- which(Seurat::Idents(all_data) == paste0("naive_no", treatment))
    percentage1 <- round(length(which(vec[idx1] > 0))/length(idx1),2)
    idx2 <- which(Seurat::Idents(all_data) == paste0("naive_survive", treatment))
    percentage2 <- round(length(which(vec[idx2] > 0))/length(idx2),2)
    idx3 <- which(Seurat::Idents(all_data) == treatment)
    percentage3 <- round(length(which(vec[idx3] > 0))/length(idx3),2)
    
    plot1 <- Seurat::VlnPlot(all_data, features = marker_genes[i],
                             idents = c(paste0("naive_no", treatment), paste0("naive_survive", treatment), treatment),
                             assay = "SCT",
                             slot = "scale.data") + ggplot2::theme_classic()
    plot1 <- plot1 + ggplot2::scale_x_discrete(labels = c(
      paste0("naive - will not survive\n(",percentage1, "% of ", length(idx1), " cells)"),
      paste0("naive - will survive\n(",percentage2, "% of ", length(idx2), " cells)"),
      paste0(treatment, "\n(",percentage3, "% of ", length(idx3), " cells)")))
    plot1 <- plot1 + Seurat::NoLegend()
    ggplot2::ggsave(filename = paste0("../../out/figures/Writeup3e/Writeup3e_sydney_naive_within", treatment, "_vln_", marker_genes[i], "_sct-scale_single_deathtoall.png"),
                    plot1, device = "png", width = 6, height = 5, units = "in")
  }
}

##########################

# experiment 2: try SCTransform, but for naive that survive vs naive that die to everything
# AND average across lineages
tabulate_mat <- .tabulate_lineages(all_data)
treatment_vec <- c("tram", "dab")

for(treatment in treatment_vec){
  res <- .select_expansion_naives(all_data2, 
                                  tabulate_mat, 
                                  treatment = treatment,
                                  threshold = 3,
                                  type = 2)
  naive_deathtoall <- .select_naive_deathtoall(all_data2,
                                               tabulate_mat,
                                               threshold = 1,
                                               type = 2)
  newgroup <- all_data@meta.data$Original_condition
  names(newgroup) <- rownames(all_data@meta.data)
  newgroup[res$naive_terminal] <- paste0("naive_survive", treatment)
  newgroup[naive_deathtoall] <- paste0("naive_no", treatment)
  table(newgroup)
  all_data2[[paste0(treatment, "Status")]] <- newgroup
  Seurat::Idents(all_data2) <- paste0(treatment, "Status")
  
  marker_genes <- c("EGFR", "NGFR", "NT5E")
  for(i in 1:length(marker_genes)){
    vec <- all_data2[["RNA"]]@counts[ marker_genes[i],]
    idx1 <- which(Seurat::Idents(all_data2) == paste0("naive_no", treatment))
    percentage1 <- round(length(which(vec[idx1] > 0))/length(idx1),2)
    idx2 <- which(Seurat::Idents(all_data2) == paste0("naive_survive", treatment))
    percentage2 <- round(length(which(vec[idx2] > 0))/length(idx2),2)
    idx3 <- which(Seurat::Idents(all_data2) == treatment)
    percentage3 <- round(length(which(vec[idx3] > 0))/length(idx3),2)
    
    plot1 <- Seurat::VlnPlot(all_data2, features = marker_genes[i],
                             idents = c(paste0("naive_no", treatment), paste0("naive_survive", treatment), treatment),
                             assay = "RNA",
                             slot = "counts") + ggplot2::theme_classic()
    plot1 <- plot1 + ggplot2::scale_x_discrete(labels = c(
      paste0("naive - will not survive\n(",percentage1, "% of ", length(idx1), " cells)"),
      paste0("naive - will survive\n(",percentage2, "% of ", length(idx2), " cells)"),
      paste0(treatment, "\n(",percentage3, "% of ", length(idx3), " cells)")))
    plot1 <- plot1 + Seurat::NoLegend()
    ggplot2::ggsave(filename = paste0("../../out/figures/Writeup3e/Writeup3e_sydney_naive_within", treatment, "_vln_", marker_genes[i], "_sct-scale_average_deathtoall.png"),
                    plot1, device = "png", width = 6, height = 5, units = "in")
  }
}

##########################

# experiment 3: try log1p, but for naive that survive vs naive that die to everything
# AND average across lineages
tabulate_mat <- .tabulate_lineages(all_data)
treatment_vec <- c("tram", "dab")

for(treatment in treatment_vec){
  res <- .select_expansion_naives(all_data3, 
                                  tabulate_mat, 
                                  treatment = treatment,
                                  threshold = 3,
                                  type = 2)
  naive_deathtoall <- .select_naive_deathtoall(all_data3,
                                               tabulate_mat,
                                               threshold = 1,
                                               type = 2)
  newgroup <- all_data3@meta.data$Original_condition
  names(newgroup) <- rownames(all_data3@meta.data)
  newgroup[res$naive_terminal] <- paste0("naive_survive", treatment)
  newgroup[naive_deathtoall] <- paste0("naive_no", treatment)
  table(newgroup)
  all_data3[[paste0(treatment, "Status")]] <- newgroup
  Seurat::Idents(all_data3) <- paste0(treatment, "Status")
  
  marker_genes <- c("EGFR", "NGFR", "NT5E")
  for(i in 1:length(marker_genes)){
    vec <- all_data3[["RNA"]]@counts[ marker_genes[i],]
    idx1 <- which(Seurat::Idents(all_data3) == paste0("naive_no", treatment))
    percentage1 <- round(length(which(vec[idx1] > 0))/length(idx1),2)
    idx2 <- which(Seurat::Idents(all_data3) == paste0("naive_survive", treatment))
    percentage2 <- round(length(which(vec[idx2] > 0))/length(idx2),2)
    idx3 <- which(Seurat::Idents(all_data3) == treatment)
    percentage3 <- round(length(which(vec[idx3] > 0))/length(idx3),2)
    
    plot1 <- Seurat::VlnPlot(all_data3, features = marker_genes[i],
                             idents = c(paste0("naive_no", treatment), paste0("naive_survive", treatment), treatment),
                             assay = "RNA",
                             slot = "counts") + ggplot2::theme_classic()
    plot1 <- plot1 + ggplot2::scale_x_discrete(labels = c(
      paste0("naive - will not survive\n(",percentage1, "% of ", length(idx1), " cells)"),
      paste0("naive - will survive\n(",percentage2, "% of ", length(idx2), " cells)"),
      paste0(treatment, "\n(",percentage3, "% of ", length(idx3), " cells)")))
    plot1 <- plot1 + Seurat::NoLegend()
    ggplot2::ggsave(filename = paste0("../../out/figures/Writeup3e/Writeup3e_sydney_naive_within", treatment, "_vln_", marker_genes[i], "_log1p_average_deathtoall.png"),
                    plot1, device = "png", width = 6, height = 5, units = "in")
  }
}


##########################

# experiment 4: try log-normalize, but for naive that survive vs naive that die to everything
# AND average across lineages
tabulate_mat <- .tabulate_lineages(all_data)
treatment_vec <- c("tram", "dab")

for(treatment in treatment_vec){
  res <- .select_expansion_naives(all_data4, 
                                  tabulate_mat, 
                                  treatment = treatment,
                                  threshold = 3,
                                  type = 2)
  naive_deathtoall <- .select_naive_deathtoall(all_data4,
                                               tabulate_mat,
                                               threshold = 1,
                                               type = 2)
  newgroup <- all_data4@meta.data$Original_condition
  names(newgroup) <- rownames(all_data4@meta.data)
  newgroup[res$naive_terminal] <- paste0("naive_survive", treatment)
  newgroup[naive_deathtoall] <- paste0("naive_no", treatment)
  table(newgroup)
  all_data4[[paste0(treatment, "Status")]] <- newgroup
  Seurat::Idents(all_data4) <- paste0(treatment, "Status")
  
  marker_genes <- c("EGFR", "NGFR", "NT5E")
  for(i in 1:length(marker_genes)){
    vec <- all_data4[["RNA"]]@counts[ marker_genes[i],]
    idx1 <- which(Seurat::Idents(all_data4) == paste0("naive_no", treatment))
    percentage1 <- round(length(which(vec[idx1] > 0))/length(idx1),2)
    idx2 <- which(Seurat::Idents(all_data4) == paste0("naive_survive", treatment))
    percentage2 <- round(length(which(vec[idx2] > 0))/length(idx2),2)
    idx3 <- which(Seurat::Idents(all_data4) == treatment)
    percentage3 <- round(length(which(vec[idx3] > 0))/length(idx3),2)
    
    plot1 <- Seurat::VlnPlot(all_data4, features = marker_genes[i],
                             idents = c(paste0("naive_no", treatment), paste0("naive_survive", treatment), treatment),
                             assay = "RNA",
                             slot = "counts") + ggplot2::theme_classic()
    plot1 <- plot1 + ggplot2::scale_x_discrete(labels = c(
      paste0("naive - will not survive\n(",percentage1, "% of ", length(idx1), " cells)"),
      paste0("naive - will survive\n(",percentage2, "% of ", length(idx2), " cells)"),
      paste0(treatment, "\n(",percentage3, "% of ", length(idx3), " cells)")))
    plot1 <- plot1 + Seurat::NoLegend()
    ggplot2::ggsave(filename = paste0("../../out/figures/Writeup3e/Writeup3e_sydney_naive_within", treatment, "_vln_", marker_genes[i], "_lognorm_average_deathtoall.png"),
                    plot1, device = "png", width = 6, height = 5, units = "in")
  }
}


