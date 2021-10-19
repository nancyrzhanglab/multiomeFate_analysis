rm(list=ls())
library(Seurat)
load("../../../../out/kevin/Writeup3d/09302021_sydney_basic_preprocess.RData")
gene_names <- rownames(all_data)
all_data <- Seurat::SCTransform(all_data,
                                residual.features = gene_names,
                                do.center = F,
                                verbose = T)
dim(all_data[["SCT"]]@scale.data)

save(all_data, 
     file = "../../../../out/kevin/Writeup3e/10122021_sydney_sctransform_preprocess.RData")

###############
load("../../../../out/kevin/Writeup3e/10122021_sydney_sctransform_preprocess.RData")
source("../Writeup3d/funcs.R")

Seurat::Idents(all_data) <- "Original_condition"
table(Seurat::Idents(all_data))

tabulate_mat <- .tabulate_lineages(all_data)
head(tabulate_mat)
quantile(tabulate_mat[,"naive"], probs = seq(0.9,1,length.out=11))
max_val <- sapply(1:nrow(tabulate_mat), function(i){max(tabulate_mat[i,-5])})
tabulate_mat <- tabulate_mat[order(max_val, decreasing = T),]

##############################

# first make sure we can separate cocl2 from naive
marker_genes <- c("AXL", "EGFR", "JUN", "VEGFC", "WNT5A", "NGFR",
                  "SERPINE1", "FGFR1", "LOXL2", "NRG1", "PDGFRB")
set.seed(10)
seurat_de <- Seurat::FindMarkers(all_data,
                                 assay = "SCT",
                                 slot = "scale.data",
                                 ident.1 = "cocl2",
                                 ident.2 = "naive",
                                 logfc.threshold = 0.1,
                                 min.pct = 0.05,
                                 verbose = T)
seurat_de[c("AXL", "NRG1", "EGFR"),]
seurat_de[1:50,]
head(seurat_de)

for(i in 1:length(marker_genes)){
  plot1 <- Seurat::VlnPlot(all_data, features = marker_genes[i],
                           idents = c("cocl2", "naive"),
                           assay = "SCT",
                           slot = "scale.data") + ggplot2::theme_classic()
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup3e/Writeup3e_sydney_Cocl2vsNaive_vln_", marker_genes[i], ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}

mean(all_data[["SCT"]]@scale.data["EGFR",])

###########################3

# first up -- cocl2
treatment_vec <- c("cocl2", "acid", "dab", "tram")
pval_thres <- c(1e-3, 1e-3, 1e-3, 1e-2)
de_list <- vector("list", 4)

table(all_data@meta.data$Original_condition)
for(kk in 1:length(treatment_vec)){
  print(treatment_vec[kk])
  treatment <- treatment_vec[kk]
  
  # set up identities
  res <- .select_expansion_naives(all_data, 
                                  tabulate_mat, 
                                  treatment = treatment)
  newgroup <- all_data@meta.data$Original_condition
  names(newgroup) <- rownames(all_data@meta.data)
  newgroup[res$naive_terminal] <- paste0("naive_survive", treatment)
  newgroup[setdiff(res$naive_all, res$naive_terminal)] <- paste0("naive_no", treatment)
  all_data[[paste0(treatment, "Status")]] <- newgroup
  Seurat::Idents(all_data) <- paste0(treatment, "Status")

  set.seed(10)
  seurat_de <- Seurat::FindMarkers(all_data,
                                   assay = "SCT",
                                   slot = "scale.data",
                                   ident.1 = paste0("naive_survive", treatment),
                                   ident.2 = paste0("naive_no", treatment),
                                   logfc.threshold = 0.1,
                                   min.pct = 0.05,
                                   verbose = T)
  de_list[[kk]] <- seurat_de
  
  tentative_genes <- rownames(seurat_de)[which(seurat_de$p_val_adj <= pval_thres[kk])]
  proportion_list <- .order_genes_by_threshold(genes = tentative_genes,
                                               naive_terminal = res$naive_terminal,
                                               seurat_obj = all_data)
  
  marker_genes <- unique(c(names(proportion_list)[1:min(length(proportion_list), 15)], "AXL", "EGFR", "NRG1"))
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
      paste0("naive - will not survive\n(",percentage1, ")"),
      paste0("naive - will survive\n(",percentage2, ")"),
      paste0(treatment, "\n(",percentage3, ")")))
    plot1 <- plot1 + Seurat::NoLegend()
    ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup3e/Writeup3e_sydney_naive_within", treatment, "_vln_", marker_genes[i], ".png"),
                    plot1, device = "png", width = 6, height = 5, units = "in")
  }
}

#############################

tabulate_mat2 <- tabulate_mat
tabulate_mat2 <- tabulate_mat2[which(tabulate_mat2[,"naive"] > 2),]
# tabulate_mat2 <- tabulate_mat2[,which(colnames(tabulate_mat2) != "naive")]
tabulate_mat2 <- tabulate_mat2[which(apply(tabulate_mat2, 1, max) >= 100),]
tabulate_mat2 <- log1p(tabulate_mat2)
tabulate_mat2 <- as.data.frame(tabulate_mat2)

plot1 <- GGally::ggpairs(tabulate_mat2)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup3e/Writeup3e_sydney_pairs.png"),
                plot1, device = "png", width = 7, height = 7, units = "in")

#############################


