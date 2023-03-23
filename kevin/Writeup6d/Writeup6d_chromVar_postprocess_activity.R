rm(list=ls())
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)

load("../../../../out/kevin/Writeup6b/Writeup6b_chromVar.RData")
source("../Writeup6b/gene_list.R")

Seurat::DefaultAssay(all_data) <- "chromvar"

treatment_vec <- c("CIS", "COCL2", "DABTRAM")
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
gene_vec <- sort(unique(unlist(keygenes)))
gene_vec <- gene_vec[which(gene_vec %in% rownames(all_data[["RNA"]]))]

for(treatment in treatment_vec){
  print("====")
  print(treatment)
  
  ident_vec <- all_data$dataset
  names(ident_vec) <- colnames(all_data)
  lineage_names <- rownames(tab_mat)[which(tab_mat[,paste0("day10_", treatment)] >= 20)]
  cell_names1 <- colnames(all_data)[which(all_data$assigned_lineage %in% lineage_names)]
  cell_names2 <- colnames(all_data)[which(all_data$dataset == "day0")]
  cell_names_winning <- intersect(cell_names1, cell_names2)
  cell_names_losing <- setdiff(cell_names2, cell_names1)
  ident_vec[cell_names_winning] <- paste0("day0_win_", treatment)
  ident_vec[cell_names_losing] <- paste0("day0_lose_", treatment)
  all_data$ident <- ident_vec
  Seurat::Idents(all_data) <- "ident"
  
  de_res <- Seurat::FindMarkers(
    object = all_data,
    ident.1 = paste0("day0_win_", treatment),
    ident.2 = paste0("day0_lose_", treatment),
    only.pos = FALSE,
    min.pct = 0,
    logfc.threshold = 0,
    mean.fxn = rowMeans,
    fc.name = "avg_diff",
    verbose = F
  )
  
  idx <- which(de_res[,"p_val_adj"] <= 1e-2)
  if(length(idx) == 0) idx <- 1
  print(length(idx))
  de_res[idx,]
  
  # hacked from https://github.com/stuart-lab/signac/blob/master/R/visualization.R
  motifs <- rownames(de_res)[idx]
  data.use <- Signac::GetMotifData(object = all_data,
                                   assay = "ATAC",
                                   slot = "pwm")
  data.use <- data.use[motifs]
  names(data.use) <- Signac::GetMotifData(
    object = all_data,
    assay = "ATAC",
    slot = "motif.names"
  )[motifs]
  names(data.use) <- sapply(1:length(names(data.use)), function(i){
    paste0(names(data.use)[i], ", ", rownames(de_res)[idx][i], 
           "\npadj=", formatC(de_res[idx[i],"p_val_adj"], format = "e", digits = 2),
           ", value=", round(de_res[idx[i],"avg_diff"], 2))
  })
  
  if(length(data.use) >= 50){
    data.use <- data.use[1:50]
  }
  
  plot1 <- ggseqlogo::ggseqlogo(data = data.use, ncol = 4)
  plot1 <- plot1 + ggplot2::theme_bw()
  
  if(length(data.use) <= 10){
    width <- 10; height <- 8
  } else if(length(data.use) <= 30){
    width <- 15; height <- 10
  } else {
    width <- 15; height <- 20
  }
  
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6d/Writeup6d_motif_day0-day10_", treatment, ".png"),
                  plot1, device = "png", width = width, height = height, units = "in")
}

#####################

# now do this for Day10->Week5

for(treatment in treatment_vec){
  print("====")
  print(treatment)
  
  ident_vec <- all_data$dataset
  names(ident_vec) <- colnames(all_data)
  lineage_names <- rownames(tab_mat)[which(tab_mat[,paste0("week5_", treatment)] >= 50)]
  cell_names1 <- colnames(all_data)[which(all_data$assigned_lineage %in% lineage_names)]
  cell_names2 <- colnames(all_data)[which(all_data$dataset == paste0("day10_", treatment))]
  cell_names_winning <- intersect(cell_names1, cell_names2)
  cell_names_losing <- setdiff(cell_names2, cell_names1)
  ident_vec[cell_names_winning] <- paste0("day10_win_", treatment)
  ident_vec[cell_names_losing] <- paste0("day10_lose_", treatment)
  all_data$ident <- ident_vec
  Seurat::Idents(all_data) <- "ident"
  
  de_res <- Seurat::FindMarkers(
    object = all_data,
    ident.1 = paste0("day10_win_", treatment),
    ident.2 = paste0("day10_lose_", treatment),
    only.pos = FALSE,
    min.pct = 0,
    logfc.threshold = 0,
    mean.fxn = rowMeans,
    fc.name = "avg_diff",
    verbose = F
  )
  
  idx <- which(de_res[,"p_val_adj"] <= 1e-2)
  if(length(idx) == 0) idx <- 1
  print(length(idx))
  
  # hacked from https://github.com/stuart-lab/signac/blob/master/R/visualization.R
  motifs <- rownames(de_res)[idx]
  data.use <- Signac::GetMotifData(object = all_data,
                                   assay = "ATAC",
                                   slot = "pwm")
  data.use <- data.use[motifs]
  names(data.use) <- Signac::GetMotifData(
    object = all_data,
    assay = "ATAC",
    slot = "motif.names"
  )[motifs]
  names(data.use) <- sapply(1:length(names(data.use)), function(i){
    paste0(names(data.use)[i], ", ", rownames(de_res)[idx][i], 
           "\npadj=", formatC(de_res[idx[i],"p_val_adj"], format = "e", digits = 2),
           ", value=", round(de_res[idx[i],"avg_diff"], 2))
  })
  
  if(length(data.use) >= 50){
    data.use <- data.use[1:50]
  }
  
  plot1 <- ggseqlogo::ggseqlogo(data = data.use, ncol = 4)
  plot1 <- plot1 + ggplot2::theme_bw()
  
  if(length(data.use) <= 10){
    width <- 10; height <- 8
  } else if(length(data.use) <= 30){
    width <- 15; height <- 10
  } else {
    width <- 15; height <- 20
  }
  
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6d/Writeup6d_motif_day10-week5_", treatment, ".png"),
                  plot1, device = "png", width = width, height = height, units = "in")
}

#######################

.nonzero_col <- function(mat, col_idx, bool_value){
  stopifnot(inherits(mat, c("dgCMatrix", "lgCMatrix")), col_idx %% 1 == 0,
            col_idx > 0, col_idx <= ncol(mat))
  
  val1 <- mat@p[col_idx]
  val2 <- mat@p[col_idx+1]
  
  if(val1 == val2) return(numeric(0))
  if(bool_value){
    # return the value
    mat@x[(val1+1):val2]
  } else {
    # return the row index
    mat@i[(val1+1):val2]+1
  }
}

all_data[["chromvar"]]@data[1:5,1:5]
all_data[["ATAC"]]@motifs@data[1:5,1:5]
all_data[["ATAC"]]@motifs@pwm[[1]]
head(all_data[["ATAC"]]@motifs@motif.names)
all_data[["ATAC"]]@motifs@positions[[1]]

zz <- all_data[["ATAC"]]@motifs@data
idx <- .nonzero_col(zz, col_idx = 1, bool_value = F)
length(idx)
rownames(all_data[["ATAC"]])[idx[1:5]]
length(all_data[["ATAC"]]@motifs@positions[[1]]) # hm... not sure
732642-732629+1
156009225-156009212+1
dim(all_data[["ATAC"]]@motifs@pwm[[1]])


