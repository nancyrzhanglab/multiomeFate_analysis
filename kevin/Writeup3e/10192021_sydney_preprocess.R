rm(list=ls())
library(Seurat); library(Signac)
source("../Writeup3d/funcs.R")

load("../../../../data/Sydney_stressors_2021-09-24/all_data_SCT.RData")

# keep only lineages with more than 5 counts
lin_mat <- all_data[["lineage"]]@counts
counts_vec <- sparseMatrixStats::rowSums2(lin_mat)
lin_mat <- lin_mat[counts_vec >= 5, ]
depth_vec <- sparseMatrixStats::colSums2(lin_mat)
lin_mat <- lin_mat[,depth_vec > 0]

cell_idx <- sapply(1:ncol(lin_mat), function(i){
  if(i %% floor(ncol(lin_mat)/10) == 0) cat('*')
  
  vec <- lin_mat[,i]
  idx <- which(vec > 0)
  if(length(idx) == 0) return(NA)
  if(length(idx) == 1) return(idx[1])
  val <- sort(vec[idx], decreasing = T)
  if(val[1] >= 2*val[2]) return(idx[which.max(vec[idx])])
  return(NA)
})
names(cell_idx) <- colnames(lin_mat)

cell_name <- colnames(lin_mat)
cell_name <- cell_name[which(!is.na(cell_idx))]
cell_idx <- cell_idx[which(!is.na(cell_idx))]
lin_mat2 <- Matrix::sparseMatrix(i = cell_idx, 
                                 j = 1:length(cell_idx),
                                 dims = c(nrow(lin_mat), length(cell_idx)))
lin_mat2 <- as(lin_mat2, "dgCMatrix")
rownames(lin_mat2) <- rownames(lin_mat)
colnames(lin_mat2) <- cell_name

while(TRUE){
  bool1 <- TRUE
  counts_vec <- sparseMatrixStats::rowSums2(lin_mat2)
  lin_mat2 <- lin_mat2[counts_vec > 0, ]
  if(any(counts_vec == 0)) bool1 <- FALSE
  
  depth_vec <- sparseMatrixStats::colSums2(lin_mat2)
  lin_mat2 <- lin_mat2[,depth_vec > 0]
  if(any(depth_vec == 0)) bool1 <- FALSE
  
  if(bool1) break()
}

keep_vec <- rep(0, ncol(all_data))
keep_vec[which(colnames(all_data) %in% colnames(lin_mat2))] <- 1
table(keep_vec)

all_data[["keep"]] <- keep_vec
all_data <- subset(all_data, keep == 1)
all_data[["lineage"]]@counts <- lin_mat2
all_data[["lineage"]]@data <- lin_mat2

all_data2 <- all_data
while(TRUE){
  print(dim(all_data2[["lineage"]]@counts))
  bool1 <- TRUE
  
  lin_mat3 <- all_data2[["lineage"]]@counts
  tabulate_mat <- .tabulate_lineages(all_data2)
  
  # make sure lineages only have between 1 and 10 naive cells
  lineage_keep <- rownames(tabulate_mat)[intersect(
    which(tabulate_mat[,"naive"] > 0), 
    which(tabulate_mat[,"naive"] <= 10)
  )]
  if(length(lineage_keep) != nrow(lin_mat3)) bool1 <- FALSE
  lin_mat3 <- lin_mat3[lineage_keep,]
  
  # make sure each cell is assigned to a lineage
  depth_vec <- sparseMatrixStats::colSums2(lin_mat3)
  lin_mat3 <- lin_mat3[,depth_vec > 0]
  if(any(depth_vec == 0)) bool1 <- FALSE
 
  # make sure each lineage has a cell
  counts_vec <- sparseMatrixStats::rowSums2(lin_mat3)
  lin_mat3 <- lin_mat3[counts_vec > 0, ]
  if(any(counts_vec == 0)) bool1 <- FALSE
  
  keep_vec <- rep(0, ncol(all_data2))
  keep_vec[which(colnames(all_data2) %in% colnames(lin_mat3))] <- 1
  all_data2[["keep"]] <- keep_vec
  all_data2 <- subset(all_data2, keep == 1)
  all_data2[["lineage"]]@counts <- lin_mat3
  if(bool1) break()
}

all_data2[["lineage"]]@data <- all_data2[["lineage"]]@counts

### now do our checks ####
tabulate_mat <- .tabulate_lineages(all_data2)
head(tabulate_mat)
dim(tabulate_mat)
# make sure lineages only have between 1 and 10 naive cells
stopifnot(all(tabulate_mat[,"naive"] > 0),
          all(tabulate_mat[,"naive"] <= 10))
# make sure each cell is assigned to a lineage
stopifnot(all(sparseMatrixStats::colSums2(all_data2[["lineage"]]@counts) > 0))
# make sure each lineage has a cell
stopifnot(all(sparseMatrixStats::rowSums2(all_data2[["lineage"]]@counts) > 0))

dim(all_data2[["RNA"]])
dim(all_data2[["lineage"]])

all_data <- all_data2
save(all_data, 
     file = "../../../../out/kevin/Writeup3e/10192021_sydney_preprocess.RData")

head(all_data@meta.data)
gene_names <- rownames(all_data)
all_data <- Seurat::SCTransform(all_data,
                                residual.features = gene_names,
                                do.center = F,
                                verbose = T)
dim(all_data[["SCT"]]@scale.data)
save(all_data, 
     file = "../../../../out/kevin/Writeup3e/10192021_sydney_preprocess.RData")


###########################

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
    ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup3e/Writeup3e_sydney_naive_within", treatment, "_vln_", marker_genes[i], ".png"),
                    plot1, device = "png", width = 6, height = 5, units = "in")
  }
}

