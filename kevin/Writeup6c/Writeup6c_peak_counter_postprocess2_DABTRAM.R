rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")
load("../../../../out/kevin/Writeup6c/Writeup6c_peak_counter_5000.RData")

countmat_nopeak <- Matrix::Matrix(countmat_nopeak, sparse = T)
quantile(countmat_nopeak@x)

# format the Seurat object appropriately
all_data[["offpeak"]] <- Seurat::CreateAssayObject(counts = Matrix::t(countmat_nopeak))

#########################

lin_mat <- table(all_data$assigned_lineage, all_data$dataset)
DABTRAM_lineage <- rownames(lin_mat)[which(lin_mat[,"day10_DABTRAM"] >= 20)]
day0_survive_idx <- intersect(which(all_data$assigned_lineage %in% DABTRAM_lineage),
                              which(all_data$dataset == "day0"))
day0_die_idx <- setdiff(intersect(which(!is.na(all_data$assigned_lineage)),
                                  which(all_data$dataset == "day0")),
                        day0_survive_idx)
day0_unknown_idx <- setdiff(which(all_data$dataset == "day0"),
                            c(day0_survive_idx, day0_die_idx))

tmp_vec <- rep(NA, ncol(all_data))
tmp_vec[day0_survive_idx] <- "day0_survive"
tmp_vec[day0_die_idx] <- "day0_die"
tmp_vec[day0_unknown_idx] <- "day0_unknown"
all_data$day0_DABTRAM_status <- tmp_vec

cell_idx <- sort(unique(c(day0_survive_idx, day0_die_idx, day0_unknown_idx)))
cell_idx2 <- sort(unique(c(day0_survive_idx, day0_die_idx)))
lineage_vec <- all_data$assigned_lineage[cell_idx2]
fitness_vec_onlyday0 <- sapply(1:length(cell_idx2), function(i){
  log1p(lin_mat[lineage_vec[i],"day10_DABTRAM"])
})
fitness_vec <- rep(NA, ncol(all_data))
fitness_vec[cell_idx2] <- fitness_vec_onlyday0
all_data$day0_DABTRAM_fitness <- fitness_vec

n <- ncol(all_data)

##########################

.nonzero_col <- function(mat, col_idx, bool_value){
  stopifnot(inherits(mat, "dgCMatrix"), col_idx %% 1 == 0,
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

cell_idx <- sort(unique(c(day0_survive_idx, day0_die_idx, day0_unknown_idx)))
survive_idx <- which(cell_idx %in% day0_survive_idx)
die_idx <- which(cell_idx %in% day0_die_idx)
unknown_idx <- which(cell_idx %in% day0_unknown_idx)

fitness_vec2 <- all_data$day0_DABTRAM_fitness[cell_idx2]
survive_idx2 <- which(cell_idx2 %in% day0_survive_idx)
die_idx2 <- which(cell_idx2 %in% day0_die_idx)

length(survive_idx) == length(survive_idx2)
length(die_idx) == length(die_idx2)

##########################

# 1) do number of active genes
active_gene_vec <- all_data$nFeature_RNA[cell_idx]
stats::cor(active_gene_vec[c(survive_idx, die_idx)], 
           fitness_vec2[c(survive_idx2, die_idx2)], 
           method = "spearman")
stats::wilcox.test(x = active_gene_vec[survive_idx],
                   y = active_gene_vec[die_idx])
stats::wilcox.test(x = active_gene_vec[survive_idx],
                   y = active_gene_vec[c(die_idx, unknown_idx)])

# pretty promising
active_gene_vec2 <- sapply(cell_idx, function(i){
  length(which(
    .nonzero_col(all_data[["Saver"]]@data, col_idx = i, bool_value = T) > 0.5
  ))
})
stats::cor(active_gene_vec2[c(survive_idx, die_idx)], 
           fitness_vec2[c(survive_idx2, die_idx2)], 
           method = "spearman")
stats::wilcox.test(x = active_gene_vec2[survive_idx],
                   y = active_gene_vec2[die_idx])
stats::wilcox.test(x = active_gene_vec2[survive_idx],
                   y = active_gene_vec2[c(die_idx, unknown_idx)])

##########################

# 2) do number of genes w/ off-peak counts
active_gene_vec <- sapply(cell_idx, function(i){
  length(
    .nonzero_col(all_data[["offpeak"]]@data, col_idx = i, bool_value = T)
  )
})
stats::cor(active_gene_vec[c(survive_idx, die_idx)], 
           fitness_vec2[c(survive_idx2, die_idx2)], 
           method = "spearman")
stats::wilcox.test(x = active_gene_vec[survive_idx],
                   y = active_gene_vec[die_idx])
stats::wilcox.test(x = active_gene_vec[survive_idx],
                   y = active_gene_vec[c(die_idx, unknown_idx)])

##########################

# 3) sum of off-peak counts
active_gene_vec <- sapply(cell_idx, function(i){
  sum(
    .nonzero_col(all_data[["offpeak"]]@counts, col_idx = i, bool_value = T)
  )
})
stats::cor(active_gene_vec[c(survive_idx, die_idx)], 
           fitness_vec2[c(survive_idx2, die_idx2)], 
           method = "spearman")
stats::wilcox.test(x = active_gene_vec[survive_idx],
                   y = active_gene_vec[die_idx])
stats::wilcox.test(x = active_gene_vec[survive_idx],
                   y = active_gene_vec[c(die_idx, unknown_idx)])

#############################

# 4) normalize the off-peak counts, and then sum
Seurat::DefaultAssay(all_data) <- "offpeak"
all_data <- Seurat::NormalizeData(all_data)

active_gene_vec <- sapply(cell_idx, function(i){
  sum(
    .nonzero_col(all_data[["offpeak"]]@data, col_idx = i, bool_value = T)
  )
})
stats::cor(active_gene_vec[c(survive_idx, die_idx)], 
           fitness_vec2[c(survive_idx2, die_idx2)], 
           method = "spearman")
stats::wilcox.test(x = active_gene_vec[survive_idx],
                   y = active_gene_vec[die_idx])
stats::wilcox.test(x = active_gene_vec[survive_idx],
                   y = active_gene_vec[c(die_idx, unknown_idx)])

# 5) TF-IDF the off-peak counts, and then sum
Seurat::DefaultAssay(all_data) <- "offpeak"
all_data <- Signac::RunTFIDF(all_data)

active_gene_vec <- sapply(cell_idx, function(i){
  sum(
    .nonzero_col(all_data[["offpeak"]]@data, col_idx = i, bool_value = T)
  )
})
stats::cor(active_gene_vec[c(survive_idx, die_idx)], 
           fitness_vec2[c(survive_idx2, die_idx2)], 
           method = "spearman")
stats::wilcox.test(x = active_gene_vec[survive_idx],
                   y = active_gene_vec[die_idx])
stats::wilcox.test(x = active_gene_vec[survive_idx],
                   y = active_gene_vec[c(die_idx, unknown_idx)])

#############################

lineage_fitness_vec <- sapply(lineage_names, function(lineage_name){
  mean(all_data$day0_DABTRAM_fitness[cell_idx[which(lineage_vec == lineage_name)]])
})
all(!is.na(lineage_fitness_vec))

# 6) Let's try summing into averaging across lineages first, and then testing?
active_gene_vec <- sapply(cell_idx, function(i){
  sum(
    .nonzero_col(all_data[["offpeak"]]@counts, col_idx = i, bool_value = T)
  )
})
lineage_vec <- all_data$assigned_lineage[cell_idx]
lineage_names <- sort(unique(lineage_vec[!is.na(lineage_vec)]))
lineage_active_gene_vec <- sapply(lineage_names, function(lineage_name){
  mean(active_gene_vec[which(lineage_vec == lineage_name)])
})

stats::cor(lineage_active_gene_vec, lineage_fitness_vec)
lineage_thrive <- rownames(lin_mat)[which(lin_mat[,"day10_DABTRAM"] >= 20)]
stats::wilcox.test(x = lineage_active_gene_vec[which(names(lineage_active_gene_vec) %in% lineage_thrive)],
                   y = lineage_active_gene_vec[which(!names(lineage_active_gene_vec) %in% lineage_thrive)])

# 7) Let's try summing the normalized-peaks into averaging across lineages first, and then testing?
Seurat::DefaultAssay(all_data) <- "offpeak"
all_data <- Seurat::NormalizeData(all_data)

active_gene_vec <- sapply(cell_idx, function(i){
  sum(
    .nonzero_col(all_data[["offpeak"]]@data, col_idx = i, bool_value = T)
  )
})
lineage_vec <- all_data$assigned_lineage[cell_idx]
lineage_names <- sort(unique(lineage_vec[!is.na(lineage_vec)]))
lineage_active_gene_vec <- sapply(lineage_names, function(lineage_name){
  mean(active_gene_vec[which(lineage_vec == lineage_name)])
})

stats::cor(lineage_active_gene_vec, lineage_fitness_vec)
lineage_thrive <- rownames(lin_mat)[which(lin_mat[,"day10_DABTRAM"] >= 20)]
stats::wilcox.test(x = lineage_active_gene_vec[which(names(lineage_active_gene_vec) %in% lineage_thrive)],
                   y = lineage_active_gene_vec[which(!names(lineage_active_gene_vec) %in% lineage_thrive)])


# 7) Let's try summing the TFIDF-peaks into averaging across lineages first, and then testing?
Seurat::DefaultAssay(all_data) <- "offpeak"
all_data <- Signac::RunTFIDF(all_data)

active_gene_vec <- sapply(cell_idx, function(i){
  sum(
    .nonzero_col(all_data[["offpeak"]]@data, col_idx = i, bool_value = T)
  )
})
lineage_vec <- all_data$assigned_lineage[cell_idx]
lineage_names <- sort(unique(lineage_vec[!is.na(lineage_vec)]))
lineage_active_gene_vec <- sapply(lineage_names, function(lineage_name){
  mean(active_gene_vec[which(lineage_vec == lineage_name)])
})

stats::cor(lineage_active_gene_vec, lineage_fitness_vec)
lineage_thrive <- rownames(lin_mat)[which(lin_mat[,"day10_DABTRAM"] >= 20)]
stats::wilcox.test(x = lineage_active_gene_vec[which(names(lineage_active_gene_vec) %in% lineage_thrive)],
                   y = lineage_active_gene_vec[which(!names(lineage_active_gene_vec) %in% lineage_thrive)])

