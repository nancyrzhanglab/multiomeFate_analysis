rm(list=ls())
library(patchwork)
library(ggplot2)
library(ggtern)
library(Seurat)
library(multiomeFate)
load("~/project/Multiome_fate/out/kevin/Writeup8/Writeup8_larry-dataset_step3_fasttopics.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

seurat_object_safe <- seurat_object

treatment_vec <- as.character(sort(unique(seurat_object$time_celltype)))
day_early <- "4"
day_later <- "6"
day_early_vec <- treatment_vec[grep(paste0("^.*-", day_early), treatment_vec)]
treatment_vec <- treatment_vec[grep(paste0("^.*-", day_later), treatment_vec)]

#####

cell_imputation_mat <- numeric(0)

for(treatment in treatment_vec){
  load(paste0("~/project/Multiome_fate/out/kevin/Writeup9/Writeup9_", treatment, "_from_day", day_early, "_postprocess.RData"))
  
  cell_imputation_mat <- cbind(cell_imputation_mat, cell_imputed_score)
}
colnames(cell_imputation_mat) <- treatment_vec

cell_imputation_mat2 <- 10^cell_imputation_mat
cellsize <- rowSums(cell_imputation_mat2)
n <- nrow(cell_imputation_mat2)
for(i in 1:n){
  tmp <- cell_imputation_mat2[i,]
  if(sum(tmp) <= 0.01){
    cell_imputation_mat2[i,] <- NA
  } else {
    cell_imputation_mat2[i,] <- tmp/sum(tmp)
  }
}
colnames(cell_imputation_mat2) <- paste0("percent_fate_", treatment_vec)
idx <- unique(unlist(apply(cell_imputation_mat2, 2, function(x){which(is.na(x))})))
if(length(idx) > 0) {
  cell_imputation_mat2 <- cell_imputation_mat2[-idx,,drop = FALSE]
  cellsize <- cellsize[-idx]
}

celltype_vec <- as.character(seurat_object$time_celltype[rownames(cell_imputation_mat2)])
celltype_vec <- sapply(celltype_vec, function(val){
  strsplit(val, split = "-")[[1]][1]
})
names(celltype_vec) <- rownames(cell_imputation_mat2)

cell_imputation_mat3 <- cell_imputation_mat2
set.seed(10)
n <- nrow(cell_imputation_mat3)
for(i in 1:n){
  if(any(is.na(cell_imputation_mat3[i,]))) next()
  cell_imputation_mat3[i,] <- cell_imputation_mat3[i,] + stats::runif(3, min = 0, max = 0.1)
  cell_imputation_mat3[i,] <- cell_imputation_mat3[i,]/sum(cell_imputation_mat3[i,])
}

# https://cran.r-project.org/web/packages/Ternary/vignettes/Ternary.html ?
# https://www.marvinschmitt.com/blog/ggsimplex-prerelease/
# https://cran.r-project.org/web/packages/ggtern/index.html
# https://rpubs.com/KDVdecisions/triadtutorial1
# http://www.ggtern.com/d/2.2.2/ggsave.html

df <- as.data.frame(cell_imputation_mat3)
colnames(df) <- c("Monocyte", "Neutrophil", "Undifferentiated")
df <- cbind(df, 
            celltype_vec[rownames(cell_imputation_mat3)], 
            cellsize[rownames(cell_imputation_mat3)])
rownames(df) <- rownames(cell_imputation_mat3)
colnames(df)[4:5] <- c("celltype", "size")
df$celltype <- factor(df$celltype)
idx <- which(is.na(df[,1]))
if(length(idx) > 0) df <- df[-idx,]

color_palette <- c("blue3", "coral2", "gray50")
names(color_palette) <- paste0(c("Monocyte", "Neutrophil", "Undifferentiated"))
tab_vec <- table(df$celltype)
tab_vec <- sort(tab_vec, decreasing = TRUE)
row_idx <- unlist(lapply(names(tab_vec), function(val){
  which(df$celltype == val)
}))
df <- df[row_idx,]

plot1 <- ggtern::ggtern(data = df,
                        mapping = ggplot2::aes(x = Monocyte, 
                                               y = Neutrophil, 
                                               z = Undifferentiated,
                                               color = celltype,
                                               size = size)) +
  ggplot2::geom_point() +
  ggplot2::scale_color_manual(values = color_palette) +
  ggtern::theme_showarrows() + 
  ggplot2::scale_size_area() + 
  ggplot2::labs(x = "Monocyte", 
                y = "Neutrophil",
                z = "Undifferentiated",
                title = paste0("Percent fate for day ", day_early))
ggtern::ggsave(filename = paste0("~/project/Multiome_fate/out/figures/Writeup9/Writeup9_percent-fate_", day_early, "_cell-pairs.png"),
               plot1, 
               device = "png", 
               width = 10, 
               height = 5, 
               units = "in")