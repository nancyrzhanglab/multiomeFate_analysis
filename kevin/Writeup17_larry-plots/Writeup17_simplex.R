rm(list=ls())
library(Seurat)
library(ggplot2)
library(reshape2)
library(multiomeFate)
library(ggtern)

out_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/out/Writeup13/"
plot_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/git/multiomeFate_analysis/fig/kevin/Writeup17/"
code_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/git/multiomeFate_analysis/kevin/Writeup17_larry-plots/"

load(paste0(out_folder, "Writeup13_larry-dataset_step2_fasttopics.RData"))

# construct the matrix of fate potentials

color_by <- "celltype"
size_by <- "cellsize"

time_early <- 4
time_later <- 6

treatment_vec <- as.character(sort(unique(seurat_object$time_celltype)))
celltype_vec <- c("Monocyte", "Neutrophil", "Undifferentiated")
day_early_vec <- treatment_vec[grep("^.*-4", treatment_vec)]
day_early <- "4"
day_later <- "6"
treatment_vec <- treatment_vec[grep("^.*-6", treatment_vec)]

cell_imputation_mat <- numeric(0)
for(treatment in treatment_vec){
  load(paste0(out_folder, "Writeup13_", treatment, "_from_day", day_early, "_lineage-imputation.RData"))
  
  cell_imputation_mat <- cbind(cell_imputation_mat, cell_imputed_score)
}
colnames(cell_imputation_mat) <- celltype_vec

stopifnot(all(!is.na(cell_imputation_mat)))

df <- multiomeFate:::compute_entropy(cell_imputation_mat = cell_imputation_mat,
                                     later_timepoint = time_later,
                                     seurat_object = seurat_object,
                                     variable_celltype = "Cell.type.annotation",
                                     variable_lineage = "assigned_lineage",
                                     variable_timepoint = "Time.point")
color_palette <- c("#A6CEE3", "#FCA5A5", "#A9A9A9")
names(color_palette) <- paste0(c("Monocyte", "Neutrophil", "Undifferentiated"))

color_by_name <- ifelse(color_by == "celltype", "celltype", "dominantFate")
color_by_plot <- ifelse(color_by == "celltype", "current cell type", "future dominant fate")
size_by_plot <- ifelse(size_by == "cellsize", "# predicted progenies", "entropy")

aes_formula <- ggplot2::aes_string(x = "Monocyte", 
                                   y = "Neutrophil", 
                                   z = "Undifferentiated", 
                                   color = color_by, 
                                   size = size_by)
# rearrange by color
tmp <- lapply(rev(names(color_palette)), function(celltype){
  which(df[,color_by] == celltype)
})
df <- df[unlist(tmp),]

plot1 <- multiomeFate:::plot_simplex(
  df = df,
  aes_formula = aes_formula, 
  col_palette = color_palette,
  xlab = "Monocyte",
  ylab = "Neutrophil",
  zlab = "Undifferentiated",
  title = paste0("Day ", time_early, " cells predicting ", time_later, ", color by ", color_by_plot, "\nSize by ", size_by_plot)
)

# create empty plot
png(file = paste0(plot_folder, "Writeup17_d46_simplex_color-", color_by_name, "_size-", size_by, ".png"),
    width = 3000, 
    height = 2000, 
    res = 300, 
    units = "px")
print(plot1)
graphics.off()

##########

df$celltype <- factor(df$celltype, levels = c("Monocyte", "Neutrophil", "Undifferentiated"))

plot1 <- ggtern::ggtern(df, aes(x = Monocyte, y = Neutrophil, z = Undifferentiated, 
                                color = celltype)) +
  stat_density_tern(aes(alpha = ..level.., fill = celltype), 
                    geom = 'polygon', 
                    bins = 20) +
  ggplot2::scale_color_manual(values = color_palette) + 
  ggplot2::scale_fill_manual(values = color_palette) + Seurat::NoLegend() +
  ggplot2::labs(x = "", 
                y = "",
                z = "",
                title = "")

png(file = paste0(plot_folder, "Writeup17_d46_simplex_color-", color_by_name, "_size-", size_by, "_cleaned.png"),
    width = 1500, 
    height = 1000, 
    res = 300, 
    units = "px")
print(plot1)
graphics.off()


