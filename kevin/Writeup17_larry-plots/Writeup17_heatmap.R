rm(list=ls())
library(Seurat)
library(ggplot2)
library(reshape2)

out_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/out/Writeup13/"
plot_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/git/multiomeFate_analysis/fig/kevin/Writeup17/"
code_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/git/multiomeFate_analysis/kevin/Writeup17_larry-plots/"

load(paste0(out_folder, "Writeup13_larry-dataset.RData"))
source(paste0(code_folder, "welch_anova.R"))

seurat_object_full <- seurat_object
seurat_object <- subset(seurat_object, Time.point == "4")

# only keep the top 10 lineages
tab_vec <- table(seurat_object$assigned_lineage)
tab_vec <- sort(tab_vec, decreasing = TRUE)
keep_lineages <- names(tab_vec)[1:10]

seurat_object <- subset(seurat_object, assigned_lineage %in% keep_lineages)
mat <- SeuratObject::LayerData(
  seurat_object,
  layer = "data",
  assay = "RNA",
  features = Seurat::VariableFeatures(seurat_object)
)
mat@x <- pmin(mat@x, 20)

rowsum_vec <- Matrix::rowSums(mat)
if(any(rowsum_vec <= 1)){
  mat <- mat[-which(rowsum_vec <= 1),]
}

anova_res <- sapply(1:nrow(mat), function(j){
  if(j %% floor(nrow(mat)/10) == 0) cat('*')
  
  df <- data.frame(
    response = mat[j,],
    group = seurat_object$assigned_lineage
  )
  
  df <- .remove_zero_variance(
    data = df, 
    response_col = "response", 
    group_col = "group"
  )
  
  if(length(unique(df$group)) <= 1) return(c(pvalue = NA, R2 = NA))
  res <- .welch_anova(
    data = df, 
    response_col = "response", 
    group_col = "group"
  )
  
  c(pvalue = res$p.value, R2 = res$R2_welch)
})
anova_res <- t(anova_res)
rownames(anova_res) <- rownames(mat)
anova_res <- anova_res[which(!is.na(anova_res[,"pvalue"])),]

anova_res[,"pvalue"] <- pmax(anova_res[,"pvalue"], 1e-15)
plot(-log10(anova_res[,"pvalue"]), anova_res[,"R2"], pch = 16) # double checking

anova_res <- anova_res[order(anova_res[,"R2"], decreasing = TRUE),]
selected_genes <- rownames(anova_res)[1:10]
mat_subset <- mat[selected_genes,]

###########

library(ComplexHeatmap)
library(circlize)
library(Matrix)

# Convert to a dense matrix for visualization
mat_dense <- as.matrix(mat_subset)

# Transpose the matrix so that cells (original columns) are now rows
mat_transposed <- t(mat_dense)

# Extract lineage labels for each cell (row in transposed matrix)
lineages <- seurat_object$assigned_lineage

# Ensure the ordering of the rows based on lineage size (descending)
lineage_counts <- table(lineages)
sorted_lineages <- names(sort(lineage_counts, decreasing = TRUE))
ordered_indices <- order(match(lineages, sorted_lineages))
mat_transposed <- mat_transposed[ordered_indices, ]
lineages <- lineages[ordered_indices]

# Define a color mapping function
col_fun <- colorRamp2(
  c(0, min(mat_transposed[mat_transposed > 0], na.rm = TRUE), max(mat_transposed, na.rm = TRUE)),
  c("white", "lightgray", "forestgreen")
)

# Create row annotation to separate lineages
row_annotation <- factor(lineages, levels = sorted_lineages)
row_split <- factor(lineages, levels = sorted_lineages)

# Plot the heatmap
png(filename = paste0(plot_folder, "heatmap.png"),
    height = 5000, width = 2000, units = "px", res = 300)
# Plot the heatmap without the legend
Heatmap(
  mat_transposed,
  name = "Expression",
  col = col_fun,
  cluster_rows = FALSE,  # Keep lineage order
  cluster_columns = FALSE,  
  show_row_dend = FALSE,  # Hide row dendrogram
  show_column_dend = FALSE,  # Hide column dendrogram
  show_row_names = FALSE,  # REMOVE row names (cell names)
  row_split = row_split,  # Group by lineage
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_gap = unit(1, "mm"),  # Ensure small space between individual rows
  border = TRUE,  # Add black lines between lineage groups
  border_gp = gpar(lwd = 5),
  show_heatmap_legend = FALSE,  # REMOVE the legend
  column_names_gp = gpar(fontsize = 8),
  row_title_rot = 0,
  row_title_side = "left"
)
dev.off()

# Plot the heatmap
png(filename = paste0(plot_folder, "heatmap_cleaned.png"),
    height = 1500, width = 800, units = "px", res = 300)
# Plot the heatmap without the legend
Heatmap(
  mat_transposed,
  name = "Expression",
  col = col_fun,
  cluster_rows = FALSE,  # Keep lineage order
  cluster_columns = FALSE,  
  show_row_dend = FALSE,  # Hide row dendrogram
  show_column_dend = FALSE,  # Hide column dendrogram
  show_row_names = FALSE,  # REMOVE row names (cell names)
  row_split = row_split,  # Group by lineage
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_gap = unit(1, "mm"),  # Ensure small space between individual rows
  border = TRUE,  # Add black lines between lineage groups
  border_gp = gpar(lwd = 3),
  show_heatmap_legend = FALSE,  # REMOVE the legend
  column_names_gp = gpar(fontsize = 10),
  row_title_rot = 0,
  row_title_side = "left"
)
dev.off()

########################

lineage_names <- levels(row_split)
tab_mat <- table(seurat_object_full$assigned_lineage, seurat_object_full$time_celltype)
lineage_size_mat <- tab_mat[lineage_names, c("Monocyte-6", "Neutrophil-6", "Undifferentiated-6")]

# Convert the matrix to a tidy data frame
df <- as.data.frame(as.table(as.matrix(lineage_size_mat)))

# Rename columns for clarity
colnames(df) <- c("Lineage", "Cell_Type", "Value")
df$Lineage <- factor(df$Lineage, levels = rev(levels(df$Lineage)))

# Define color mapping for each cell type
color_map <- c("Monocyte-6" = "#4472C4", 
               "Neutrophil-6" = "#D55E00" , 
               "Undifferentiated-6" = "#7F7F7F")

# Create the dot plot
plot1 <- ggplot(df, aes(x = Cell_Type, y = Lineage, size = Value, color = Cell_Type)) +
  geom_point() +
  scale_size(range = c(0, 8)) +  # Adjust dot size (0 for no values)
  scale_color_manual(values = color_map) +  # Set colors per column
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    axis.title = element_blank(),
    legend.position = "right"
  ) +
  guides(color = "none")  # Hide color legend (since it's redundant)
ggplot2::ggsave(plot1, filename = paste0(plot_folder, "day6_lineage_size.png"),
                height = 8, width = 8)

plot1 <- plot1 + theme_minimal() +
  theme(
    axis.text.x = element_blank(),  # Remove x-axis labels
    axis.text.y = element_blank(),  # Remove y-axis labels
    axis.title = element_blank(),  # Remove axis titles
    legend.position = "right"
  ) +
  guides(color = "none") + Seurat::NoLegend()
ggplot2::ggsave(plot1, filename = paste0(plot_folder, "day6_lineage_size_cleaned.png"),
                height = 4, width = 3)

