cell_imputed_score <- all_data@meta.data[,paste0("fatepotential_", file)]
names(cell_imputed_score) <- Seurat::Cells(all_data)

bool_anova = TRUE
bool_mark_mean = TRUE
bool_mark_max = FALSE
min_lineage_size = 2
num_lineages_top = 10
num_lineages_bottom = 10

seurat_object = all_data
cell_imputed_score = cell_imputed_score
assigned_lineage_variable = "assigned_lineage"
time_celltype_variable = "dataset"
day_later = day_later_full
bool_anova = FALSE
ylab = paste0(day_later_vec, " growth potential")
ylim = c(max(stats::quantile(cell_imputed_score, na.rm = TRUE, prob = 0.05), -5), 
         NA)

# grab the vector of which celltype-time each cell is
assigned_lineage <- seurat_object@meta.data[,assigned_lineage_variable]
names(assigned_lineage) <- Seurat::Cells(seurat_object)

time_celltype <- seurat_object@meta.data[,time_celltype_variable]
names(time_celltype) <- Seurat::Cells(seurat_object)
stopifnot(day_later %in% time_celltype)

# determine which lineages qualify to be in the plot
lineage_vec <- assigned_lineage[names(cell_imputed_score)]
tab_mat <- table(assigned_lineage, time_celltype)
lineage_future_size <- tab_mat[, day_later]
names(lineage_future_size) <- rownames(tab_mat)

##########

seurat_object = seurat_object
cell_imputed_score = cell_imputed_score
assigned_lineage_variable = assigned_lineage_variable
lineage_future_size = lineage_future_size
bool_anova = bool_anova
bool_mark_mean = bool_mark_mean
bool_mark_max = bool_mark_max
min_lineage_size = min_lineage_size
num_lineages_top = num_lineages_top
num_lineages_bottom = num_lineages_bottom
ylab = ylab
ylim = ylim

col_all_lineages = "#E69F00"
col_lineages = "#999999"


stopifnot(length(names(cell_imputed_score)) == length(cell_imputed_score))

if(any(is.na(cell_imputed_score))){
  cell_imputed_score <- cell_imputed_score[!is.na(cell_imputed_score)]
}

# grab the vector of which celltype-time each cell is
assigned_lineage <- seurat_object@meta.data[,assigned_lineage_variable]
names(assigned_lineage) <- Seurat::Cells(seurat_object)
assigned_lineage <- assigned_lineage[names(cell_imputed_score)]

# filter out lineages that too small
lineage_vec <- assigned_lineage[names(cell_imputed_score)]
tab_vec <- table(assigned_lineage)
tab_vec <- tab_vec[tab_vec >= min_lineage_size] # current size needs to be big enough
passed_lineages <- names(tab_vec)
passed_cells_names <- names(assigned_lineage)[which(assigned_lineage %in% passed_lineages)]
cell_imputed_score <- cell_imputed_score[passed_cells_names]
lineage_vec <- assigned_lineage[names(cell_imputed_score)]
lineage_future_size <- lineage_future_size[which(names(lineage_future_size) %in% passed_lineages)]

# determine which lineages qualify to be in the plot
lineage_names_ordered <- names(lineage_future_size)[order(lineage_future_size, decreasing = TRUE)]
lineage_names_top <- lineage_names_ordered[1:num_lineages_top]
lineage_names_bottom <- lineage_names_ordered[(length(lineage_names_ordered)-num_lineages_bottom+1):length(lineage_names_ordered)]
lineage_names <- unique(c(lineage_names_top, lineage_names_bottom))
idx <- which(lineage_vec %in% lineage_names)
