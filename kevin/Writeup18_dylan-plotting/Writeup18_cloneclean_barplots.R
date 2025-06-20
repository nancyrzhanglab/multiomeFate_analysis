rm(list=ls())
library(Seurat)
library(Signac)
library(multiomeFate)
library(clue)

out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup10a/"
load(paste0(out_folder, "Writeup10a_ppStep4_lineage.RData"))
plot_folder <- "~/project/Multiome_fate/git/multiomeFate_analysis_kevin/fig/kevin/Writeup18/"


lin_mat <- SeuratObject::LayerData(all_data,
                                   assay = "Lineage",
                                   layer = "counts")
lin_mat <- lin_mat[Matrix::rowSums(lin_mat) > 0,]
lin_mat <- lin_mat[,Matrix::colSums(lin_mat) > 0]

res_clusters <- multiomeFate:::barcode_clustering(lin_mat)
lin_mat <- as.matrix(lin_mat)
lin_mat2 <- lin_mat
lin_mat <- multiomeFate:::barcode_combine(lin_mat = lin_mat,
                                          lineage_clusters = res_clusters$lineage_clusters,
                                          verbose = 1)

lin_mat_safe <- lin_mat 
dim(lin_mat_safe)# dimension 7223 lineages x 65256 cells
lin_mat <- Matrix::Matrix(lin_mat, sparse = TRUE)

##############################################

all_data <- multiomeFate:::data_loader(which_files = c("fasttopics"))
posterior_mat <- data.frame(
  lineage = all_data$assigned_lineage,
  posterior = all_data$assigned_posterior
)
rownames(posterior_mat) <- Seurat::Cells(all_data)

# 44358 assigned cells
# 3286 lineages

##############################################

# strategy 1: assign to a unique maximum
assignment_vec_1 <- sapply(1:ncol(lin_mat), function(i){
  if(i %% floor(ncol(lin_mat)/10) == 0) cat('*')
  
  indices <-  values <- multiomeFate:::.nonzero_col(lin_mat, col_idx = i, bool_value = FALSE)
  values <- multiomeFate:::.nonzero_col(lin_mat, col_idx = i, bool_value = TRUE)
  stopifnot(length(indices) == length(values))
  
  if(length(which(values == max(values))) == 1){
    idx <- which.max(values)
    indices[idx]
  } else {
    return(NA)
  }
})
names(assignment_vec_1) <- colnames(lin_mat)
table(!is.na(assignment_vec_1)) # 52293 assigned cells 
length(unique(assignment_vec_1[!is.na(assignment_vec_1)])) # 2932 lineages

# strategy 2: assign to maximum larger than a count of 2 than second-largest
assignment_vec_2 <- sapply(1:ncol(lin_mat), function(i){
  if(i %% floor(ncol(lin_mat)/10) == 0) cat('*')
  
  indices <-  values <- multiomeFate:::.nonzero_col(lin_mat, col_idx = i, bool_value = FALSE)
  values <- multiomeFate:::.nonzero_col(lin_mat, col_idx = i, bool_value = TRUE)
  stopifnot(length(indices) == length(values))
  
  sorted_tmp_idx <- order(values, decreasing = TRUE)
  values <- values[sorted_tmp_idx]
  indices <- indices[sorted_tmp_idx]
  if(length(values) == 1 || values[1] >= values[2]+2){
    indices[1]
  } else {
    return(NA)
  }
})
names(assignment_vec_2) <- colnames(lin_mat)
table(!is.na(assignment_vec_2)) # 44455 assigned cells
length(unique(assignment_vec_2[!is.na(assignment_vec_2)])) # 2802 lineages

# strategy 3: assign to maximum larger than 2x than second-largest
assignment_vec_3 <- sapply(1:ncol(lin_mat), function(i){
  if(i %% floor(ncol(lin_mat)/10) == 0) cat('*')
  
  indices <-  values <- multiomeFate:::.nonzero_col(lin_mat, col_idx = i, bool_value = FALSE)
  values <- multiomeFate:::.nonzero_col(lin_mat, col_idx = i, bool_value = TRUE)
  stopifnot(length(indices) == length(values))
  
  sorted_tmp_idx <- order(values, decreasing = TRUE)
  values <- values[sorted_tmp_idx]
  indices <- indices[sorted_tmp_idx]
  if(length(values) == 1 || values[1] >= 2*values[2]){
    indices[1]
  } else {
    return(NA)
  }
})
names(assignment_vec_3) <- colnames(lin_mat)
table(!is.na(assignment_vec_3)) # 44190 assigned cells
length(unique(assignment_vec_3[!is.na(assignment_vec_3)])) # 2787 lineages

# strategy 4: assign to maximum larger than 3x than second-largest
assignment_vec_4 <- sapply(1:ncol(lin_mat), function(i){
  if(i %% floor(ncol(lin_mat)/10) == 0) cat('*')
  
  indices <-  values <- multiomeFate:::.nonzero_col(lin_mat, col_idx = i, bool_value = FALSE)
  values <- multiomeFate:::.nonzero_col(lin_mat, col_idx = i, bool_value = TRUE)
  stopifnot(length(indices) == length(values))
  
  sorted_tmp_idx <- order(values, decreasing = TRUE)
  values <- values[sorted_tmp_idx]
  indices <- indices[sorted_tmp_idx]
  if(length(values) == 1 || values[1] >= 3*values[2]){
    indices[1]
  } else {
    return(NA)
  }
})
names(assignment_vec_4) <- colnames(lin_mat)
table(!is.na(assignment_vec_4)) # 34731 assigned cells
length(unique(assignment_vec_4[!is.na(assignment_vec_4)])) # 2592 lineages


###########################

# what is the disagreement between assignment_vec_3 and posterior_mat?
assignment_vec_3_subset <- assignment_vec_3[rownames(posterior_mat)]
tab_mat <- table(assignment_vec_3_subset, posterior_mat$lineage)

tab_mat <- tab_mat[order(rowSums(tab_mat), decreasing = TRUE),]
tab_mat <- tab_mat[,order(colSums(tab_mat), decreasing = TRUE)]

col_ordering <- as.numeric(clue::solve_LSAP(tab_mat, maximum = TRUE))
tab_mat_new <- tab_mat[,col_ordering]
tab_mat_new <- tab_mat_new[1:min(dim(tab_mat_new)), 1:min(dim(tab_mat_new))]
sum(diag(tab_mat_new))/sum(tab_mat_new)

###########################

# make barplots
df <- data.frame(
  "Assign_plus1" = sum(!is.na(assignment_vec_1)),
  "Assign_plus2" = sum(!is.na(assignment_vec_2)),
  "Assign_times2" = sum(!is.na(assignment_vec_3)),
  "Assign_times3" = sum(!is.na(assignment_vec_4)),
  "Assign_posterior" = length(posterior_mat$lineage[!is.na(posterior_mat$lineage)])
)

library(tidyverse)      # ggplot2, dplyr, tidyr, readr, etc.
library(scales)         # nicer number formatting

# 2. Tidy ── convert to long format
df_long <- df %>%
  pivot_longer(everything(),
               names_to   = "Method",
               values_to  = "Count") %>%
  mutate(Method = str_replace(Method, "Assign_", "Assign "))  # nicer labels
df_long$Method <- factor(df_long$Method,
                         levels = c("Assign plus1",
                                    "Assign plus2",
                                    "Assign times2",
                                    "Assign times3",
                                    "Assign posterior"))

# 3. Plot ── barplot with minimal, Cell-style theme
p <- ggplot(df_long, aes(x = Method, y = Count, fill = Method)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.4) +
  # Add counts above bars
  geom_text(aes(label = comma(Count)),
            vjust = -0.4, size = 3.5, family = "Helvetica") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                     labels = comma_format()) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(y = "Number of Cells Assigned", x = NULL) +
  coord_cartesian(clip = "off") +
  theme_classic(base_family = "Helvetica", base_size = 10) +
  theme(
    axis.text.x  = element_text(angle = 25, hjust = 1),
    axis.ticks.x = element_blank(),
    axis.line.x  = element_blank(),
    plot.margin  = margin(10, 10, 15, 10),  # extra space for labels
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3)
  )

# 4. Save ── high-resolution export (300 dpi for print)
ggsave(paste0(plot_folder, "Writeup18_cloneclean_number-cells.pdf"),  
       p, width = 90,  height = 70,  units = "mm")
ggsave(paste0(plot_folder, "Writeup18_cloneclean_number-cells.png"),  
       p, width = 90,  height = 70,  units = "mm", dpi = 300)

###########################

# make barplots
df <- data.frame(
  "Assign_plus1" = length(unique(assignment_vec_1[!is.na(assignment_vec_1)])) ,
  "Assign_plus2" = length(unique(assignment_vec_2[!is.na(assignment_vec_2)])) ,
  "Assign_times2" = length(unique(assignment_vec_3[!is.na(assignment_vec_3)])) ,
  "Assign_times3" = length(unique(assignment_vec_4[!is.na(assignment_vec_4)])) ,
  "Assign_posterior" = length(unique(posterior_mat$lineage[!is.na(posterior_mat$lineage)]))
)

# 2. Tidy ── convert to long format
df_long <- df %>%
  pivot_longer(everything(),
               names_to   = "Method",
               values_to  = "Count") %>%
  mutate(Method = str_replace(Method, "Assign_", "Assign "))  # nicer labels
df_long$Method <- factor(df_long$Method,
                         levels = c("Assign plus1",
                                    "Assign plus2",
                                    "Assign times2",
                                    "Assign times3",
                                    "Assign posterior"))

# 3. Plot ── barplot with minimal, Cell-style theme
p <- ggplot(df_long, aes(x = Method, y = Count, fill = Method)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.4) +
  # Add counts above bars
  geom_text(aes(label = comma(Count)),
            vjust = -0.4, size = 3.5, family = "Helvetica") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                     labels = comma_format()) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(y = "Number of Lineages Detected", x = NULL) +
  coord_cartesian(clip = "off") +
  theme_classic(base_family = "Helvetica", base_size = 10) +
  theme(
    axis.text.x  = element_text(angle = 25, hjust = 1),
    axis.ticks.x = element_blank(),
    axis.line.x  = element_blank(),
    plot.margin  = margin(10, 10, 15, 10),  # extra space for labels
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3)
  )

# 4. Save ── high-resolution export (300 dpi for print)
ggsave(paste0(plot_folder, "Writeup18_cloneclean_number-lineages.pdf"),  
       p, width = 90,  height = 70,  units = "mm")
ggsave(paste0(plot_folder, "Writeup18_cloneclean_number-lineages.png"),  
       p, width = 90,  height = 70,  units = "mm", dpi = 300)

