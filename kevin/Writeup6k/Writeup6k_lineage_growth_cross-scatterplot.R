rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")

lineage_names <- sort(unique(all_data$assigned_lineage))

load("../../../../out/kevin/Writeup6j/Writeup6j_DABTRAM_day0_lineage-imputation_stepdown-LOOCV_concise-postprocessed.RData")
dabtram_day0 <- sapply(lineage_names, function(lineage_name){
  idx <- which(all_data$assigned_lineage == lineage_name)
  cell_names <- colnames(all_data)[idx]
  cell_names <- intersect(cell_names, names(cell_imputed_count))
    
  if(length(cell_names) == 0) return(NA)
  if(all(is.na(cell_imputed_count[cell_names]))) return(NA)
  mean(cell_imputed_count[cell_names], na.rm = T)
})
names(dabtram_day0) <- lineage_names

load("../../../../out/kevin/Writeup6j/Writeup6j_COCL2_day0_lineage-imputation_stepdown-LOOCV_concise-postprocessed.RData")
cocl2_day0 <- sapply(lineage_names, function(lineage_name){
  idx <- which(all_data$assigned_lineage == lineage_name)
  cell_names <- colnames(all_data)[idx]
  cell_names <- intersect(cell_names, names(cell_imputed_count))
  
  if(length(cell_names) == 0) return(NA)
  if(all(is.na(cell_imputed_count[cell_names]))) return(NA)
  mean(cell_imputed_count[cell_names], na.rm = T)
})
names(cocl2_day0) <- lineage_names

load("../../../../out/kevin/Writeup6j/Writeup6j_COCL2_day10_lineage-imputation_stepdown_concise-postprocessed.RData")
cocl2_day10 <- sapply(lineage_names, function(lineage_name){
  idx <- which(all_data$assigned_lineage == lineage_name)
  cell_names <- colnames(all_data)[idx]
  cell_names <- intersect(cell_names, names(cell_imputed_count))
  
  if(length(cell_names) == 0) return(NA)
  if(all(is.na(cell_imputed_count[cell_names]))) return(NA)
  mean(cell_imputed_count[cell_names], na.rm = T)
})
names(cocl2_day10) <- lineage_names

load("../../../../out/kevin/Writeup6k/Writeup6k_DABTRAM_day10_lineage-imputation_stepdown_concise-postprocessed.RData")
dabtram_day10 <- sapply(lineage_names, function(lineage_name){
  idx <- which(all_data$assigned_lineage == lineage_name)
  cell_names <- colnames(all_data)[idx]
  cell_names <- intersect(cell_names, names(cell_imputed_count))
  
  if(length(cell_names) == 0) return(NA)
  if(all(is.na(cell_imputed_count[cell_names]))) return(NA)
  mean(cell_imputed_count[cell_names], na.rm = T)
})
names(dabtram_day10) <- lineage_names

df <- cbind(cocl2_day0, dabtram_day0, cocl2_day10, dabtram_day10)
colnames(df) <- c("day0_COCL2", "day0_DABTRAM", "day10_COCL2", "day10_DABTRAM")
df <- as.data.frame(df)

bool_vec <- sapply(1:nrow(df), function(i){
  all(!is.na(df[i,]))
})
table(bool_vec)

#####

plot1 <- GGally::ggpairs(df)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6k/Writeup6k_lineage-count_growth-potential.png"),
                plot1, device = "png", width = 7, height = 7, units = "in")

################################
################################
################################


load("../../../../out/kevin/Writeup6j/Writeup6j_DABTRAM_day0_lineage-imputation_stepdown-LOOCV_concise-postprocessed.RData")
dabtram_day0 <- cell_imputed_count

load("../../../../out/kevin/Writeup6j/Writeup6j_COCL2_day0_lineage-imputation_stepdown-LOOCV_concise-postprocessed.RData")
cocl2_day0 <- cell_imputed_count

df <- cbind(cocl2_day0, dabtram_day0)
colnames(df) <- c("day0_COCL2", "day0_DABTRAM")
df <- as.data.frame(df)

plot1 <- GGally::ggpairs(df)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6k/Writeup6k_cell-count_growth-potential.png"),
                plot1, device = "png", width = 7, height = 7, units = "in")

df <- exp(df)

plot1 <- GGally::ggpairs(df)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6k/Writeup6k_cell-count_growth-potential_exponentiated.png"),
                plot1, device = "png", width = 7, height = 7, units = "in")

