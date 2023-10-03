rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

tp_early <- "day10"
treatment <- "CIS"
load(paste0("../../../../out/kevin/Writeup6n/Writeup6n_", treatment, "_", tp_early, "_lineage-imputation_stepdown_concise-postprocessed.RData"))
cis_imputed <- cell_imputed_count[!is.na(cell_imputed_count)]
treatment <- "COCL2"
load(paste0("../../../../out/kevin/Writeup6n/Writeup6n_", treatment, "_", tp_early, "_lineage-imputation_stepdown_concise-postprocessed.RData"))
cocl2_imputed <- cell_imputed_count[!is.na(cell_imputed_count)]
treatment <- "DABTRAM"
load(paste0("../../../../out/kevin/Writeup6n/Writeup6n_", treatment, "_", tp_early, "_lineage-imputation_stepdown_concise-postprocessed.RData"))
dabtram_imputed <- cell_imputed_count[!is.na(cell_imputed_count)]

imputed_list <- list(CIS = cis_imputed,
                     COCL2 = cocl2_imputed,
                     DABTRAM = dabtram_imputed)

load("../../../../out/kevin/Writeup6l/Writeup6l_chromVar_rna-chromvar_lightweight.RData")

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
assigned_lineage <- all_data$assigned_lineage
dataset <- all_data$dataset

treatment_vec <- c("CIS", "COCL2", "DABTRAM")
proportion_list <- sapply(treatment_vec, function(treatment){
  desired_lineages <- rownames(tab_mat)[intersect(which(tab_mat[,paste0("week5_", treatment)] >= 10),
                                                  which(tab_mat[,paste0("day10_", treatment)] >= 10))]
  
  # find all the cells in that lineage
  proportion_vec <- sapply(desired_lineages, function(desired_lineage){
    cell_names <- names(assigned_lineage)[intersect(which(assigned_lineage == desired_lineage),
                                                    which(dataset == paste0("day10_", treatment)))]
    cell_expansion <- 10^(imputed_list[[treatment]][cell_names])
    cell_expansion <- sort(cell_expansion, decreasing = T)
    
    sum_val <- sum(cell_expansion)
    cum_sum <- cumsum(cell_expansion)
    
    idx <- min(which(cum_sum >= 0.75*sum_val))
    idx/length(cell_expansion)
  })
  
  sort(round(proportion_vec,2))
})
proportion_list

png("../../../../out/figures/Writeup6n/Writeup6n_lineage-imputation_day10-expansion.png",
    height = 1000, width = 1500, units = "px", res = 300)
hist(unlist(proportion_list), 
     main = "Considering only Day10 >= 10, Week5 >= 30",
     xlab = "Percentage of Day10 cells that yield >= 75% of Week5",
     ylab = "Frequency",
     cex.lab = 0.8)
mean_val <- mean(unlist(proportion_list))
lines(x = rep(mean_val, 2), y = c(-100,100), col = 2, lwd = 2, lty = 2)
graphics.off()
