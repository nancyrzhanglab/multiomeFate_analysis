rm(list=ls())
library(Seurat)
library(multiomeFate)

out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup15/"
plot_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/git/multiomeFate_analysis_kevin/fig/kevin/Writeup15/"
all_data <- multiomeFate:::data_loader(which_files = c("fasttopics"))

treatment_vec <- c("CIS", "COCL2", "DABTRAM")
day_early_vec <- c("d0", "d10")
modality_vec <- c("rna", "atac")

for(treatment in treatment_vec){
  for(day_early in day_early_vec){
    
    if(day_early == "d0"){
      day_later <- "d10"
      day_later_full <- paste0("day10_", treatment)
    } else {
      day_later <- "w5"
      day_later_full <- paste0("week5_", treatment)
    }
    
    for(modality in modality_vec){
      load(paste0(out_folder, "Writeup15_", treatment, "-from-", day_early, "_fatepotential-", modality, ".RData"))
      
      tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
      lineage_imputed_count <- final_fit$lineage_imputed_count
      lineage_future_count <- tab_mat[names(lineage_imputed_count), day_later_full]
      
      lineage_future_count <- lineage_future_count[names(lineage_imputed_count)]
      
      plot1 <- multiomeFate:::plot_lineageScatterplot(
        lineage_future_count = lineage_future_count,
        lineage_imputed_count = lineage_imputed_count,
        title = paste0(
          treatment, " ", day_later, "(", modality, ") growth potential of ", day_early, 
          " cells\n(Ridge for RNA fasttopics, ATAC PeakVI), (Log-scale)")
      )
      
      ggplot2::ggsave(filename = paste0(plot_folder, "Writeup15_lineage_", treatment, "_", day_early, "_", modality, ".png"),
                      plot1, device = "png", width = 10, height = 10, units = "in")
      
    }
  }
}
