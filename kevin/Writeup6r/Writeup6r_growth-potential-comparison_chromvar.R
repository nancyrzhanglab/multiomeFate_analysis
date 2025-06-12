rm(list=ls())

treatment_vec <- c("CIS", "COCL2", "DABTRAM")
day_early_vec <- c("day0", "day10")

motif_corr_list <- vector("list", length = 6)
names(motif_corr_list) <- c(paste0(day_early_vec[1], "_", treatment_vec),
                           paste0(day_early_vec[2], "_", treatment_vec))

cell_growth_potential_list <- vector("list", length = 6)
names(cell_growth_potential_list) <- names(motif_corr_list)

for(day_early in day_early_vec){
  if(day_early == "day0"){
    load("~/nzhanglab/project/emiliac/lineage_trace/task_identify_features_predictive_growth_and_lineage_specific/data/chromVar_day0_data.RData")
  } else {
    load("~/nzhanglab/project/emiliac/lineage_trace/task_identify_features_predictive_growth_and_lineage_specific/data/chromVar_day10_data.RData")
    chromvar_list <- list(
      CIS = chromvar_results_cis,
      COCL2 = chromvar_results_cocl2,
      DABTRAM = chromvar_results_dabtram
    )
  }
  
  for(treatment in treatment_vec){
    day_treatment <- paste0(day_early[1], "_", treatment)
    print(day_treatment)
    
    load(paste0("../../../../out/kevin/Writeup6r/Writeup6r_", treatment, "_", day_early, "_lineage-imputation_postprocess.RData"))
    
    if(day_early == "day10") {
      chromvar_mat <- t(chromvar_list[[treatment]])
    } else {
      chromvar_mat <- t(chromvar_results_day0)
    }
    
    # intersect cell names
    cell_names <- intersect(rownames(chromvar_mat), names(cell_imputed_score))
    cell_imputed_score <- cell_imputed_score[cell_names]
    chromvar_mat <- chromvar_mat[cell_names,]
    
    corr_vec <- sapply(1:ncol(chromvar_mat), function(j){
      stats::cor(cell_imputed_score, chromvar_mat[,j])
    })
    names(corr_vec) <- colnames(chromvar_mat)
    corr_vec[is.na(corr_vec)] <- 0
    
    motif_corr_list[[day_treatment]] <- corr_vec
    cell_growth_potential_list[[day_treatment]] <- cell_imputed_score
  }
}

##########################

round(stats::cor(cbind(motif_corr_list[["day0_CIS"]], 
                       motif_corr_list[["day0_COCL2"]], 
                       motif_corr_list[["day0_DABTRAM"]])), 2)

round(stats::cor(cbind(motif_corr_list[["day10_CIS"]], 
                       motif_corr_list[["day10_COCL2"]], 
                       motif_corr_list[["day10_DABTRAM"]])), 2)

round(stats::cor(cbind(motif_corr_list[["day0_CIS"]], 
                       motif_corr_list[["day10_CIS"]])), 2)
round(stats::cor(cbind(motif_corr_list[["day0_COCL2"]], 
                       motif_corr_list[["day10_COCL2"]])), 2)
round(stats::cor(cbind(motif_corr_list[["day0_DABTRAM"]], 
                       motif_corr_list[["day10_DABTRAM"]])), 2)

#########################

df <- data.frame(cis_cor = motif_corr_list[["day0_CIS"]],
                 cocl2_cor = motif_corr_list[["day0_COCL2"]],
                 dabtram_cor = motif_corr_list[["day0_DABTRAM"]])
p1 <- GGally::ggpairs(df, 
                      lower = list(continuous = GGally::wrap("points", alpha = 0.2, shape = 16)),
                      progress = FALSE) 
ggplot2::ggsave(filename = "../../../../out/figures/kevin/Writeup6r/Writeup6r_growth-potential_day0-pairs_motif-corr.png",
                p1, device = "png", width = 6, height = 6, units = "in")

df <- data.frame(cis_cor = motif_corr_list[["day10_CIS"]],
                 cocl2_cor = motif_corr_list[["day10_COCL2"]],
                 dabtram_cor = motif_corr_list[["day10_DABTRAM"]])
p1 <- GGally::ggpairs(df, 
                      lower = list(continuous = GGally::wrap("points", alpha = 0.2, shape = 16)),
                      progress = FALSE) 
ggplot2::ggsave(filename = "../../../../out/figures/kevin/Writeup6r/Writeup6r_growth-potential_day10-pairs_motif-corr.png",
                p1, device = "png", width = 6, height = 6, units = "in")

####

df <- data.frame(day0 = motif_corr_list[["day0_CIS"]],
                 day10 = motif_corr_list[["day10_CIS"]])
p1 <- GGally::ggpairs(df, 
                      lower = list(continuous = GGally::wrap("points", alpha = 0.2, shape = 16)),
                      progress = FALSE) 
ggplot2::ggsave(filename = "../../../../out/figures/kevin/Writeup6r/Writeup6r_growth-potential_cis-pairs_motif-corr.png",
                p1, device = "png", width = 6, height = 6, units = "in")


df <- data.frame(day0 = motif_corr_list[["day0_COCL2"]],
                 day10 = motif_corr_list[["day10_COCL2"]])
p1 <- GGally::ggpairs(df, 
                      lower = list(continuous = GGally::wrap("points", alpha = 0.2, shape = 16)),
                      progress = FALSE) 
ggplot2::ggsave(filename = "../../../../out/figures/kevin/Writeup6r/Writeup6r_growth-potential_cocl2-pairs_motif-corr.png",
                p1, device = "png", width = 6, height = 6, units = "in")


df <- data.frame(day0 = motif_corr_list[["day0_DABTRAM"]],
                 day10 = motif_corr_list[["day10_DABTRAM"]])
p1 <- GGally::ggpairs(df, 
                      lower = list(continuous = GGally::wrap("points", alpha = 0.2, shape = 16)),
                      progress = FALSE) 
ggplot2::ggsave(filename = "../../../../out/figures/kevin/Writeup6r/Writeup6r_growth-potential_dabtram-pairs_motif-corr.png",
                p1, device = "png", width = 6, height = 6, units = "in")

##############################

save(cell_growth_potential_list,
     motif_corr_list,
     date_of_run, session_info,
     file = paste0("../../../../out/kevin/Writeup6r/Writeup6r_all_cell-growth-potential_motif-correlation_output.RData"))


