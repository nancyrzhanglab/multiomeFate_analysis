rm(list=ls())
library(Seurat)

load("../../../../out/kevin/Writeup6l/Writeup6l_chromVar_rna-chromvar_lightweight.RData")

treatment_vec <- c("CIS", "COCL2", "DABTRAM")
day_early_vec <- c("day0", "day10")

gene_corr_list <- vector("list", length = 6)
names(gene_corr_list) <- c(paste0(day_early_vec[1], "_", treatment_vec),
                           paste0(day_early_vec[2], "_", treatment_vec))

cell_growth_potential_list <- vector("list", length = 6)
names(cell_growth_potential_list) <- names(gene_corr_list)

for(day_early in day_early_vec){
  for(treatment in treatment_vec){
    day_treatment <- paste0(day_early[1], "_", treatment)
    print(day_treatment)
    
    load(paste0("../../../../out/kevin/Writeup6r/Writeup6r_", treatment, "_", day_early, "_lineage-imputation_postprocess.RData"))
    
    rna_mat <- t(all_data[["Saver"]]@scale.data[,names(cell_imputed_score)])
    corr_vec <- sapply(1:ncol(rna_mat), function(j){
      stats::cor(cell_imputed_score, rna_mat[,j])
    })
    names(corr_vec) <- colnames(rna_mat)
    corr_vec[is.na(corr_vec)] <- 0
    
    gene_corr_list[[day_treatment]] <- corr_vec
    cell_growth_potential_list[[day_treatment]] <- cell_imputed_score
  }
}

##########################

round(stats::cor(cbind(gene_corr_list[["day0_CIS"]], 
                       gene_corr_list[["day0_COCL2"]], 
                       gene_corr_list[["day0_DABTRAM"]])), 2)

round(stats::cor(cbind(gene_corr_list[["day10_CIS"]], 
                       gene_corr_list[["day10_COCL2"]], 
                       gene_corr_list[["day10_DABTRAM"]])), 2)

round(stats::cor(cbind(gene_corr_list[["day0_CIS"]], 
                       gene_corr_list[["day10_CIS"]])), 2)
round(stats::cor(cbind(gene_corr_list[["day0_COCL2"]], 
                       gene_corr_list[["day10_COCL2"]])), 2)
round(stats::cor(cbind(gene_corr_list[["day0_DABTRAM"]], 
                       gene_corr_list[["day10_DABTRAM"]])), 2)

#########################

df <- data.frame(cis_cor = gene_corr_list[["day0_CIS"]],
                 cocl2_cor = gene_corr_list[["day0_COCL2"]],
                 dabtram_cor = gene_corr_list[["day0_DABTRAM"]])
p1 <- GGally::ggpairs(df, 
                      mapping = ggplot2::aes(color = "coral1"),
                      lower = list(continuous = GGally::wrap("points", alpha = 0.2, shape = 16)),
                      progress = FALSE) 
ggplot2::ggsave(filename = "../../../../out/figures/kevin/Writeup6r/Writeup6r_growth-potential_day0-pairs_gene-corr.png",
                p1, device = "png", width = 6, height = 6, units = "in")

df <- data.frame(cis_cor = gene_corr_list[["day10_CIS"]],
                 cocl2_cor = gene_corr_list[["day10_COCL2"]],
                 dabtram_cor = gene_corr_list[["day10_DABTRAM"]])
p1 <- GGally::ggpairs(df, 
                      mapping = ggplot2::aes(color = "coral1"),
                      lower = list(continuous = GGally::wrap("points", alpha = 0.2, shape = 16)),
                      progress = FALSE) 
ggplot2::ggsave(filename = "../../../../out/figures/kevin/Writeup6r/Writeup6r_growth-potential_day10-pairs_gene-corr.png",
                p1, device = "png", width = 6, height = 6, units = "in")

####

df <- data.frame(day0 = gene_corr_list[["day0_CIS"]],
                 day10 = gene_corr_list[["day10_CIS"]])
p1 <- GGally::ggpairs(df, 
                      mapping = ggplot2::aes(color = "coral1"),
                      lower = list(continuous = GGally::wrap("points", alpha = 0.2, shape = 16)),
                      progress = FALSE) 
ggplot2::ggsave(filename = "../../../../out/figures/kevin/Writeup6r/Writeup6r_growth-potential_cis-pairs_gene-corr.png",
                p1, device = "png", width = 6, height = 6, units = "in")


df <- data.frame(day0 = gene_corr_list[["day0_COCL2"]],
                 day10 = gene_corr_list[["day10_COCL2"]])
p1 <- GGally::ggpairs(df, 
                      mapping = ggplot2::aes(color = "coral1"),
                      lower = list(continuous = GGally::wrap("points", alpha = 0.2, shape = 16)),
                      progress = FALSE) 
ggplot2::ggsave(filename = "../../../../out/figures/kevin/Writeup6r/Writeup6r_growth-potential_cocl2-pairs_gene-corr.png",
                p1, device = "png", width = 6, height = 6, units = "in")


df <- data.frame(day0 = gene_corr_list[["day0_DABTRAM"]],
                 day10 = gene_corr_list[["day10_DABTRAM"]])
p1 <- GGally::ggpairs(df, 
                      mapping = ggplot2::aes(color = "coral1"),
                      lower = list(continuous = GGally::wrap("points", alpha = 0.2, shape = 16)),
                      progress = FALSE) 
ggplot2::ggsave(filename = "../../../../out/figures/kevin/Writeup6r/Writeup6r_growth-potential_dabtram-pairs_gene-corr.png",
                p1, device = "png", width = 6, height = 6, units = "in")

##############################

save(cell_growth_potential_list,
     gene_corr_list,
     date_of_run, session_info,
     file = paste0("../../../../out/kevin/Writeup6r/Writeup6r_all_cell-growth-potential_gene-correlation_output.RData"))


