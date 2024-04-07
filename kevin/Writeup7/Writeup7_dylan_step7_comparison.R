rm(list=ls())
library(Seurat)

load("~/project/Multiome_fate/out/kevin/Writeup7/Writeup7_dylan_step3_fasttopics.RData")

treatment_vec <- sort(unique(all_data$OG_condition))
treatment_vec <- treatment_vec[treatment_vec != "naive"]
day_early <- "naive"

gene_corr_list <- vector("list", length = length(treatment_vec))
names(gene_corr_list) <- paste0(day_early, "_", treatment_vec)

cell_growth_potential_list <- vector("list", length = length(treatment_vec))
names(cell_growth_potential_list) <- names(gene_corr_list)

for(treatment in treatment_vec){
  day_treatment <- paste0(day_early, "_", treatment)
  print(day_treatment)
  
  load(paste0("~/project/Multiome_fate/out/kevin/Writeup7/Writeup7_dylan_step6_growth-potential_", treatment, "_postprocess.RData"))
  
  rna_mat <- t(all_data[["Saver"]]@scale.data[,names(cell_imputed_score)])
  corr_vec <- sapply(1:ncol(rna_mat), function(j){
    stats::cor(cell_imputed_score, rna_mat[,j])
  })
  names(corr_vec) <- colnames(rna_mat)
  corr_vec[is.na(corr_vec)] <- 0
  
  gene_corr_list[[day_treatment]] <- corr_vec
  cell_growth_potential_list[[day_treatment]] <- cell_imputed_score
}

#########################

df <- data.frame(do.call(cbind, gene_corr_list))
p1 <- GGally::ggpairs(df, 
                      mapping = ggplot2::aes(color = "coral1"),
                      lower = list(continuous = GGally::wrap("points", alpha = 0.2, shape = 16)),
                      progress = FALSE) 
ggplot2::ggsave(filename = "~/project/Multiome_fate/out/figures/Writeup7/Writeup7_growth-potential_pairs_gene-corr.png",
                p1, device = "png", width = 8, height = 8, units = "in")

write.csv(df,
          file = "~/project/Multiome_fate/out/kevin/Writeup7/Writeup7_growth-potential_pairs_gene-corr.csv" )

####

df <- data.frame(do.call(cbind, gene_corr_list))
col_vec <- rep("coral1", nrow(df)); names(col_vec) <- rownames(df)
cex_vec <- rep(1, nrow(df)); names(cex_vec) <- rownames(df)

gene_idx <- which(rownames(df) == "CD44")
col_vec[gene_idx] <- "black"
cex_vec[gene_idx] <- 3
gene_ordering <- c(setdiff(1:nrow(df), gene_idx), gene_idx)

df <- df[gene_ordering,]
col_vec <- col_vec[gene_ordering]
cex_vec <- cex_vec[gene_ordering]

png("~/project/Multiome_fate/out/figures/Writeup7/Writeup7_growth-potential_pairs_gene-corr_highlight_CD44.png",
    height = 2500, width = 2500, units = "px", res = 300)
graphics::pairs(df, col = col_vec, pch = 16, cex = cex_vec)
graphics.off()

####

file_vec <- c("acid_markers.csv",
              "cis_markers.csv",
              "cocl2_markers.csv",
              "dab_markers.csv",
              "dox_markers.csv",
              "markers_in_three_or_more.csv",
              "tram_markers.csv")
names(file_vec) <- c("acid", "cis", "cocl2", "dab", "dox", "three", "tram")

for(i in 1:length(file_vec)){
  print(names(file_vec)[i])
  dylan_vec <- read.csv(paste0("~/project/Multiome_fate/data/", file_vec[i]))[,2]
  
  df <- data.frame(do.call(cbind, gene_corr_list))
  col_vec <- rep("coral1", nrow(df)); names(col_vec) <- rownames(df)
  cex_vec <- rep(1, nrow(df)); names(cex_vec) <- rownames(df)
  
  gene_idx <- which(rownames(df) %in% dylan_vec)
  col_vec[gene_idx] <- "black"
  cex_vec[gene_idx] <- 2
  gene_ordering <- c(setdiff(1:nrow(df), gene_idx), gene_idx)
  
  df <- df[gene_ordering,]
  col_vec <- col_vec[gene_ordering]
  cex_vec <- cex_vec[gene_ordering]
  
  png(paste0("~/project/Multiome_fate/out/figures/Writeup7/Writeup7_growth-potential_pairs_gene-corr_highlight_",
             names(file_vec)[i],
             ".png"),
      height = 2500, width = 2500, units = "px", res = 300)
  graphics::pairs(df, col = col_vec, pch = 16, cex = cex_vec,
                  main = names(file_vec)[i])
  graphics.off()
}

######################################

cell_imputed_mat <- numeric(0)

for(treatment in treatment_vec){
  load(paste0("~/project/Multiome_fate/out/kevin/Writeup7/Writeup7_dylan_step6_growth-potential_", treatment, "_postprocess.RData"))
  cell_imputed_mat <- cbind(cell_imputed_mat, cell_imputed_score)
}

df <- data.frame(cell_imputed_mat)
colnames(df) <- treatment_vec

p1 <- GGally::ggpairs(df, 
                      lower = list(continuous = GGally::wrap("points", alpha = 0.2, shape = 16)),
                      progress = FALSE) 
ggplot2::ggsave(filename = "~/project/Multiome_fate/out/figures/Writeup7/Writeup7_growth-potential_pairs_cell-corr_logscale.png",
                p1, device = "png", width = 8, height = 8, units = "in")


df <- data.frame(pmin(10^cell_imputed_mat, 100))
colnames(df) <- treatment_vec

p1 <- GGally::ggpairs(df, 
                      lower = list(continuous = GGally::wrap("points", alpha = 0.2, shape = 16)),
                      progress = FALSE) 
ggplot2::ggsave(filename = "~/project/Multiome_fate/out/figures/Writeup7/Writeup7_growth-potential_pairs_cell-corr.png",
                p1, device = "png", width = 8, height = 8, units = "in")

