rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

tp_early <- "day0"
treatment <- "CIS"
load(paste0("../../../../out/kevin/Writeup6n/Writeup6n_", treatment, "_", tp_early, "_lineage-imputation_stepdown_concise-postprocessed.RData"))
cis_d0_imputed <- cell_imputed_count[!is.na(cell_imputed_count)]
treatment <- "COCL2"
load(paste0("../../../../out/kevin/Writeup6n/Writeup6n_", treatment, "_", tp_early, "_lineage-imputation_stepdown_concise-postprocessed.RData"))
cocl2_d0_imputed <- cell_imputed_count[!is.na(cell_imputed_count)]
treatment <- "DABTRAM"
load(paste0("../../../../out/kevin/Writeup6n/Writeup6n_", treatment, "_", tp_early, "_lineage-imputation_stepdown_concise-postprocessed.RData"))
dabtram_d0_imputed <- cell_imputed_count[!is.na(cell_imputed_count)]

tp_early <- "day10"
treatment <- "CIS"
load(paste0("../../../../out/kevin/Writeup6n/Writeup6n_", treatment, "_", tp_early, "_lineage-imputation_stepdown_concise-postprocessed.RData"))
cis_d10_imputed <- cell_imputed_count[!is.na(cell_imputed_count)]
treatment <- "COCL2"
load(paste0("../../../../out/kevin/Writeup6n/Writeup6n_", treatment, "_", tp_early, "_lineage-imputation_stepdown_concise-postprocessed.RData"))
cocl2_d10_imputed <- cell_imputed_count[!is.na(cell_imputed_count)]
treatment <- "DABTRAM"
load(paste0("../../../../out/kevin/Writeup6n/Writeup6n_", treatment, "_", tp_early, "_lineage-imputation_stepdown_concise-postprocessed.RData"))
dabtram_d10_imputed <- cell_imputed_count[!is.na(cell_imputed_count)]

load("../../../../out/kevin/Writeup6l/Writeup6l_chromVar_rna-chromvar_lightweight.RData")

imputed_list <- list(
  day0_CIS = cis_d0_imputed,
  day0_COCL2 = cocl2_d0_imputed,
  day0_DABTRAM = dabtram_d0_imputed,
  day10_CIS = cis_d10_imputed,
  day10_COCL2 = cocl2_d10_imputed,
  day10_DABTRAM = dabtram_d10_imputed
)

########

correlation_list <- lapply(imputed_list, function(vec){
  rna_mat <- t(all_data[["Saver"]]@data[,names(vec)])
  cor_vec <- sapply(1:ncol(rna_mat), function(j){
    stats::cor(vec, rna_mat[,j])
  })
  names(cor_vec) <- colnames(rna_mat)
  cor_vec[is.na(cor_vec)] <- 0
  
  cor_vec
})

gene_names <- names(correlation_list[[1]])
for(i in 1:length(correlation_list)){
  gene_names <- intersect(gene_names, names(correlation_list[[i]]))
}

for(i in 1:length(correlation_list)){
  correlation_list[[i]] <- correlation_list[[i]][gene_names]
}

tmp <- do.call(cbind, correlation_list)
round(100*cor(tmp))

treatment_vec <- c("CIS", "COCL2", "DABTRAM")
for(treatment in treatment_vec){
  day0 <- correlation_list[[paste0("day0_", treatment)]]
  day10 <- correlation_list[[paste0("day10_", treatment)]]
  
  df <- data.frame(day0 = day0, day10 = day10)
  p1 <- ggplot2::ggplot(df, ggplot2::aes(x = day0, 
                                         y = day10))
  p1 <- p1 + ggplot2::geom_point(color = "palevioletred")
  p1 <- p1 + ggplot2::ggtitle(paste0(treatment, ": Correlation ", round(stats::cor(day0, day10), 2)))
  p1 <- p1 + ggplot2::labs(y = "D10 gene corr w/ W5 growth potent.", x = "D0 gene corr w/ D10 growth potent.")
  
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6n/Writeup6n_lineage-imputation_", treatment, "_d0-d10_gene-corr.png"),
                  p1, device = "png", width = 6, height = 6, units = "in")
}

#################

source("../Writeup6b/gene_list.R")
keygenes <- unique(unlist(keygenes))

for(treatment in treatment_vec){
  day0 <- correlation_list[[paste0("day0_", treatment)]]
  day10 <- correlation_list[[paste0("day10_", treatment)]]
  
  labeling <- rep(0, length(day0))
  labeling[names(day0) %in% keygenes] <- 1
  
  df <- data.frame(day0 = day0,
                   day10 = day10,
                   name = names(day0),
                   labeling = labeling)
  df[,"labeling"] <- as.factor(df[,"labeling"])
  # put all the labeling == TRUE on bottom
  df <- df[c(which(df[,"labeling"] == 0),  which(df[,"labeling"] == 1)),]
  
  p1 <- ggplot2::ggplot(df, ggplot2::aes(x = day0,
                                         y = day10))
  p1 <- p1 + ggplot2::geom_point(ggplot2::aes(color = labeling))
  p1 <- p1 + ggplot2::scale_colour_manual(values=c("black", "red"))
  p1 <- p1 + ggrepel::geom_text_repel(data = subset(df, labeling == "1"),
                                      ggplot2::aes(label = name, color = labeling),
                                      box.padding = ggplot2::unit(0.5, 'lines'),
                                      point.padding = ggplot2::unit(1.6, 'lines'),
                                      max.overlaps = 50)
  p1 <- p1 + ggplot2::ggtitle(paste0(treatment, ": Correlation ", round(stats::cor(day0, day10), 2)))
  p1 <- p1 + ggplot2::labs(y = "D10 gene corr w/ W5 growth potent.", x = "D0 gene corr w/ D10 growth potent.")
  p1 <- p1 + Seurat::NoLegend()
  
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6n/Writeup6n_lineage-imputation_", treatment, "_d0-d10_gene-corr_ggrepel.png"),
                  p1, device = "png", width = 6, height = 6, units = "in")
}
