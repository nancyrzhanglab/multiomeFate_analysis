rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6k/Writeup6k_DABTRAM_day10_lineage-imputation_stepdown_concise-postprocessed.RData")
dabtram_imputed <- cell_imputed_count[!is.na(cell_imputed_count)]
dabtram_fit <- fit
load("../../../../out/kevin/Writeup6j/Writeup6j_COCL2_day10_lineage-imputation_stepdown_concise-postprocessed.RData")
cocl2_imputed <- cell_imputed_count[!is.na(cell_imputed_count)]
cocl2_fit <- fit

length(intersect(names(dabtram_imputed)[!is.na(dabtram_imputed)],
                 names(cocl2_imputed)[!is.na(cocl2_imputed)]))

load("../../../../out/kevin/Writeup6n/Writeup6n_topics.RData")

cocl2_loading <- cocl2_fasttopics@feature.loadings
dabtram_loading <- dabtram_fasttopics@feature.loadings

dabtram_coef <- dabtram_fit$coefficient_vec
dabtram_coef <- abs(dabtram_coef[grep("fastTopic", names(dabtram_coef))])
cocl2_coef <- cocl2_fit$coefficient_vec
cocl2_coef <- abs(cocl2_coef[grep("fastTopic", names(cocl2_coef))])

dabtram_weights <- (dabtram_loading[,names(dabtram_coef)] %*% dabtram_coef)[,1]
cocl2_weights <- (cocl2_loading[,names(cocl2_coef)] %*% cocl2_coef)[,1]
common_genes <- intersect(names(dabtram_weights), names(cocl2_weights))
dabtram_weights <- dabtram_weights[common_genes]
cocl2_weights <- cocl2_weights[common_genes]

stats::cor(dabtram_weights, cocl2_weights)

df <- data.frame(cocl2_geneweight = cocl2_weights,
                 dabtram_geneweight = dabtram_weights)
p1 <- GGally::ggpairs(df, 
                      lower = list(continuous = GGally::wrap("points", alpha = 0.2, shape = 16)),
                      progress = FALSE) + ggplot2::theme_bw()
ggplot2::ggsave(filename = "../../../../out/figures/Writeup6n/Writeup6n_lineage-imputation_day10-comparison_geneweights.png",
                p1, device = "png", width = 6, height = 6, units = "in")

###########

labeling <- rep(0, length(cocl2_weights))
labeling[order(cocl2_weights, decreasing = T)[1:20]] <- 1
labeling[order(dabtram_weights, decreasing = T)[1:20]] <- 1

df <- data.frame(cocl2_geneweight = cocl2_weights,
                 dabtram_geneweight = dabtram_weights,
                 name = names(cocl2_weights),
                 labeling = labeling)
df[,"labeling"] <- as.factor(df[,"labeling"])
# put all the labeling == TRUE on bottom
df <- df[c(which(df[,"labeling"] == 0),  which(df[,"labeling"] == 1)),]

p1 <- ggplot2::ggplot(df, ggplot2::aes(x = cocl2_geneweight, 
                                       y = dabtram_geneweight))
p1 <- p1 + ggplot2::geom_point(ggplot2::aes(color = labeling))
p1 <- p1 + ggplot2::scale_colour_manual(values=c("black", "red"))
p1 <- p1 + ggrepel::geom_text_repel(data = subset(df, labeling == "1"),
                                    ggplot2::aes(label = name, color = labeling),
                                    box.padding = ggplot2::unit(0.5, 'lines'),
                                    point.padding = ggplot2::unit(1.6, 'lines'),
                                    max.overlaps = 50)
p1 <- p1 + ggplot2::ggtitle(paste0("DABTRAM vs COCL2 (Day10, Gene weights, RNA)")) 
p1 <- p1 + Seurat::NoLegend()
p1 <- p1 + ggplot2::labs(y = "DABTRAM", x = "COCL2")

ggplot2::ggsave(filename = "../../../../out/figures/Writeup6n/Writeup6n_lineage-imputation_day10-comparison_geneweights_ggrepel_COCL2-DABTRAM.png",
                p1, device = "png", width = 6, height = 6, units = "in")

###########################

load("../../../../out/kevin/Writeup6l/Writeup6l_chromVar_rna-chromvar_lightweight.RData")

rna_mat <- t(all_data[["Saver"]]@data[,names(dabtram_imputed)])
dabtram_cor_vec <- sapply(1:ncol(rna_mat), function(j){
  stats::cor(dabtram_imputed, rna_mat[,j])
})
names(dabtram_cor_vec) <- colnames(rna_mat)
dabtram_cor_vec[is.na(dabtram_cor_vec)] <- 0

rna_mat <- t(all_data[["Saver"]]@data[,names(cocl2_imputed)])
cocl2_cor_vec <- sapply(1:ncol(rna_mat), function(j){
  stats::cor(cocl2_imputed, rna_mat[,j])
})
names(cocl2_cor_vec) <- colnames(rna_mat)
cocl2_cor_vec[is.na(cocl2_cor_vec)] <- 0

stats::cor(cocl2_cor_vec, dabtram_cor_vec)

labeling <- rep(0, length(cocl2_cor_vec))
labeling[order(cocl2_cor_vec, decreasing = T)[1:15]] <- 1
labeling[order(dabtram_cor_vec, decreasing = T)[1:15]] <- 1
labeling[order(-cocl2_cor_vec, decreasing = T)[1:15]] <- 1
labeling[order(-dabtram_cor_vec, decreasing = T)[1:15]] <- 1

df <- data.frame(cocl2_cor = cocl2_cor_vec,
                 dabtram_cor = dabtram_cor_vec,
                 name = names(cocl2_cor_vec),
                 labeling = labeling)
df[,"labeling"] <- as.factor(df[,"labeling"])
# put all the labeling == TRUE on bottom
df <- df[c(which(df[,"labeling"] == 0),  which(df[,"labeling"] == 1)),]

p1 <- ggplot2::ggplot(df, ggplot2::aes(x = cocl2_cor, 
                                       y = dabtram_cor))
p1 <- p1 + ggplot2::geom_point(ggplot2::aes(color = labeling))
p1 <- p1 + ggplot2::scale_colour_manual(values=c("black", "red"))
p1 <- p1 + ggrepel::geom_text_repel(data = subset(df, labeling == "1"),
                                    ggplot2::aes(label = name, color = labeling),
                                    box.padding = ggplot2::unit(0.5, 'lines'),
                                    point.padding = ggplot2::unit(1.6, 'lines'),
                                    max.overlaps = 50)
p1 <- p1 + ggplot2::ggtitle(paste0("DABTRAM vs COCL2\n(RNA corr w/ Week5 growth potential computed at Day10 cells)\nCorrelation: ", 
                                   round(stats::cor(cocl2_cor_vec, dabtram_cor_vec), 2))) 
p1 <- p1 + Seurat::NoLegend()
p1 <- p1 + ggplot2::labs(y = "DABTRAM", x = "COCL2")

ggplot2::ggsave(filename = "../../../../out/figures/Writeup6n/Writeup6n_lineage-imputation_day10-comparison_gene-corr_ggrepel_COCL2-DABTRAM.png",
                p1, device = "png", width = 6, height = 6, units = "in")

