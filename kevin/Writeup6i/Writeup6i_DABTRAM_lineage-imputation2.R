rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6h/Writeup6h_DABTRAM-day10_pseudotime.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

cell_features <- cbind(1, scale(all_data$pseudotime[colnames(all_data2)]))
colnames(cell_features) <- c("Intercept", "Pseudotime")
cell_lineage <- all_data2$assigned_lineage
cell_lineage <- cell_lineage[!is.na(cell_lineage)]
cell_features <- cell_features[names(cell_lineage),]
uniq_lineage <- sort(unique(cell_lineage))
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
lineage_current_count <- tab_mat[uniq_lineage,"day10_DABTRAM"]
lineage_future_count <- tab_mat[uniq_lineage,"week5_DABTRAM"]

# now try using the somewhat correlated features of the topic model
topic_mat <- all_data[["fasttopic_DABTRAM"]]@cell.embeddings[names(cell_lineage),]
pseudotime_vec <- all_data$pseudotime[names(cell_lineage)]

# let's try using all the topics #Clueless
cell_features <- cbind(1, scale(topic_mat), scale(pseudotime_vec))
p <- ncol(cell_features)
colnames(cell_features)[1] <- "Intercept"
colnames(cell_features)[p] <- "Pseudotime"

tmp <- quantile(abs(cell_features[,-1]), probs = 0.95)
coef_val <- 2*log(max(lineage_future_count/lineage_current_count))/((p-1)*tmp)
coefficient_initial <- c(0, rep(coef_val, p-1))
names(coefficient_initial) <- colnames(cell_features)

set.seed(10)
lineage_res <- multiomeFate:::lineage_imputation(cell_features = cell_features,
                                                 cell_lineage = cell_lineage,
                                                 coefficient_initial = coefficient_initial,
                                                 lineage_future_count = lineage_future_count,
                                                 random_initializations = 10,
                                                 verbose = 2)

cell_imputed_count <- as.numeric(exp(cell_features %*% lineage_res$fit$coefficient_vec))
names(cell_imputed_count) <- rownames(cell_features)
quantile(round(cell_imputed_count), probs = seq(0.9,1,length.out=11))
round(lineage_res$fit$coefficient_vec, 3)

lineage_imputed_count <- sapply(uniq_lineage, function(lineage){
  sum(cell_imputed_count[which(cell_lineage == lineage)])
})
stats::cor(lineage_imputed_count, lineage_future_count)
stats::cor(log1p(lineage_imputed_count), log1p(lineage_future_count))

imputed_vec <- rep(NA, ncol(all_data))
names(imputed_vec) <- colnames(all_data)
imputed_vec[names(cell_imputed_count)] <- log(cell_imputed_count)
all_data$imputed_count <- imputed_vec

# grDevices::colorRampPalette("lightgray", "blue")
p1 <- scCustomize::FeaturePlot_scCustom(all_data, 
                                        colors_use = list("lightgray", "blue"),
                                        na_cutoff = stats::quantile(all_data$imputed_count, probs = 0.05, na.rm = T),
                                        na_color = "bisque",
                                        reduction = "umap", 
                                        features = "imputed_count")
p1 <- p1 + ggplot2::ggtitle("DABTRAM Day10 imputed counts (Log-scale)\nAll RNA topics")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6i/Writeup6i_DABTRAM-day10_imputed-lineage_umap_all-fasttopics.png"),
                p1, device = "png", width = 7, height = 5, units = "in")

lineage_imputed_count2 <- log10(lineage_imputed_count+1)
lineage_future_count2 <- log10(lineage_future_count+1)

labeling_vec <- rep(FALSE, length(lineage_imputed_count2))
labeling_vec[which(lineage_imputed_count2 >= 1.5)] <- TRUE
labeling_vec[which(lineage_future_count2 >= 1.5)] <- TRUE
table(labeling_vec)

df <- data.frame(lineage_imputed_count = lineage_imputed_count2,
                 lineage_future_count = lineage_future_count2,
                 name = names(lineage_imputed_count2),
                 labeling = labeling_vec)
# put all the labeling == TRUE on bottom
df <- df[c(which(!df[,"labeling"]), which(df[,"labeling"])),]

p1 <- ggplot2::ggplot(df, ggplot2::aes(x = lineage_future_count, y = lineage_imputed_count))
p1 <- p1 + ggplot2::geom_point(ggplot2::aes(color = labeling))
p1 <- p1 + ggplot2::scale_colour_manual(values=c("black", "red"))
p1 <- p1 + ggrepel::geom_text_repel(data = subset(df, labeling == TRUE),
                                    ggplot2::aes(label = name, color = labeling),
                                    box.padding = ggplot2::unit(0.5, 'lines'),
                                    point.padding = ggplot2::unit(1.6, 'lines'),
                                    max.overlaps = 50)
p1 <- p1 + ggplot2::ggtitle(paste0("DABTRAM Day10 imputed counts (Log-scale)\nAll RNA topics, Corr: ", 
                                   round(stats::cor(lineage_imputed_count2, lineage_future_count2), 3))) +
  ggplot2::xlab("Observed lineage count (Log10)") + ggplot2::ylab("Predicted lineage count (Log10)")
p1 <- p1 + Seurat::NoLegend()

ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6i/Writeup6i_DABTRAM-day10_lineage-level_counts_all-fasttopics_labeled.png"),
                p1, device = "png", width = 10, height = 10, units = "in")


##################################

# Select only the ATAC vectors that are correlated with the tiers

# now try using the somewhat correlated features of the topic model
atac_mat <- all_data2[["lsi"]]@cell.embeddings

log10pval_vec <- sapply(1:ncol(atac_mat), function(j){
  x <- atac_mat[names(cell_lineage),j]
  y <- as.factor(all_data2$tier_vec[names(cell_lineage)])
  res <- stats::oneway.test(x ~ y)
  -log10(res$p.value)
})
names(log10pval_vec) <- colnames(atac_mat)
round(quantile(log10pval_vec))

# let's try using all the topics #Clueless
cell_features <- cbind(1, scale(topic_mat), 
                       scale(atac_mat[names(cell_lineage),which(log10pval_vec>=18)]), 
                       scale(pseudotime_vec))
p <- ncol(cell_features)
colnames(cell_features)[1] <- "Intercept"
colnames(cell_features)[p] <- "Pseudotime"

tmp <- quantile(abs(cell_features[,-1]), probs = 0.95)
coef_val <- 2*log(max(lineage_future_count/lineage_current_count))/((p-1)*tmp)
coefficient_initial <- c(0, rep(coef_val, p-1))
names(coefficient_initial) <- colnames(cell_features)

set.seed(10)
lineage_res <- multiomeFate:::lineage_imputation(cell_features = cell_features,
                                                 cell_lineage = cell_lineage,
                                                 coefficient_initial = coefficient_initial,
                                                 lineage_future_count = lineage_future_count,
                                                 random_initializations = 5,
                                                 verbose = 2)

cell_imputed_count <- as.numeric(exp(cell_features %*% lineage_res$fit$coefficient_vec))
names(cell_imputed_count) <- rownames(cell_features)
quantile(round(cell_imputed_count), probs = seq(0.9,1,length.out=11))
round(lineage_res$fit$coefficient_vec, 3)

lineage_imputed_count <- sapply(uniq_lineage, function(lineage){
  sum(cell_imputed_count[which(cell_lineage == lineage)])
})
stats::cor(lineage_imputed_count, lineage_future_count)
stats::cor(log1p(lineage_imputed_count), log1p(lineage_future_count))

imputed_vec <- rep(NA, ncol(all_data))
names(imputed_vec) <- colnames(all_data)
imputed_vec[names(cell_imputed_count)] <- log(cell_imputed_count)
all_data$imputed_count <- imputed_vec

# grDevices::colorRampPalette("lightgray", "blue")
p1 <- scCustomize::FeaturePlot_scCustom(all_data, 
                                        colors_use = list("lightgray", "blue"),
                                        na_cutoff = stats::quantile(all_data$imputed_count, probs = 0.05, na.rm = T),
                                        na_color = "bisque",
                                        reduction = "umap", 
                                        features = "imputed_count")
p1 <- p1 + ggplot2::ggtitle("DABTRAM Day10 imputed counts (Log-scale)\nAll RNA topics and some ATAC LSI")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6i/Writeup6i_DABTRAM-day10_imputed-lineage_umap_all-fasttopics-and-some-lsi.png"),
                p1, device = "png", width = 7, height = 5, units = "in")

lineage_imputed_count2 <- log10(lineage_imputed_count+1)
lineage_future_count2 <- log10(lineage_future_count+1)

labeling_vec <- rep(FALSE, length(lineage_imputed_count2))
labeling_vec[which(lineage_imputed_count2 >= 1.5)] <- TRUE
labeling_vec[which(lineage_future_count2 >= 1.5)] <- TRUE
table(labeling_vec)

df <- data.frame(lineage_imputed_count = lineage_imputed_count2,
                 lineage_future_count = lineage_future_count2,
                 name = names(lineage_imputed_count2),
                 labeling = labeling_vec)
# put all the labeling == TRUE on bottom
df <- df[c(which(!df[,"labeling"]), which(df[,"labeling"])),]

p1 <- ggplot2::ggplot(df, ggplot2::aes(x = lineage_future_count, y = lineage_imputed_count))
p1 <- p1 + ggplot2::geom_point(ggplot2::aes(color = labeling))
p1 <- p1 + ggplot2::scale_colour_manual(values=c("black", "red"))
p1 <- p1 + ggrepel::geom_text_repel(data = subset(df, labeling == TRUE),
                                    ggplot2::aes(label = name, color = labeling),
                                    box.padding = ggplot2::unit(0.5, 'lines'),
                                    point.padding = ggplot2::unit(1.6, 'lines'),
                                    max.overlaps = 50)
p1 <- p1 + ggplot2::ggtitle(paste0("DABTRAM Day10 imputed counts (Log-scale)\nAll RNA topics and some ATAC LSI, Corr: ", 
                                   round(stats::cor(lineage_imputed_count2, lineage_future_count2), 3))) +
  ggplot2::xlab("Observed lineage count (Log10)") + ggplot2::ylab("Predicted lineage count (Log10)")
p1 <- p1 + Seurat::NoLegend()

ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6i/Writeup6i_DABTRAM-day10_lineage-level_counts_all-fasttopics-and-some-lsi_labeled.png"),
                p1, device = "png", width = 10, height = 10, units = "in")

