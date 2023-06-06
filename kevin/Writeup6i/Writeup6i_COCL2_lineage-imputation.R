rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

Seurat::DefaultAssay(all_data) <- "RNA"
all_data[["geneActivity"]] <- NULL
all_data[["fasttopic_CIS"]] <- NULL
all_data[["fasttopic_DABTRAM"]] <- NULL
all_data[["common_tcca"]] <- NULL
all_data[["distinct1_tcca"]] <- NULL
all_data[["distinct2_tcca"]] <- NULL
all_data[["activity.umap"]] <- NULL

keep_vec <- rep(NA, ncol(all_data))
idx1 <- which(all_data$dataset %in% c("day0", "day10_COCL2", "week5_COCL2"))
idx2 <- which(all_data$assigned_posterior >= 0.5)
idx3 <- which(!is.na(all_data$assigned_lineage))
keep_vec[intersect(intersect(idx1, idx2), idx3)] <- TRUE
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == TRUE)

all_data[["umap"]] <- NULL
set.seed(10)
all_data <- Seurat::RunUMAP(all_data, 
                            reduction = "fasttopic_COCL2",
                            dims = 1:30)


tab_mat <- table(all_data$assigned_lineage, all_data$dataset)

tier3_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("week5_", treatment)] >= 50)]
tier2_lineages <- rownames(tab_mat)[intersect(which(tab_mat[,paste0("week5_", treatment)] >= 3),
                                              which(tab_mat[,paste0("week5_", treatment)] <= 49))]
tier1_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("week5_", treatment)] <= 2)]
length(tier1_lineages); length(tier2_lineages); length(tier3_lineages)

tier3_idx <- intersect(
  which(all_data$assigned_lineage %in% tier3_lineages),
  which(all_data$dataset == paste0("day10_", treatment))
)
tier2_idx <- intersect(
  which(all_data$assigned_lineage %in% tier2_lineages),
  which(all_data$dataset == paste0("day10_", treatment))
)
tier1_idx <- intersect(
  which(all_data$assigned_lineage %in% tier1_lineages),
  which(all_data$dataset == paste0("day10_", treatment))
)
keep_vec <- rep(NA, ncol(all_data))
keep_vec[tier1_idx] <- paste0("3high_winner_", treatment)
keep_vec[tier2_idx] <- paste0("2mid_winner_", treatment)
keep_vec[tier3_idx] <- paste0("1loser_", treatment)
table(keep_vec)
all_data$keep <- keep_vec
all_data2 <- subset(all_data, keep %in% c(paste0("3high_winner_", treatment),
                                          paste0("2mid_winner_", treatment),
                                          paste0("1loser_", treatment)))

save(all_data2, all_data,
     date_of_run, session_info, treatment,
     file = "../../../../out/kevin/Writeup6i/Writeup6i_COCL2-day10_extracted.RData")

###############

# now try using the somewhat correlated features of the topic model
topic_mat <- all_data2[["fasttopic_COCL2"]]@cell.embeddings
# fish out all the topics that are correlated with growth
all_data2$tier <- all_data2$keep

log10pval_vec <- sapply(1:ncol(topic_mat), function(j){
  x <- topic_mat[,j]
  y <- as.factor(all_data2$tier)
  res <- stats::oneway.test(x ~ y)
  -log10(res$p.value)
})
names(log10pval_vec) <- colnames(topic_mat)
round(quantile(log10pval_vec))

cell_features <- cbind(1, scale(topic_mat[,which(log10pval_vec >= 29)]))
p <- ncol(cell_features)
colnames(cell_features)[1] <- "Intercept"

cell_lineage <- all_data2$assigned_lineage
cell_lineage <- cell_lineage[!is.na(cell_lineage)]
cell_features <- cell_features[names(cell_lineage),]
uniq_lineage <- sort(unique(cell_lineage))
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
lineage_current_count <- tab_mat[uniq_lineage,"day10_COCL2"]
lineage_future_count <- tab_mat[uniq_lineage,"week5_COCL2"]

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
round(lineage_res$fit$coefficient_vec, 2)

cell_imputed_count <- as.numeric(exp(cell_features %*% lineage_res$fit$coefficient_vec))
names(cell_imputed_count) <- rownames(cell_features)
quantile(round(cell_imputed_count), probs = seq(0.9,1,length.out=11))

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
                                        na_cutoff = min(all_data$imputed_count, na.rm = T),
                                        na_color = "bisque",
                                        reduction = "umap", 
                                        features = "imputed_count")
p1 <- p1 + ggplot2::ggtitle("COCL2 Day10 imputed counts (Log-scale)")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6i/Writeup6i_COCL2-day10_imputed-lineage_umap.png"),
                p1, device = "png", width = 7, height = 5, units = "in")

# p1 <- Seurat::FeaturePlot(all_data, reduction = "umap", 
#                           features = "imputed_count")
# p1 <- p1 + ggplot2::ggtitle("COCL2 Day10 imputed counts (Log-scale)")
# ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6i/Writeup6i_COCL2-day10_imputed-lineage_umap2.png"),
#                 p1, device = "png", width = 7, height = 5, units = "in")


p1 <- Seurat::DimPlot(all_data, reduction = "umap",
                      group.by = "dataset")
p1 <- p1 + ggplot2::ggtitle("COCL2 (Only confident-assigned cells)")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6i/Writeup6i_COCL2_umap.png"),
                p1, device = "png", width = 7, height = 5, units = "in")


png("../../../../out/figures/Writeup6i/Writeup6i_COCL2-day10_lineage-level_counts.png",
    height = 1500, width = 1500, res = 300, units = "px")
plot(log10(lineage_future_count+1), 
     log10(lineage_imputed_count+1),
     xlab = "Observed lineage count", ylab = "Predicted lineage count",
     main = paste0("COCL2 Day10 imputed counts\n(Log-scale), Cor: ", 
                   round(stats::cor(log10(lineage_imputed_count+1), 
                                    log10(lineage_future_count+1)),2)),
     pch = 16, asp = T, col = rgb(0.5,0.5,0.5,0.5))
graphics.off()
