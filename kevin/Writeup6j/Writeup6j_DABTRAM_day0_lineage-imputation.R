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

treatment <- "DABTRAM"

Seurat::DefaultAssay(all_data) <- "RNA"
all_data[["geneActivity"]] <- NULL
all_data[["fasttopic_CIS"]] <- NULL
all_data[["fasttopic_COCL2"]] <- NULL
all_data[["common_tcca"]] <- NULL
all_data[["distinct1_tcca"]] <- NULL
all_data[["distinct2_tcca"]] <- NULL
all_data[["activity.umap"]] <- NULL

keep_vec <- rep(NA, ncol(all_data))
idx1 <- which(all_data$dataset %in% c("day0", "day10_DABTRAM", "week5_DABTRAM"))
idx2 <- which(all_data$assigned_posterior >= 0.5)
idx3 <- which(!is.na(all_data$assigned_lineage))
keep_vec[intersect(intersect(idx1, idx2), idx3)] <- TRUE
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == TRUE)

all_data[["umap"]] <- NULL
set.seed(10)
all_data <- Seurat::RunUMAP(all_data, 
                            reduction = "fasttopic_DABTRAM",
                            dims = 1:30)

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)

tier3_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("day10_", treatment)] >= 20)]
tier2_lineages <- rownames(tab_mat)[intersect(which(tab_mat[,paste0("day10_", treatment)] >= 3),
                                              which(tab_mat[,paste0("day10_", treatment)] <= 19))]
tier1_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("day10_", treatment)] <= 2)]
length(tier1_lineages); length(tier2_lineages); length(tier3_lineages)

tier3_idx <- intersect(
  which(all_data$assigned_lineage %in% tier3_lineages),
  which(all_data$dataset == paste0("day0"))
)
tier2_idx <- intersect(
  which(all_data$assigned_lineage %in% tier2_lineages),
  which(all_data$dataset == paste0("day0"))
)
tier1_idx <- intersect(
  which(all_data$assigned_lineage %in% tier1_lineages),
  which(all_data$dataset == paste0("day0"))
)
keep_vec <- rep(NA, ncol(all_data))
keep_vec[tier1_idx] <- paste0("1loser_", treatment)
keep_vec[tier2_idx] <- paste0("2mid_winner_", treatment)
keep_vec[tier3_idx] <- paste0("3high_winner_", treatment)
table(keep_vec)
all_data$keep <- keep_vec
all_data2 <- subset(all_data, keep %in% c(paste0("3high_winner_", treatment),
                                          paste0("2mid_winner_", treatment),
                                          paste0("1loser_", treatment)))

save(all_data2, all_data,
     date_of_run, session_info, treatment,
     file = "../../../../out/kevin/Writeup6j/Writeup6j_DABTRAM-day0_extracted.RData")

################################


cell_lineage <- all_data2$assigned_lineage
cell_lineage <- cell_lineage[!is.na(cell_lineage)]

# construct cell_features matrix
topic_mat <- all_data[["fasttopic_DABTRAM"]]@cell.embeddings[names(cell_lineage),]
atac_mat <- all_data2[["lsi"]]@cell.embeddings[names(cell_lineage),]

# let's try using all the topics #Clueless
cell_features <- cbind(1, scale(topic_mat), 
                       scale(atac_mat))
p <- ncol(cell_features)
colnames(cell_features)[1] <- "Intercept"

cell_features <- cell_features[names(cell_lineage),]
uniq_lineage <- sort(unique(cell_lineage))
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
lineage_current_count <- tab_mat[uniq_lineage,"day0"]
lineage_future_count <- tab_mat[uniq_lineage,"day10_DABTRAM"]

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
                                                 verbose = 1)

#############


round(lineage_res$fit$coefficient_vec, 2)
round(quantile(abs(lineage_res$fit$coefficient_vec)), 2)

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
                                        na_cutoff = quantile(all_data$imputed_count, probs = 0.05, na.rm = T),
                                        na_color = "bisque",
                                        reduction = "umap", 
                                        features = "imputed_count")
p1 <- p1 + ggplot2::ggtitle("COCL2 Day0 imputed counts (Log-scale)")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6j/Writeup6j_COCL2-day0_imputed-lineage_umap.png"),
                p1, device = "png", width = 7, height = 5, units = "in")

p1 <- Seurat::DimPlot(all_data, reduction = "umap",
                      group.by = "dataset")
p1 <- p1 + ggplot2::ggtitle("COCL2 (Only confident-assigned cells)")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6j/Writeup6j_COCL2_umap.png"),
                p1, device = "png", width = 7, height = 5, units = "in")

####

all(names(lineage_imputed_count) == names(lineage_future_count))

lineage_imputed_count2 <- log10(lineage_imputed_count+1)
lineage_future_count2 <- log10(lineage_future_count+1)

labeling_vec <- rep(FALSE, length(lineage_imputed_count2))
labeling_vec[which(lineage_imputed_count2 >= 1.25)] <- TRUE
labeling_vec[which(lineage_future_count2 >= 1.25)] <- TRUE
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
p1 <- p1 + ggplot2::ggtitle(paste0(
  "COCL2 Day0 imputed counts",
  "\n(RNA fasttopics, ATAC LSI), (Log-scale)",
  "\nCorrelation:", round(stats::cor(lineage_imputed_count2, lineage_future_count2), 2))
) +
  ggplot2::xlab("Observed lineage count (Log10)") + ggplot2::ylab("Predicted lineage count (Log10)")
p1 <- p1 + Seurat::NoLegend()

ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6j/Writeup6j_COCL2-day0_imputed-lineage_all-RNA-and-ATAC.png"),
                p1, device = "png", width = 10, height = 10, units = "in")


