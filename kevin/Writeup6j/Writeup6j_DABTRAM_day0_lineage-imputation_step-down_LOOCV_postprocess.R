# try data fission
rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6j/Writeup6j_DABTRAM_day0_lineage-imputation_step-down_LOOCV.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

# extract the average test and train
len_vec <- sapply(coefficient_list_list, function(lis){
  length(lis$fit$fit$coefficient_vec)
})
len <- length(coefficient_list_list)
train_vec <- sapply(coefficient_list_list, function(x){
  x$fit$fit$objective_val
})
test_vec <- sapply(1:len, function(i){
  if(i == 1) return(NA)
  removed_var <- setdiff(colnames(coefficient_list_list[[i]]$loocv_mat),
                         names(coefficient_list_list[[i]]$fit$fit$coefficient_vec))
  mean(coefficient_list_list[[i]]$loocv_mat[,removed_var])
})

png("../../../../out/figures/Writeup6j/Writeup6j_DABTRAM-day0_imputation_stepdown-LOOCV_loglikelihood.png", 
    width = 3000, height = 1500, units = "px", res = 300)
par(mfrow = c(1,2))
ylim <- range(train_vec)
x_vec <- len_vec
plot(NA, xlim = range(len_vec), ylim = ylim, 
     xlab = "# variables", ylab = "-Loglikelihood", 
     main = "Stepdown variable selection: DABTRAM\n(Training)")
points(x = len_vec, y = train_vec, pch = 16, cex = 0.5)
lines(x = len_vec, y = train_vec, lwd = 2, lty = 2)

ylim <- range(test_vec, na.rm = T)
plot(NA, xlim = range(len_vec), ylim = ylim, 
     xlab = "# variables", ylab = "-Loglikelihood", 
     main = "Stepdown variable selection: DABTRAM\n(Testing)")
points(x = len_vec, y = test_vec, pch = 16, cex = 0.5, col = "red")
lines(x = len_vec, y = test_vec, lwd = 2, lty = 2, col = "red")

graphics.off()


####################################################

# let's try fitting the model with the selected model

load("../../../../out/kevin/Writeup6j/Writeup6j_DABTRAM-day0_extracted.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

all_data[["umap"]] <- NULL
set.seed(10)
all_data <- Seurat::RunUMAP(all_data, 
                            reduction = "fasttopic_DABTRAM",
                            dims = 1:30)

idx <- which.min(test_vec)
lineage_res <- coefficient_list_list[[idx]]$fit
var_names <- names(lineage_res$fit$coefficient_vec)

fasttopic_mat <- all_data2[["fasttopic_DABTRAM"]]@cell.embeddings
lsi_mat <- all_data2[["lsi"]]@cell.embeddings

cell_features <- cbind(1, 
                       scale(fasttopic_mat), 
                       scale(lsi_mat))
colnames(cell_features)[1] <- "Intercept"
cell_features <- cell_features[,var_names]
cell_lineage <- all_data2$assigned_lineage
cell_lineage <- cell_lineage[!is.na(cell_lineage)]
cell_features <- cell_features[names(cell_lineage),]
uniq_lineage <- sort(unique(cell_lineage))
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
lineage_current_count <- tab_mat[uniq_lineage,"day0"]
lineage_future_count <- tab_mat[uniq_lineage,"day10_DABTRAM"]

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
                                        colors_use = list("red", "lightgray", "blue"),
                                        na_cutoff = quantile(all_data$imputed_count, probs = 0.05, na.rm = T),
                                        na_color = "bisque",
                                        reduction = "umap", 
                                        features = "imputed_count")
p1 <- p1 + ggplot2::ggtitle("DABTRAM Day0 imputed counts\n(Stepdown from RNA fasttopics, ATAC LSI), (Log-scale)")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6j/Writeup6j_DABTRAM-day0_imputation_stepdown-LOOCV_umap.png"),
                p1, device = "png", width = 7, height = 5, units = "in")

###

all(names(lineage_imputed_count) == names(lineage_future_count))

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
p1 <- p1 + ggplot2::ggtitle(paste0(
  "DABTRAM Day0 imputed counts",
  "\n(Stepdown from RNA fasttopics, ATAC LSI), (Log-scale)",
  "\nCorrelation:", round(stats::cor(lineage_imputed_count2, lineage_future_count2), 2), " with ", length(var_names)-1, " non-Intercept variables")
) +
  ggplot2::xlab("Observed lineage count (Log10)") + ggplot2::ylab("Predicted lineage count (Log10)")
p1 <- p1 + Seurat::NoLegend()

ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6j/Writeup6j_DABTRAM-day0_imputation_stepdown-LOOCV_lineage-level_counts_labeled.png"),
                p1, device = "png", width = 10, height = 10, units = "in")


