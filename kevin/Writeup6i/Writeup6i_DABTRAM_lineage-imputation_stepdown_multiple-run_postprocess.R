# try data fission
rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6i/Writeup6i_DABTRAM_lineage-imputation_stepdown_multiple-run.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

# extract the average test and train
len_vec <- sapply(coefficient_list_list, function(lis){
  length(lis[[1]]$coefficient_vec)
})
len <- length(coefficient_list_list)
train_mat <- matrix(NA, nrow = 3, ncol = len)
test_mat <- matrix(NA, nrow = 3, ncol = len)
for(i in 1:len){
  train_tmp <- sapply(1:10, function(j){
    coefficient_list_list[[i]][[j]]$training_obj
  })
  test_tmp <- sapply(1:10, function(j){
    coefficient_list_list[[i]][[j]]$testing_obj
  })
  
  train_mat[,i] <- stats::quantile(train_tmp, probs = c(0.25,0.5,0.75))
  test_mat[,i] <- stats::quantile(test_tmp, probs = c(0.25,0.5,0.75))
}
rownames(train_mat) <- c("quant0.25", "quant0.5", "quant0.75")
rownames(test_mat) <- c("quant0.25", "quant0.5", "quant0.75")
colnames(train_mat) <- paste0("Num", len_vec)
colnames(test_mat) <- paste0("Num", len_vec)

png("../../../../out/figures/Writeup6i/Writeup6i_DABTRAM-day10_imputation3_stepdown-variable_loglikelihood_banded.png", 
    width = 1500, height = 1500, units = "px", res = 300)
ylim <- range(c(train_mat), c(test_mat))
x_vec <- len_vec
plot(NA, xlim = range(len_vec), ylim = ylim, 
     xlab = "# variables", ylab = "-Loglikelihood", 
     main = "Stepdown variable selection: DABTRAM\nBlack:Train, Red:Test")
polygon(x = c(len_vec, rev(len_vec)),
        y = c(train_mat[1,], rev(train_mat[3,])),
        col = rgb(0.5, 0.5, 0.5, 0.5),
        border = "black")
polygon(x = c(len_vec, rev(len_vec)),
        y = c(test_mat[1,], rev(test_mat[3,])),
        col = rgb(0.5, 0, 0, 0.5),
        border = "red")
points(x = len_vec, y = train_mat[2,], pch = 16, cex = 0.5)
lines(x = len_vec, y = train_mat[2,], lwd = 2, lty = 2)
points(x = len_vec, y = test_mat[2,], pch = 16, cex = 0.5, col = "red")
lines(x = len_vec, y = test_mat[2,], lwd = 2, lty = 2, col = "red")
graphics.off()

####################################################

# let's try fitting the model with the selected model

library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6h/Writeup6h_DABTRAM-day10_pseudotime.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

selected_vars <- names(coefficient_list_list[[which(colnames(test_mat) == "Num5")]][[1]]$coefficient_vec)

# > selected_vars
# [1] "Intercept"          "fastTopicDABTRAM_2" "LSI_2"             
# [4] "LSI_24"             "Pseudotime"

fasttopic_mat <- all_data2[["fasttopic_DABTRAM"]]@cell.embeddings
lsi_mat <- all_data2[["lsi"]]@cell.embeddings
pseudotime_vec <- all_data$pseudotime[colnames(all_data2)]

cell_features <- cbind(1, 
                       scale(fasttopic_mat[,2,drop=F]), 
                       scale(lsi_mat[,c(2,24)]), 
                       scale(pseudotime_vec))
colnames(cell_features)[1] <- "Intercept"
colnames(cell_features)[5] <- "Pseudotime"
cell_lineage <- all_data2$assigned_lineage
cell_lineage <- cell_lineage[!is.na(cell_lineage)]
cell_features <- cell_features[names(cell_lineage),]
uniq_lineage <- sort(unique(cell_lineage))
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
lineage_current_count <- tab_mat[uniq_lineage,"day10_DABTRAM"]
lineage_future_count <- tab_mat[uniq_lineage,"week5_DABTRAM"]

round(stats::cor(cell_features[,-1]),2)

p <- ncol(cell_features)
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

######

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
                                        na_cutoff = quantile(all_data$imputed_count, probs = 0.05, na.rm = T),
                                        na_color = "bisque",
                                        reduction = "umap", 
                                        features = "imputed_count")
p1 <- p1 + ggplot2::ggtitle("DABTRAM Day10 imputed counts\n(Stepdown from RNA fasttopics, ATAC LSI, Pseudotime), (Log-scale)")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6i/Writeup6i_DABTRAM-day10_imputed-lineage_stepdown_umap.png"),
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
  "DABTRAM Day10 imputed counts",
  "\n(Stepdown from RNA fasttopics, ATAC LSI, Pseudotime), (Log-scale)",
  "\nCorrelation:", round(stats::cor(lineage_imputed_count2, lineage_future_count2), 2))
) +
  ggplot2::xlab("Observed lineage count (Log10)") + ggplot2::ylab("Predicted lineage count (Log10)")
p1 <- p1 + Seurat::NoLegend()

ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6i/Writeup6i_DABTRAM-day10_imputed-lineage_stepdown_lineage-level_counts_labeled.png"),
                p1, device = "png", width = 10, height = 10, units = "in")



