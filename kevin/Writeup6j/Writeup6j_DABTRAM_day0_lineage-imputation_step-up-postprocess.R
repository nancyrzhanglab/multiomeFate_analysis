rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6j/Writeup6j_DABTRAM_day0_lineage-imputation_step-up_step2.RData")
load("../../../../out/kevin/Writeup6j/Writeup6j_DABTRAM-day0_extracted.RData")
all_data2$tier_vec <- all_data2$keep

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

cell_lineage <- all_data2$assigned_lineage
cell_lineage <- cell_lineage[!is.na(cell_lineage)]

# construct cell_features matrix
topic_mat <- all_data[["fasttopic_DABTRAM"]]@cell.embeddings[names(cell_lineage),]
atac_mat <- all_data2[["lsi"]]@cell.embeddings[names(cell_lineage),]

# let's try using all the topics 
cell_features_full <- cbind(1, scale(topic_mat), 
                            scale(atac_mat))
p <- ncol(cell_features_full)
colnames(cell_features_full)[1] <- "Intercept"

cell_features_full <- cell_features_full[names(cell_lineage),]
uniq_lineage <- sort(unique(cell_lineage))
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
lineage_current_count <- tab_mat[uniq_lineage,"day0"]
lineage_future_count <- tab_mat[uniq_lineage,"day10_DABTRAM"]

##########################

# fit each training 
training_list <- vector("list", length = length(coefficient_list_list))
for(i in 1:length(training_list)){
  print(paste0("Working on Model ", i, " out of ", length(training_list)))
  var_names <- coefficient_list_list[[i]]$var_next_iteration
  
  cell_features <- cell_features_full[,var_names,drop=F]
  p <- ncol(cell_features)
  tmp <- quantile(abs(cell_features[,-which(colnames(cell_features) == "Intercept"),drop=F]), probs = 0.95)
  coef_val <- 2*log(max(lineage_future_count/lineage_current_count))/((p-1)*tmp)
  coefficient_initial <- pmin(rep(coef_val, p), 5)
  names(coefficient_initial) <- colnames(cell_features)
  coefficient_initial["Intercept"] <- 0
  
  set.seed(10)
  tmp_res <- multiomeFate:::lineage_imputation(cell_features = cell_features,
                                               cell_lineage = cell_lineage,
                                               coefficient_initial = coefficient_initial,
                                               lineage_future_count = lineage_future_count,
                                               random_initializations = 10,
                                               verbose = 0)
  
  print(tmp_res$fit$objective_val)
  training_list[[i]] <- tmp_res
}


save(date_of_run, session_info,
     training_list,
     file = "../../../../out/kevin/Writeup6j/Writeup6j_DABTRAM_day0_lineage-imputation_stepwise-up_training.RData")

##########################

# load("../../../../out/kevin/Writeup6j/Writeup6j_DABTRAM_day0_lineage-imputation_stepwise-up_training.RData")

# extract the average test and train
len_vec <- sapply(coefficient_list_list, function(lis){
  length(lis$var_next_iteration)
})
len <- length(coefficient_list_list)
train_vec <- sapply(training_list, function(x){
  x$fit$objective_val
})
test_vec <- colMeans(loocv_mat)

png("../../../../out/figures/Writeup6j/Writeup6j_DABTRAM-day0_imputation_stepup_loglikelihood.png", 
    width = 3000, height = 1500, units = "px", res = 300)
par(mfrow = c(1,2))
ylim <- range(train_vec)
x_vec <- len_vec
plot(NA, xlim = range(len_vec), ylim = ylim, 
     xlab = "# variables", ylab = "-Loglikelihood", 
     main = "Stepup variable selection: DABTRAM\n(Day0->Day10, Stepup, Training)")
points(x = len_vec, y = train_vec, pch = 16, cex = 0.5)
lines(x = len_vec, y = train_vec, lwd = 2, lty = 2)

ylim <- range(test_vec, na.rm = T)
plot(NA, xlim = range(len_vec), ylim = ylim, 
     xlab = "# variables", ylab = "-Loglikelihood", 
     main = "Stepup variable selection: DABTRAM\n(Day0->Day10, Stepup, Testing)")
points(x = len_vec, y = test_vec, pch = 16, cex = 0.5, col = "red")
lines(x = len_vec, y = test_vec, lwd = 2, lty = 2, col = "red")

graphics.off()


#######################

all_data[["umap"]] <- NULL
set.seed(10)
all_data <- Seurat::RunUMAP(all_data, 
                            reduction = "fasttopic_DABTRAM",
                            dims = 1:30)

idx <- which.min(test_vec)
lineage_res <- training_list[[idx]]$fit
var_names <- names(lineage_res$coefficient_vec)

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

cell_imputed_count <- as.numeric(exp(cell_features %*% lineage_res$coefficient_vec))
names(cell_imputed_count) <- rownames(cell_features)
quantile(round(cell_imputed_count), probs = seq(0.9,1,length.out=11))
round(lineage_res$coefficient_vec, 3)

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
p1 <- p1 + ggplot2::ggtitle("DABTRAM Day0 imputed counts\n(Stepup from RNA fasttopics, ATAC LSI), (Log-scale)")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6j/Writeup6j_DABTRAM-day0_imputation_stepup_umap.png"),
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
  "\n(Stepup from RNA fasttopics, ATAC LSI), (Log-scale)",
  "\nCorrelation:", round(stats::cor(lineage_imputed_count2, lineage_future_count2), 2), " with ", length(var_names)-1, " non-Intercept variables")
) +
  ggplot2::xlab("Observed lineage count (Log10)") + ggplot2::ylab("Predicted lineage count (Log10)")
p1 <- p1 + Seurat::NoLegend()

ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6j/Writeup6j_DABTRAM-day0_imputation_stepup_lineage-level_counts_labeled.png"),
                p1, device = "png", width = 10, height = 10, units = "in")

