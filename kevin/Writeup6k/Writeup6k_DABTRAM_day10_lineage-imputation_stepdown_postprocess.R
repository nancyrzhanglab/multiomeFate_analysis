rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6k/Writeup6k_DABTRAM_day10_lineage-imputation_stepdown_step2.RData")
load("../../../../out/kevin/Writeup6h/Writeup6h_DABTRAM-day10_pseudotime.RData")
all_data2$keep <- !is.na(all_data2$assigned_lineage)
all_data2 <- subset(all_data2, keep == TRUE)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

# extract the average test and train
len_vec <- sapply(coefficient_list_list, function(lis){
  length(lis$var_next_iteration)
})
len <- length(coefficient_list_list)
train_vec <- sapply(coefficient_list_list, function(x){
  x$fit$objective_val
})
test_vec <- colMeans(loocv_mat)

png("../../../../out/figures/Writeup6k/Writeup6k_DABTRAM-day10_imputation_stepdown_loglikelihood.png", 
    width = 3000, height = 1500, units = "px", res = 300)
par(mfrow = c(1,2))
ylim <- range(train_vec)
x_vec <- len_vec
plot(NA, xlim = range(len_vec), ylim = ylim, 
     xlab = "# variables", ylab = "-Loglikelihood", 
     main = "Stepdown variable selection: DABTRAM\n(Day10->Week5, Stepdown, Training)")
points(x = len_vec, y = train_vec, pch = 16, cex = 0.5)
lines(x = len_vec, y = train_vec, lwd = 2, lty = 2)

ylim <- range(test_vec, na.rm = T)
plot(NA, xlim = range(len_vec), ylim = ylim, 
     xlab = "# variables", ylab = "-Loglikelihood", 
     main = "Stepdown variable selection: DABTRAM\n(Day10->Week5, Stepdown, Testing)")
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
test_sd_val <- mean(matrixStats::rowSds(loocv_mat[,pmax(1, idx-2):pmin(ncol(loocv_mat), idx+2)]))
test_idx_vec <- which(test_vec <= test_vec[idx] + test_sd_val)

#######################

fasttopic_mat <- all_data2[["fasttopic_DABTRAM"]]@cell.embeddings
lsi_mat <- all_data2[["lsi"]]@cell.embeddings

cell_features <- cbind(1, 
                       scale(fasttopic_mat), 
                       scale(lsi_mat))
colnames(cell_features)[1] <- "Intercept"
cell_lineage <- all_data2$assigned_lineage
cell_lineage <- cell_lineage[!is.na(cell_lineage)]
cell_features <- cell_features[names(cell_lineage),]
uniq_lineage <- sort(unique(cell_lineage))
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
lineage_future_count <- tab_mat[uniq_lineage,"week5_DABTRAM"]

# now to select the best model among all these reasonably good models:
cor_list <- lapply(test_idx_vec, function(test_idx){
  lineage_res <- coefficient_list_list[[test_idx]]$fit
  var_names <- names(lineage_res$coefficient_vec)
  
  fasttopic_mat <- all_data2[["fasttopic_DABTRAM"]]@cell.embeddings
  lsi_mat <- all_data2[["lsi"]]@cell.embeddings
  
  cell_features <- cell_features[,var_names]
  cell_imputed_count <- as.numeric(exp(cell_features %*% lineage_res$coefficient_vec))
  names(cell_imputed_count) <- rownames(cell_features)
  lineage_imputed_count <- sapply(uniq_lineage, function(lineage){
    sum(cell_imputed_count[which(cell_lineage == lineage)])
  })
  
  imputed_vec <- rep(NA, ncol(all_data))
  names(imputed_vec) <- colnames(all_data)
  imputed_vec[names(cell_imputed_count)] <- log10(cell_imputed_count)
  cor_val <- stats::cor(log10(lineage_imputed_count+1), log10(lineage_future_count+1))
  
  list(lineage_imputed_count = lineage_imputed_count,
       imputed_vec = imputed_vec,
       cor_val = cor_val)
})
cor_vec <- sapply(cor_list, function(x){x$cor_val})

all_data$imputed_count <- cor_list[[which.max(cor_vec)]]$imputed_vec
fit <- coefficient_list_list[[test_idx_vec[which.max(cor_vec)]]]$fit
var_names <- names(fit$coefficient_vec)
lineage_imputed_count <- cor_list[[which.max(cor_vec)]]$lineage_imputed_count

# grDevices::colorRampPalette("lightgray", "blue")
p1 <- scCustomize::FeaturePlot_scCustom(all_data, 
                                        colors_use = list("red", "lightgray", "blue"),
                                        na_cutoff = quantile(all_data$imputed_count, probs = 0.05, na.rm = T),
                                        na_color = "bisque",
                                        reduction = "umap", 
                                        features = "imputed_count")
p1 <- p1 + ggplot2::ggtitle("DABTRAM Day10 imputed counts\n(Stepdown from RNA fasttopics, ATAC LSI), (Log-scale)")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6k/Writeup6k_DABTRAM-day10_imputation_stepdown_umap.png"),
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
  "\n(Stepdown from RNA fasttopics, ATAC LSI), (Log-scale)",
  "\nCorrelation:", round(stats::cor(lineage_imputed_count2, lineage_future_count2), 2), " with ", length(var_names)-1, " non-Intercept variables")
) +
  ggplot2::xlab("Observed lineage count (Log10)") + ggplot2::ylab("Predicted lineage count (Log10)")
p1 <- p1 + Seurat::NoLegend()

ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6k/Writeup6k_DABTRAM-day10_imputation_stepdown_lineage-level_counts_labeled.png"),
                p1, device = "png", width = 10, height = 10, units = "in")

####################

save(all_data, all_data2, fit, lineage_future_count,
     date_of_run, session_info, 
     file = "../../../../out/kevin/Writeup6k/Writeup6k_DABTRAM_day10_lineage-imputation_stepdown_postprocessed.RData")
