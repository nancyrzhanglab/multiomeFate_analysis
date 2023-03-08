rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")
load("../../../../out/kevin/Writeup6c/Writeup6c_peak_counter_5000.RData")

countmat_nopeak <- Matrix::Matrix(countmat_nopeak, sparse = T)
quantile(countmat_nopeak@x)

lin_mat <- table(all_data$assigned_lineage, all_data$dataset)
DABTRAM_lineage <- rownames(lin_mat)[which(lin_mat[,"day10_DABTRAM"] >= 20)]
day0_survive_idx <- intersect(which(all_data$assigned_lineage %in% DABTRAM_lineage),
                              which(all_data$dataset == "day0"))

day0_die_idx <- intersect(which(!all_data$assigned_lineage %in% DABTRAM_lineage),
                          which(all_data$dataset == "day0"))

p <- ncol(countmat_nopeak)
pval_vec <- sapply(1:p, function(j){
  if(j %% floor(p/10) == 0) cat('*')
  
  vec_1 <- as.numeric(countmat_nopeak[day0_survive_idx,j])
  vec_2 <- as.numeric(countmat_nopeak[day0_die_idx,j])
  if(diff(range(vec_1)) == 0 & diff(range(vec_2)) == 0) return(1)
  
  set.seed(10)
  tmp <- stats::wilcox.test(x = vec_1,
                            y = vec_2)
  tmp$p.value
})
names(pval_vec) <- colnames(countmat_nopeak)
# nan_idx <- which(is.nan(pval_vec))

logpval_vec <- -log10(pval_vec)
sort(names(logpval_vec)[which(logpval_vec >= 3)])

source("../Writeup6b/gene_list.R")
gene_vec <- sort(unique(unlist(keygenes)))
logpval_vec[intersect(gene_vec, names(logpval_vec))]

length(which(logpval_vec >= 1))/length(logpval_vec)
length(which(logpval_vec >= 2))/length(logpval_vec)

# make the ggplot with repealing text
gene_idx <- which(names(logpval_vec) %in% gene_vec)
idx <- rank(logpval_vec)
labeling_vec <- rep(0, length(logpval_vec))
labeling_vec[gene_idx] <- 1
text_vec <- names(logpval_vec)

logpval_vec <- c(logpval_vec[-gene_idx], logpval_vec[gene_idx])
labeling_vec <- c(labeling_vec[-gene_idx], labeling_vec[gene_idx])
labeling_vec <- as.factor(labeling_vec)
text_vec <- c(text_vec[-gene_idx], text_vec[gene_idx])

df <- data.frame(logpval = logpval_vec, Gene_rank = idx, 
                 labeling = labeling_vec,
                 text = text_vec)

plot1 <- ggplot2::ggplot(df, ggplot2::aes(x = Gene_rank, y = logpval)) 
plot1 <- plot1 + ggplot2::geom_point(ggplot2::aes(color = labeling))
plot1 <- plot1 + ggplot2::scale_colour_manual(values=c("black", "red"))
plot1 <- plot1 + ggrepel::geom_text_repel(data = subset(df, labeling == 1), 
                                          ggplot2::aes(label = text, color = labeling),
                                          box.padding = ggplot2::unit(0.5, 'lines'),
                                          size = 2,
                                          nudge_y = 1,
                                          max.overlaps = Inf,
                                          min.segment.length = 0)
plot1 <- plot1 + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6c/Writeup6c_peak_counter_all-genes.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")

## plot the differentiability for specific genes
genes <- c("EGFR", "MITF", "SOX10")
for(gene in genes){
  j <- which(colnames(countmat_nopeak) == gene)
  vec_1 <- as.numeric(countmat_nopeak[day0_survive_idx,j])
  vec_2 <- as.numeric(countmat_nopeak[day0_die_idx,j])
  
  range_vec <- c(0, quantile(c(vec_1, vec_2), probs = 0.99))
  if(range_vec[2] <= 15 & max(c(vec_1, vec_2)) >= 15){
    range_vec[2] <- (15+max(c(vec_1, vec_2)))/2
  }
  break_vec <- seq(min(range_vec), max(range_vec), length.out = 21)
  vec_1_tmp <- vec_1[vec_1 <= range_vec[2]]
  vec_2_tmp <- vec_2[vec_2 <= range_vec[2]]
  
  png(paste0("../../../../out/figures/Writeup6c/Writeup6c_peak_counter_DABTRAM_day0_gene-", gene,".png"),
      height = 1500, width = 3000, units = "px", res = 500)
  par(mfrow = c(1,2), mar = c(4,4,4,0.5))
  
  tmp <- hist(vec_1_tmp, breaks = break_vec, plot = F)
  # tmp$counts <- log10(tmp$counts+1)
  plot(tmp, col = "gray", 
       main = paste0("DABTRAM, for gene ", gene,
                     "\n% non-zero: ", round(100*length(which(vec_1 != 0))/length(vec_1)),
                     ", Wilcox -Log10 pval: ", round(logpval_vec[gene],3)),
       ylab = "Frequency of non-peak counts", 
       xlab = "# Day0 cells in lineages that survive by Day10",
       cex.main = 0.6, cex.lab = 0.6)
  mean_val <- mean(vec_1)
  median_val <- median(vec_1)
  lines(rep(mean_val, 2), c(-1e6,1e6), col = 2, lwd = 2)
  lines(rep(median_val, 2), c(-1e6,1e6), col = 3, lwd = 2, lty = 2)
  
  tmp <- hist(vec_2_tmp, breaks = break_vec, plot = F)
  # tmp$counts <- log10(tmp$counts+1)
  plot(tmp, col = "gray", 
       main = paste0("DABTRAM, for gene ", gene,
                     "\n% non-zero: ", round(100*length(which(vec_2 != 0))/length(vec_2))),
       ylab = "Frequency of non-peak counts", 
       xlab = "# Day0 cells in lineages that don't survive by Day10",
       cex.main = 0.6, cex.lab = 0.6)
  mean_val <- mean(vec_2)
  median_val <- median(vec_2)
  lines(rep(mean_val, 2), c(-1e6,1e6), col = 2, lwd = 2)
  lines(rep(median_val, 2), c(-1e6,1e6), col = 3, lwd = 2, lty = 2)
  
  graphics.off()
}

###################################

dim(countmat_nopeak)
quantile(as.numeric(countmat_nopeak[day0_survive_idx,]))
quantile(as.numeric(countmat_nopeak[day0_die_idx,]))

tmp <- countmat_nopeak[day0_survive_idx,]
quantile(tmp@x, probs = seq(0,1,length.out=11)); length(tmp@x)/prod(dim(tmp))

tmp <- countmat_nopeak[day0_die_idx,]
quantile(tmp@x, probs = seq(0,1,length.out=11)); length(tmp@x)/prod(dim(tmp))

###################################

# count the number of genes with non-peak counts

.nonzero_col <- function(mat, col_idx, bool_value){
  stopifnot(inherits(mat, "dgCMatrix"), col_idx %% 1 == 0,
            col_idx > 0, col_idx <= ncol(mat))
  
  val1 <- mat@p[col_idx]
  val2 <- mat@p[col_idx+1]
  
  if(val1 == val2) return(numeric(0))
  if(bool_value){
    # return the value
    mat@x[(val1+1):val2]
  } else {
    # return the row index
    mat@i[(val1+1):val2]+1
  }
}

countmat_nopeak_t <- Matrix::t(countmat_nopeak)
nonzero_survive <- sapply(day0_survive_idx, function(i){
  length(.nonzero_col(countmat_nopeak_t, i, bool_value = F))
})
nonzero_die <- sapply(day0_die_idx, function(i){
  length(.nonzero_col(countmat_nopeak_t, i, bool_value = F))
})

quantile(nonzero_survive, probs = seq(0,1,length.out=11))
quantile(nonzero_die, probs = seq(0,1,length.out=11))

set.seed(10)
tmp <- stats::wilcox.test(x = nonzero_survive,
                          y = nonzero_die)

range_vec <- c(0, quantile(c(nonzero_survive, nonzero_die), probs = 0.99))
break_vec <- seq(min(range_vec), max(range_vec), length.out = 21)
vec_1 <- nonzero_survive[nonzero_survive <= range_vec[2]]
vec_2 <- nonzero_die[nonzero_die <= range_vec[2]]

png(paste0("../../../../out/figures/Writeup6c/Writeup6c_peak_counter_DABTRAM_day0_nonzero-counts.png"),
    height = 1500, width = 3000, units = "px", res = 500)
par(mfrow = c(1,2), mar = c(4,4,4,0.5))

hist(vec_1, breaks = break_vec, col = "gray", 
     main = paste0("DABTRAM, surviving cells", 
                   "\nWilcox -Log10 pval: ", round(-log10(tmp$p.value),3)),
     ylab = "Frequency", 
     xlab = "# genes with non-peak counts, across all surviving cells",
     cex.main = 0.6, cex.lab = 0.6)
mean_val <- mean(nonzero_survive)
median_val <- median(nonzero_survive)
lines(rep(mean_val, 2), c(-1e6,1e6), col = 2, lwd = 2)
lines(rep(median_val, 2), c(-1e6,1e6), col = 3, lwd = 2, lty = 2)

hist(vec_2, breaks = break_vec, col = "gray", 
     main = paste0("DABTRAM, dying cells"),
     ylab = "Frequency", 
     xlab = "# genes with non-peak counts, across all dying cells",
     cex.main = 0.6, cex.lab = 0.6)
mean_val <- mean(nonzero_die)
median_val <- median(nonzero_die)
lines(rep(mean_val, 2), c(-1e6,1e6), col = 2, lwd = 2)
lines(rep(median_val, 2), c(-1e6,1e6), col = 3, lwd = 2, lty = 2)

graphics.off()


###################################

# look at sum of all ATAC peaks

mat <- all_data[["ATAC"]]@counts
mat_t <- Matrix::t(mat)
nonzero_survive2 <- sapply(day0_survive_idx, function(i){
  sum(.nonzero_col(mat_t, i, bool_value = T))
})
nonzero_die2 <- sapply(day0_die_idx, function(i){
  sum(.nonzero_col(mat_t, i, bool_value = T))
})

quantile(nonzero_survive2, probs = seq(0,1,length.out=11))
quantile(nonzero_die2, probs = seq(0,1,length.out=11))


set.seed(10)
tmp <- stats::wilcox.test(x = nonzero_survive2,
                          y = nonzero_die2,
                          alternative = "greater")

range_vec <- c(0, quantile(c(nonzero_survive2, nonzero_die2), probs = 0.9))
break_vec <- seq(min(range_vec), max(range_vec), length.out = 21)
vec_1 <- nonzero_survive2[nonzero_survive2 <= range_vec[2]]
vec_2 <- nonzero_die2[nonzero_die2 <= range_vec[2]]

png(paste0("../../../../out/figures/Writeup6c/Writeup6c_peak_counter_DABTRAM_day0_nonzero-counts_allpeaks.png"),
    height = 1500, width = 3000, units = "px", res = 500)
par(mfrow = c(1,2), mar = c(4,4,4,0.5))

hist(vec_1, breaks = break_vec, col = "gray", 
     main = paste0("DABTRAM, surviving cells", 
                   "\nOne-sided Wilcox -Log10 pval: ", round(-log10(tmp$p.value),3)),
     ylab = "Frequency", 
     xlab = "# genes with peak counts, across all surviving cells",
     cex.main = 0.6, cex.lab = 0.6)
mean_val <- mean(nonzero_survive2)
median_val <- median(nonzero_survive2)
lines(rep(mean_val, 2), c(-1e6,1e6), col = 2, lwd = 2)
lines(rep(median_val, 2), c(-1e6,1e6), col = 3, lwd = 2, lty = 2)

hist(vec_2, breaks = break_vec, col = "gray", 
     main = paste0("DABTRAM, dying cells"),
     ylab = "Frequency", 
     xlab = "# genes with peak counts, across all dying cells",
     cex.main = 0.6, cex.lab = 0.6)
mean_val <- mean(nonzero_die2)
median_val <- median(nonzero_die2)
lines(rep(mean_val, 2), c(-1e6,1e6), col = 2, lwd = 2)
lines(rep(median_val, 2), c(-1e6,1e6), col = 3, lwd = 2, lty = 2)

graphics.off()


###################################

# sum the number of non-peak counts

countmat_nopeak_t <- Matrix::t(countmat_nopeak)
nonzero_survive <- sapply(day0_survive_idx, function(i){
  sum(.nonzero_col(countmat_nopeak_t, i, bool_value = T))
})
nonzero_die <- sapply(day0_die_idx, function(i){
  sum(.nonzero_col(countmat_nopeak_t, i, bool_value = T))
})

round(quantile(nonzero_survive, probs = seq(0,1,length.out=11)))
round(quantile(nonzero_die, probs = seq(0,1,length.out=11)))

mean(nonzero_survive); sd(nonzero_survive)
mean(nonzero_die); sd(nonzero_die)

set.seed(10)
tmp <- stats::wilcox.test(x = nonzero_survive,
                          y = nonzero_die)

range_vec <- c(0, quantile(c(nonzero_survive, nonzero_die), probs = 0.99))
break_vec <- seq(min(range_vec), max(range_vec), length.out = 21)
vec_1 <- nonzero_survive[nonzero_survive <= range_vec[2]]
vec_2 <- nonzero_die[nonzero_die <= range_vec[2]]

png(paste0("../../../../out/figures/Writeup6c/Writeup6c_peak_counter_DABTRAM_day0_nonzero-counts-sum.png"),
    height = 1500, width = 3000, units = "px", res = 500)
par(mfrow = c(1,2), mar = c(4,4,4,0.5))

hist(vec_1, breaks = break_vec, col = "gray", 
     main = paste0("DABTRAM, surviving cells", 
                   "\nWilcox -Log10 pval: ", round(-log10(tmp$p.value),3)),
     ylab = "Frequency", 
     xlab = "sum of non-peak counts, across all surviving cells",
     cex.main = 0.6, cex.lab = 0.4)
mean_val <- mean(nonzero_survive)
median_val <- median(nonzero_survive)
lines(rep(mean_val, 2), c(-1e6,1e6), col = 2, lwd = 2)
lines(rep(median_val, 2), c(-1e6,1e6), col = 3, lwd = 2, lty = 2)

hist(vec_2, breaks = break_vec, col = "gray", 
     main = paste0("DABTRAM, dying cells"),
     ylab = "Frequency", 
     xlab = "sum of non-peak counts, across all dying cells",
     cex.main = 0.6, cex.lab = 0.4)
mean_val <- mean(nonzero_die)
median_val <- median(nonzero_die)
lines(rep(mean_val, 2), c(-1e6,1e6), col = 2, lwd = 2)
lines(rep(median_val, 2), c(-1e6,1e6), col = 3, lwd = 2, lty = 2)

graphics.off()

#############################################
#############################################

# try some regressions
# first, w/o NN smoothing
cell_idx <- intersect(intersect(which(all_data$dataset == "day0"),
                                which(all_data$assigned_posterior > 0.5)),
                      which(!is.na(all_data$assigned_lineage)))
fitness_vec <- sapply(cell_idx, function(i){
  sum(.nonzero_col(countmat_nopeak_t, i, bool_value = T))
})

lineage_vec <- all_data$assigned_lineage[cell_idx]
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
future_num_vec <- sapply(lineage_vec, function(lineage){
  tab_mat[lineage,"day10_DABTRAM"]
})

order_vec <- order(future_num_vec, decreasing = F)
x_vec <- log1p(future_num_vec[order_vec])
y_vec <- fitness_vec[order_vec]
tmp_df <- data.frame(y = y_vec, x = x_vec)
reg_res <- npregfast::frfast(y ~ x, data = tmp_df)
y_vec_pred <- stats::predict(object = reg_res, newdata = data.frame(x = x_vec))$Estimation[,"Pred"]

# compute nearest-neighbors, just among the valid day0 cells
set.seed(10)
mat <- all_data[["pca"]]@cell.embeddings[cell_idx,]
num_neigh <- 10
nn_mat <- RANN::nn2(mat, k = num_neigh+1)$nn.idx

# form smoothing nearest-neighbor matrix
n <- nrow(nn_mat)
i_vec <- rep(1:n, each = ncol(nn_mat))
j_vec <- unlist(lapply(1:n, function(i){nn_mat[i,]}))
avg_mat <- Matrix::sparseMatrix(i = i_vec, 
                                j = j_vec, 
                                x = rep(1/ncol(nn_mat), length(i_vec)), 
                                dims = c(n,n))
avg_mat <- Matrix::t(avg_mat)
future_num_vec_smoothed <- as.numeric(future_num_vec %*% avg_mat)

order_vec <- order(future_num_vec_smoothed, decreasing = F)
x_vec2 <- log1p(future_num_vec_smoothed[order_vec])
y_vec2 <- fitness_vec[order_vec]
tmp_df2 <- data.frame(y = y_vec2, x = x_vec2)
reg_res <- npregfast::frfast(y ~ x, data = tmp_df2)
y_vec_pred2 <- stats::predict(object = reg_res, newdata = data.frame(x = x_vec2))$Estimation[,"Pred"]

png(paste0("../../../../out/figures/Writeup6c/Writeup6c_peak_counter_DABTRAM_day0_nonzero-counts-sum_regression.png"),
    height = 1500, width = 3000, units = "px", res = 300)
n <- nrow(pred_diff_mat)
par(mar = c(4,4,4,0.5), mfrow = c(1,2))

tmp_lm <- stats::lm(y ~ x, data = tmp_df)
coef_val <- stats::coef(tmp_lm)["x"]
p_val <- -log10(summary(tmp_lm)$coefficients["x",4])
plot(x = x_vec, y = y_vec, pch = 16, col = rgb(0.75, 0.75, 0.75, 0.3),
     xlab = "Log(future lineage size at day10)", ylab = "Sum of non-peak counts", 
     main = paste0("No smoothing, LR coef: ", round(coef_val,2),
                   "\nCoef -Log10(pval): ", round(p_val,2)),
     cex.lab = 0.8)
lines(x = rep(log1p(20),2), y = c(-1e7,1e7), col = 2)
lines(x = x_vec, y = y_vec_pred, col = "white", lwd = 6)
lines(x = x_vec, y = y_vec_pred, col = "black", lwd = 4)

tmp_lm <- stats::lm(y ~ x, data = tmp_df2)
coef_val <- stats::coef(tmp_lm)["x"]
p_val <- -log10(summary(tmp_lm)$coefficients["x",4])
plot(x = x_vec2, y = y_vec2, pch = 16, col = rgb(0.75, 0.75, 0.75, 0.3),
     xlab = "Log(future lineage size at day10, smoothed by 10 NN)", ylab = "Sum of non-peak counts", 
     main = paste0("Smoothing the fitness, LR coef: ", round(coef_val,2),
                   "\nCoef -Log10(pval): ", round(p_val,2)),
     cex.lab = 0.8)
lines(x = x_vec2, y = y_vec_pred2, col = "white", lwd = 6)
lines(x = x_vec2, y = y_vec_pred2, col = "black", lwd = 4)

graphics.off()

#########################################

# let's try the same calculation above, but only for previously-discovered DABTRAM-related genes
source("../Writeup6b/gene_list.R")
gene_vec <- keygenes$DABTRAM
gene_vec <- intersect(gene_vec, rownames(countmat_nopeak_t))

cell_idx <- intersect(intersect(which(all_data$dataset == "day0"),
                                which(all_data$assigned_posterior > 0.5)),
                      which(!is.na(all_data$assigned_lineage)))
tmp <- countmat_nopeak_t[gene_vec,]
fitness_vec <- sapply(cell_idx, function(i){
  sum(.nonzero_col(tmp, i, bool_value = T))
})

lineage_vec <- all_data$assigned_lineage[cell_idx]
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
future_num_vec <- sapply(lineage_vec, function(lineage){
  tab_mat[lineage,"day10_DABTRAM"]
})

order_vec <- order(future_num_vec, decreasing = F)
x_vec <- log1p(future_num_vec[order_vec])
y_vec <- fitness_vec[order_vec]
tmp_df <- data.frame(y = y_vec, x = x_vec)
reg_res <- npregfast::frfast(y ~ x, data = tmp_df)
y_vec_pred <- stats::predict(object = reg_res, newdata = data.frame(x = x_vec))$Estimation[,"Pred"]

# compute nearest-neighbors, just among the valid day0 cells
set.seed(10)
mat <- all_data[["pca"]]@cell.embeddings[cell_idx,]
num_neigh <- 10
nn_mat <- RANN::nn2(mat, k = num_neigh+1)$nn.idx

# form smoothing nearest-neighbor matrix
n <- nrow(nn_mat)
i_vec <- rep(1:n, each = ncol(nn_mat))
j_vec <- unlist(lapply(1:n, function(i){nn_mat[i,]}))
avg_mat <- Matrix::sparseMatrix(i = i_vec, 
                                j = j_vec, 
                                x = rep(1/ncol(nn_mat), length(i_vec)), 
                                dims = c(n,n))
avg_mat <- Matrix::t(avg_mat)
future_num_vec_smoothed <- as.numeric(future_num_vec %*% avg_mat)

order_vec <- order(future_num_vec_smoothed, decreasing = F)
x_vec2 <- log1p(future_num_vec_smoothed[order_vec])
y_vec2 <- fitness_vec[order_vec]
tmp_df2 <- data.frame(y = y_vec2, x = x_vec2)
reg_res <- npregfast::frfast(y ~ x, data = tmp_df2)
y_vec_pred2 <- stats::predict(object = reg_res, newdata = data.frame(x = x_vec2))$Estimation[,"Pred"]

png(paste0("../../../../out/figures/Writeup6c/Writeup6c_peak_counter_DABTRAM_day0_nonzero-counts-sum_regression_keygenes.png"),
    height = 1500, width = 3000, units = "px", res = 300)
n <- nrow(pred_diff_mat)
par(mar = c(4,4,4,0.5), mfrow = c(1,2))

tmp_lm <- stats::lm(y ~ x, data = tmp_df)
coef_val <- stats::coef(tmp_lm)["x"]
p_val <- -log10(summary(tmp_lm)$coefficients["x",4])
plot(x = x_vec, y = y_vec, pch = 16, col = rgb(0.75, 0.75, 0.75, 0.3),
     xlab = "Log(future lineage size at day10)", ylab = "Sum of non-peak counts (only key genes)", 
     main = paste0("No smoothing, LR coef: ", round(coef_val,2),
                   "\nCoef -Log10(pval): ", round(p_val,2)),
     cex.lab = 0.8)
lines(x = rep(log1p(20),2), y = c(-1e7,1e7), col = 2)
lines(x = x_vec, y = y_vec_pred, col = "white", lwd = 6)
lines(x = x_vec, y = y_vec_pred, col = "black", lwd = 4)

tmp_lm <- stats::lm(y ~ x, data = tmp_df2)
coef_val <- stats::coef(tmp_lm)["x"]
p_val <- -log10(summary(tmp_lm)$coefficients["x",4])
plot(x = x_vec2, y = y_vec2, pch = 16, col = rgb(0.75, 0.75, 0.75, 0.3),
     xlab = "Log(future lineage size at day10, smoothed by 10 NN)", ylab = "Sum of non-peak counts (only key genes)", 
     main = paste0("Smoothing the fitness, LR coef: ", round(coef_val,2),
                   "\nCoef -Log10(pval): ", round(p_val,2)),
     cex.lab = 0.8)
lines(x = x_vec2, y = y_vec_pred2, col = "white", lwd = 6)
lines(x = x_vec2, y = y_vec_pred2, col = "black", lwd = 4)

graphics.off()