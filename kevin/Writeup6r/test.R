rm(list=ls())
library(Seurat)
library(CCA)

treatment_vec <- c("CIS", "COCL2", "DABTRAM")
day_early_vec <- c("day0", "day10")
day_later_vec <- c("day10", "week5")

load("../../../../out/kevin/Writeup6p/Writeup6p_all-data_lightweight_noATAC.RData")

treatment <- "CIS"
keep_vec <- rep(FALSE, ncol(all_data))
idx <- which(all_data$dataset %in% c("day0", paste0("day10_", treatment), paste0("week5_", treatment)))
keep_vec[idx] <- TRUE
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == TRUE)

ii <- 1
day_early <- day_early_vec[ii]
day_later <- day_later_vec[ii]
load(paste0("../../../../out/kevin/Writeup6r/Writeup6r_", treatment, "_", day_early, "_lineage-imputation.RData"))

if(day_early == "day0") {
  day_early_full <- day_early
} else {
  day_early_full <- paste0(day_early, "_", treatment)
}
day_later_full <- paste0(day_later, "_", treatment)

rna_mat <- all_data[[paste0("fasttopic_", treatment)]]@cell.embeddings[rownames(cell_features),]
atac_mat <- all_data[[paste0("peakVI_", treatment)]]@cell.embeddings[rownames(cell_features),]

rm(list = "all_data")
gc()

###################################

train_mat <- sapply(cv_fit_list, function(x){
  x$train_loglik
})

test_mat <- sapply(cv_fit_list, function(x){
  x$test_loglik
})

train_quantile <- apply(train_mat, 1, function(vec){stats::quantile(vec, probs = c(0.1, 0.5, 0.9))})
test_quantile <- apply(test_mat, 1, function(vec){stats::quantile(vec, probs = c(0.1, 0.5, 0.9))})

lambda_sequence <- cv_fit_list[[1]]$train_fit$lambda_sequence
lambda <- lambda_sequence[which.min(test_quantile[2,])]

final_fit <- multiomeFate:::lineage_imputation(
  cell_features = cell_features,
  cell_lineage = cell_lineage,
  coefficient_initial_list = cv_fit_list[[1]]$train_fit$fit_list[[which.min(test_quantile[2,])]]$coefficient_vec,
  lambda = lambda,
  lineage_future_count = lineage_future_count,
  verbose = 0
)

cell_features <- cbind(1, cell_features)
colnames(cell_features)[1] <- "Intercept"
stopifnot(all(colnames(cell_features) == names(final_fit$fit$coefficient_vec)))
cell_imputed_score <- as.numeric(cell_features %*% final_fit$fit$coefficient_vec)
names(cell_imputed_score) <- rownames(cell_features)
cell_imputed_count <- exp(cell_imputed_score)
uniq_lineage <- sort(unique(cell_lineage))
lineage_imputed_count <- sapply(uniq_lineage, function(lineage){
  sum(cell_imputed_count[which(cell_lineage == lineage)])
})
cell_imputed_score2 <- log10(exp(cell_imputed_score)) # this one is on the log10 scale

###############################

df <- as.data.frame(cbind(cell_imputed_score2, rna_mat, atac_mat))
colnames(df)[1] <- "y"
tmp_lm <- lm(y ~ ., data = df)
quantile(tmp_lm$residuals)

rna_mat2 <- scale(rna_mat)
atac_mat2 <- scale(atac_mat)

n <- nrow(cell_features)
d <- min(ncol(rna_mat2), ncol(atac_mat2))

rna_svd <- svd(rna_mat2)
atac_svd <- svd(atac_mat2)

tol <- 1e-6
rna_idx <- which(rna_svd$d >= tol)
atac_idx <- which(atac_svd$d >= tol)

rna_basis <- rna_svd$u[,rna_idx]
atac_basis <- atac_svd$u[,atac_idx]

df <- as.data.frame(cbind(cell_imputed_score2, rna_basis, atac_basis))
colnames(df)[1] <- "y"
tmp_lm <- lm(y ~ ., data = df)
quantile(tmp_lm$residuals)

tmp_cc <- CCA::cc(rna_basis, atac_basis)
x_coef <- tmp_cc$xcoef
y_coef <- tmp_cc$ycoef
rna_z <- rna_basis %*% x_coef
atac_z <- atac_basis %*% y_coef
stats::cor(rna_z[,1], atac_z[,1])

cca_threshold <- 0.75
common_idx <- which(tmp_cc$cor >= cca_threshold)
common_mat <- cbind(rna_z[,common_idx,drop = F], atac_z[,common_idx,drop = F])
common_svd <- svd(common_mat)
common_basis <- common_svd$u
common_proj <- tcrossprod(common_basis)

common_and_atac_basis <- cbind(1, common_basis, atac_basis)
tmp <- svd(common_and_atac_basis)
common_and_atac_proj <- tcrossprod(tmp$u[,which(tmp$d >= tol),drop=F])
rna_uniq_basis <- rna_basis - common_and_atac_proj %*% rna_basis
tmp <- svd(rna_uniq_basis)
rna_uniq_basis <- tmp$u[,which(tmp$d >= tol),drop=F]

common_and_rna_basis <- cbind(1, common_basis, rna_basis)
tmp <- svd(common_and_rna_basis)
common_and_rna_proj <- tcrossprod(tmp$u[,which(tmp$d >= tol),drop=F])
atac_uniq_basis <- atac_basis - common_and_rna_proj %*% atac_basis
tmp <- svd(atac_uniq_basis)
atac_uniq_basis <- tmp$u[,which(tmp$d >= tol),drop=F]

round(cor(rna_uniq_basis, atac_uniq_basis),2) # ... i think these are the CCA values...
round(cor(rna_uniq_basis, atac_basis),2)



