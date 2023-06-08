# try data fission
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

topic_mat <- all_data[["fasttopic_DABTRAM"]]@cell.embeddings[names(cell_lineage),]
pseudotime_vec <- all_data$pseudotime[names(cell_lineage)]
atac_mat <- all_data2[["lsi"]]@cell.embeddings

log10pval_vec <- sapply(1:ncol(atac_mat), function(j){
  x <- atac_mat[names(cell_lineage),j]
  y <- as.factor(all_data2$tier_vec[names(cell_lineage)])
  res <- stats::oneway.test(x ~ y)
  -log10(res$p.value)
})
names(log10pval_vec) <- colnames(atac_mat)

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

##################

# try to fission the cell_features
# https://arxiv.org/abs/2112.11079
set.seed(10)
n <- nrow(cell_features)
z_mat <- matrix(0, nrow = n, ncol = ncol(cell_features))
for(j in 2:ncol(z_mat)){
  sd_val <- stats::sd(cell_features[,j])
  z_mat[,j] <- stats::rnorm(n = n, mean = 0, sd = sd_val)
}
training_mat <- cell_features + z_mat
testing_mat <- cell_features - z_mat
training_mat2 <- training_mat
testing_mat2 <- testing_mat
coefficient_initial2 <- coefficient_initial

coefficient_list <- vector("list", length = 1)
iteration <- 1
while(ncol(training_mat2) > 1){
  print(paste0("On iteration ", iteration))
  
  if(iteration > 1){
    rm_var <- coefficient_list[[iteration-1]]$variable_to_be_removed
    var_idx <- which(colnames(training_mat2) == rm_var)
    training_mat2 <- training_mat2[,-var_idx,drop = F]
    testing_mat2 <- testing_mat2[,-var_idx,drop = F]
    coefficient_initial <- coefficient_initial[-var_idx]
  }
  
  set.seed(10)
  lineage_res <- multiomeFate:::lineage_imputation(cell_features = training_mat2,
                                                   cell_lineage = cell_lineage,
                                                   coefficient_initial = coefficient_initial,
                                                   lineage_future_count = lineage_future_count,
                                                   random_initializations = 10,
                                                   verbose = 1)
  coefficient_vec <- lineage_res$fit$coefficient_vec
  variable_to_be_removed <- names(coefficient_vec)[which.min(abs(coefficient_vec))]
  
  loglik_val <- multiomeFate:::evaluate_loglikelihood(cell_features = testing_mat2,
                                                      cell_lineage = cell_lineage,
                                                      coefficient_vec = lineage_res$fit$coefficient_vec,
                                                      lineage_future_count = lineage_future_count)
  print(paste0("Testing log-likelihood: ", loglik_val))
  lis <- list(coefficient_vec = coefficient_vec,
              testing_obj = loglik_val,
              training_obj = lineage_res$fit$objective_val,
              variable_to_be_removed = variable_to_be_removed)
  
  coefficient_list[[iteration]] <- lis
  iteration <- iteration+1
}

save(coefficient_list, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup6i/Writeup6i_DABTRAM_lineage-imputation3_stepdown-variable.RData")

#########

len_vec <- sapply(coefficient_list, function(x){length(x$coefficient_vec)})
test_vec <- sapply(coefficient_list, function(x){x$testing_obj})
train_vec <- sapply(coefficient_list, function(x){x$training_obj})


png("../../../../out/figures/Writeup6i/Writeup6i_DABTRAM-day10_imputation3_stepdown-variable_loglikelihood.png", 
    width = 1500, height = 1500, units = "px", res = 300)
ylim <- range(c(test_vec, train_vec))

plot(x = len_vec, y = train_vec, col = "black", 
     lwd = 2, lty = 2,
     ylim = ylim,
     type = "o", pch = 16, cex = 2,
     xlab = "# variables", ylab = "-Loglikelihood", 
     main = "Stepdown variable selection: DABTRAM\nBlack:Train, Red:Test")
points(x = len_vec, y = test_vec, col = "red", pch = 16, cex = 2)
lines(x = len_vec, y = test_vec, col = "red", lwd = 2)
graphics.off()
