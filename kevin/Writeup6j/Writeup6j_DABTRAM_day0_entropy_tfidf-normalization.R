rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

treatment <- "DABTRAM"
load(paste0("../../../../out/kevin/Writeup6i/day0_cutmatrix_extract-", treatment, "_tmp.RData"))

source("entropy_normalization.R")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

########

# idx <- which(names(cutmat_list) == "CD44")

len <- length(cutmat_list)
preprocessed_list <- vector("list", len)
names(preprocessed_list) <- names(cutmat_list)
for(idx in 1:len){
  if(idx %% floor(len/10) == 0) {
    cat('*')
    save(date_of_run, session_info,
         preprocessed_list,
         file = "../../../../out/kevin/Writeup6j/Writeup6j_DABTRAM_day0_entropy_tfidf-normalization.RData")
  }
  
  cutmat <- rbind(cutmat_list[[idx]]$cutmat_dying,
                  cutmat_list[[idx]]$cutmat_winning)
  peak_mat <- cutmat_list[[idx]]$peak_mat
  peak_prior <- cutmat_list[[idx]]$peak_prior
  thres1 <- 0.1; thres2 <- 0.05
  if(any(peak_prior >= thres1)){
    peak_mat <- peak_mat[peak_prior >= thres1,,drop = F]
    peak_prior <- peak_prior[peak_prior >= thres1]
    peak_prior <- peak_prior/sum(peak_prior)
  } else if(any(peak_prior >= thres2)){
    peak_mat <- peak_mat[peak_prior >= thres2,,drop = F]
    peak_prior <- peak_prior[peak_prior >= thres2]
    peak_prior <- peak_prior/sum(peak_prior)
  }
  
  peak_width <- stats::median(apply(peak_mat, 1, diff))
  data_mat <- .normalize_fragments(cutmat = cutmat,
                                   peak_mat = peak_mat,
                                   normalize_by_unique_cells = T)
  
  preprocessed_list[[idx]] <- list(
    idx_dying = which(data_mat$cell_name %in% rownames(cutmat_list[[idx]]$cutmat_dying)),
    idx_winning = which(data_mat$cell_name %in% rownames(cutmat_list[[idx]]$cutmat_winning)),
    data_mat = data_mat,
    peak_width = peak_width
  )
}

save(date_of_run, session_info,
     preprocessed_list,
     file = "../../../../out/kevin/Writeup6j/Writeup6j_DABTRAM_day0_entropy_tfidf-normalization.RData")

# break_vec <- sort(c(0, peak_width, 2*peak_width, 5*peak_width, max(data_mat$dist_to_peak)))
# table(cut(data_mat$dist_to_peak, breaks = break_vec),
#       cut(data_mat$tfidf, breaks = quantile(data_mat$tfidf)))
# 
# table(cut(data_mat$nearby_frags, breaks = quantile(data_mat$nearby_frags)),
#       cut(data_mat$tfidf, breaks = quantile(data_mat$tfidf)))


lrt_vec <- sapply(1:len, function(idx){
  if(idx %% floor(len/10) == 0) cat('*')
  
  data_mat <- preprocessed_list[[idx]]$data_mat
  idx_dying <- preprocessed_list[[idx]]$idx_dying
  idx_winning <- preprocessed_list[[idx]]$idx_winning
  peak_width <- preprocessed_list[[idx]]$peak_width
  
  grenander_win <- multiomeFate:::estimate_grenander(
    values = data_mat$dist_to_peak[idx_winning],
    weights = data_mat$tfidf[idx_winning],
    # weights = rep(1, length(idx_winning)),
    scaling_factor = peak_width
  )
  loglik_win <- sum(sapply(data_mat$dist_to_peak[idx_winning], function(x){
    multiomeFate:::evaluate_grenander(
      obj = grenander_win,
      x = x,
      bool_log = T
    )
  }))
  
  grenander_die <- multiomeFate:::estimate_grenander(
    values = data_mat$dist_to_peak[idx_dying],
    weights = data_mat$tfidf[idx_dying],
    # weights = rep(1, length(idx_dying)),
    scaling_factor = peak_width
  )
  loglik_die <- sum(sapply(data_mat$dist_to_peak[idx_dying], function(x){
    multiomeFate:::evaluate_grenander(
      obj = grenander_die,
      x = x,
      bool_log = T
    )
  }))
  
  grenander_all <- multiomeFate:::estimate_grenander(
    values = data_mat$dist_to_peak,
    weights = data_mat$tfidf,
    # weights = rep(1, length(data_mat$dist_to_peak)),
    scaling_factor = peak_width
  )
  loglik_all <- sum(sapply(data_mat$dist_to_peak, function(x){
    multiomeFate:::evaluate_grenander(
      obj = grenander_all,
      x = x,
      bool_log = T
    )
  }))
  
  -2*(loglik_all - (loglik_win+loglik_die))
})
names(lrt_vec) <- names(cutmat_list) 

save(date_of_run, session_info,
     preprocessed_list, lrt_vec,
     file = "../../../../out/kevin/Writeup6j/Writeup6j_DABTRAM_day0_entropy_tfidf-normalization.RData")

