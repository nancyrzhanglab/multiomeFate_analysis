rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

treatment <- "DABTRAM"
load(paste0("../../../../out/kevin/Writeup6i/day0_cutmatrix_extract-", treatment, "_tmp.RData"))

rm_idx <- which(sapply(cutmat_list, function(x){
  all(is.null(x)) | all(is.na(x)) | length(x) == 0
}))
if(length(rm_idx) > 0) cutmat_list <- cutmat_list[-rm_idx]

peak_max <- sapply(cutmat_list, function(x){
  max(x$peak_prior)
})
length(which(peak_max <= 0.1))
length(which(peak_max <= 0.1))/length(peak_max)
names(peak_max)[which(peak_max <= 0.1)]

half_prior_10 <- sapply(cutmat_list, function(x){
  y <- x$peak_prior
  sum(y[y >= 0.1]) >= 0.5
})
table(half_prior_10)

half_prior_5 <- sapply(cutmat_list, function(x){
  y <- x$peak_prior
  sum(y[y >= 0.05]) >= 0.5
})
table(half_prior_5)

###################

# let's focus on half_prior_10 -- the easier genes
cutmat_list2 <- cutmat_list[which(half_prior_10)]

# clean up the peaks
len <- length(cutmat_list2)
for(i in 1:len){
  peak_mat <- cutmat_list2[[i]]$peak_mat
  peak_prior <- cutmat_list2[[i]]$peak_prior
  
  rm_idx <- which(peak_prior <= 0.1)
  if(length(rm_idx) > 0){
    peak_mat <- peak_mat[-rm_idx,,drop = F]
    peak_prior <- peak_prior[-rm_idx]
    peak_prior <- peak_prior/sum(peak_prior)
  }
  
  peak_midpoint <- sapply(1:nrow(peak_mat), function(j){
    mean(peak_mat[j,])
  })
  
  cutmat_list2[[i]]$peak_mat <- peak_mat
  cutmat_list2[[i]]$peak_prior <- peak_prior
  cutmat_list2[[i]]$peak_midpoint <- peak_midpoint
}

# for each gene, for each fragment, compute the absolute distance to the nearest peak

len <- length(cutmat_list2)
dist_list <- vector("list", length = len)
names(dist_list) <- names(cutmat_list2)

.extract_distances <- function(cutmat,
                               norm_val,
                               peak_midpoint){
  frag_locations <- multiomeFate:::.extract_fragment_from_cutmat(cutmat)
  if(length(frag_locations) == 0){
    return(NA)
  }
  dist_vec <- sapply(frag_locations, function(x){
    min(abs(x - peak_midpoint))
  })
  dist_vec <- dist_vec/norm_val
  dist_vec
}

for(i in 1:len){
  if(i %% floor(len/10) == 0) cat('*')
  
  peak_mat <- cutmat_list2[[i]]$peak_mat
  peak_midpoint <- cutmat_list2[[i]]$peak_midpoint
  
  dist_winning <- .extract_distances(
    cutmat = cutmat_list2[[i]]$cutmat_winning,
    norm_val = stats::median(apply(peak_mat, 1, diff)),
    peak_midpoint = peak_midpoint
  )
  dist_dying <- .extract_distances(
    cutmat = cutmat_list2[[i]]$cutmat_dying,
    norm_val = stats::median(apply(peak_mat, 1, diff)),
    peak_midpoint = peak_midpoint
  )
  
  dist_list[[i]] <- list(dist_winning = dist_winning,
                         dist_dying = dist_dying)
}

rm_idx <- which(sapply(dist_list, function(lis){
  all(is.na(lis$dist_winning)) | all(is.na(lis$dist_dying))
}))
if(length(rm_idx) > 0) dist_list <- dist_list[-rm_idx]

gene <- "CD44"
quantile(dist_list[[gene]]$dist_winning)
quantile(dist_list[[gene]]$dist_dying)

###############

.ecdf <- function(values,
                  tol = 1e-6){
  stopifnot(all(values >= 0))
  
  weights <- rep(1, length(values))
  mat <- cbind(c(0, values), c(0, weights))
  colnames(mat) <- c("x", "w")
  
  # compute cdf
  mat <- mat[order(mat[,"x"]),]
  cdf_vec <- cumsum(mat[,"w"])
  cdf_vec <- cdf_vec/max(cdf_vec)
  
  # clean up duplicated values
  tmp <- .remove_duplicates(associated_vec = cdf_vec,
                            target_vec = mat[,"x"])
  
  list(cdf = tmp$associated_vec,
       x = tmp$target_vec)
}

# assumes target_vec is sorted already, and associated_vec is either NULL or of equal length to target_vec
.remove_duplicates <- function(associated_vec,
                               target_vec,
                               tol = 1e-6){
  stopifnot(length(associated_vec) == length(target_vec))
  
  diff_vec <- diff(target_vec)
  idx <- which(diff_vec <= tol)
  if(length(idx) > 0){
    target_vec <- target_vec[-idx]
    associated_vec <- associated_vec[-idx]
  }
  
  list(associated_vec = associated_vec,
       target_vec = target_vec)
}

# run LCM from https://search.r-project.org/CRAN/refmans/fdrtool/html/gcmlcm.html
# outputs are: x.knots, y.knots, slope.knots
.compute_decreasing_density <- function(cdf, x){
  res <- fdrtool::gcmlcm(x = x,
                         y = cdf,
                         type = "lcm")
  
  list(x = res$x,
       pdf = c(res$slope.knots,0))
}

###############

len <- length(dist_list)
lrt_vec <- sapply(1:len, function(i){
  if(i %% floor(len/10) == 0) cat('*')
  
  tol <- 1e-6
  
  val_win <- dist_list[[i]]$dist_winning
  val_win <- val_win[val_win >= tol]
  ecdf_res <- .ecdf(val_win)
  decr_density_win <- .compute_decreasing_density(cdf = ecdf_res$cdf,
                                                  x = ecdf_res$x)
  # density is left-continuous, meaning you want to index right above x-tol
  loglik_win <- sum(sapply(val_win, function(x){
    idx <- max(which(decr_density_win$x < x))
    log(decr_density_win$pdf[idx])
  }))
  
  val_die <- dist_list[[i]]$dist_dying
  val_die <- val_die[val_die >= tol]
  ecdf_res <- .ecdf(val_die)
  decr_density_die <- .compute_decreasing_density(cdf = ecdf_res$cdf,
                                                  x = ecdf_res$x)
  loglik_die <- sum(sapply(val_die, function(x){
    idx <- max(which(decr_density_die$x < x))
    log(decr_density_die$pdf[idx])
  }))
  
  val_all <- c(val_win, val_die)
  ecdf_res <- .ecdf(val_all)
  decr_density_all <- .compute_decreasing_density(cdf = ecdf_res$cdf,
                                                  x = ecdf_res$x)
  loglik_all <- sum(sapply(val_all, function(x){
    idx <- max(which(decr_density_all$x < x))
    log(decr_density_all$pdf[idx])
  }))
  
  -2*(loglik_all - (loglik_win+loglik_die))
})
names(lrt_vec) <- names(dist_list) 

quantile(lrt_vec)

df <- mean(lrt_vec[lrt_vec <= stats::quantile(lrt_vec, probs = 0.95)])

x_vec <- seq(0, max(lrt_vec), length.out = 1000)
y_vec <- stats::dchisq(x_vec, df = df)
png("../../../../out/figures/Writeup6i/Writeup6i_DABTRAM-day0_grenander_lrt-teststat.png", 
    width = 1500, height = 1500, units = "px", res = 300)
res <- hist(lrt_vec, breaks = 50, col = "gray", 
            xlab = "LRT test stat", 
            main = "DABTRAM Day0 Grenander LRT", plot = T)
y_vec <- y_vec*(max(res$counts)/max(y_vec))
lines(x = x_vec, y = y_vec, lwd = 2, lty = 2, col = "red")
graphics.off()

pval_vec <- stats::pchisq(lrt_vec, df = df, lower.tail = F, log.p = F)
names(pval_vec) <- names(dist_list)
pval_adj_vec <- stats::p.adjust(pval_vec)
quantile(pval_adj_vec)
length(which(pval_adj_vec <= 0.05))
pval_adj_vec[which(pval_adj_vec <= 0.05)]



