rm(list=ls())

result_list_total <- vector("list", 0)
bulk_total <- 5
for(bulk_number in 1:bulk_total){
  load(paste0("../../../../out/kevin/Writeup6f/COCL_loglikelihood_bulk", bulk_number, "of", bulk_total, "_tmp.RData"))
  result_list_total[[bulk_number]] <- result_list
}
result_list_total <- do.call(c, result_list_total)
lrt_vec_all <- rep(NA, length(result_list_total))
names(lrt_vec_all) <- names(result_list_total)

result_list_total <- result_list_total[which(sapply(result_list_total, length) > 0)]
lrt_vec_bulk <- sapply(1:length(result_list_total), function(i){
  ll_win <- result_list_total[[i]]$res_win$loglikelihood_val
  ll_die <- result_list_total[[i]]$res_die$loglikelihood_val
  ll_both <- result_list_total[[i]]$res_both$loglikelihood_val
  
  lrt <- -2 * (ll_both - (ll_win + ll_die))
})
quantile(lrt_vec_bulk, probs = seq(0,1,length.out=11))
lrt_vec_all[names(result_list_total)] <- lrt_vec_bulk

############################

load("../../../../out/kevin/Writeup6f/COCL_loglikelihood_tmp.RData")
result_list <- result_list[which(sapply(result_list, length) > 0)]
length(result_list)
lrt_vec <- sapply(1:length(result_list), function(i){
  ll_win <- result_list[[i]]$res_win$loglikelihood_val
  ll_die <- result_list[[i]]$res_die$loglikelihood_val
  ll_both <- result_list[[i]]$res_both$loglikelihood_val
  
  lrt <- -2 * (ll_both - (ll_win + ll_die))
})
quantile(lrt_vec, probs = seq(0,1,length.out=11))
lrt_vec_all[names(result_list)] <- lrt_vec

length(which(!is.na(lrt_vec_all)))

save(lrt_vec_all, result_list,
     file = "../../../../out/kevin/Writeup6f/tmp.RData")

zz <- cbind(quantile(lrt_vec, probs = seq(0,1,length.out=11)), 
            quantile(lrt_vec_bulk, probs = seq(0,1,length.out=11)))
round(zz)

################################

rm(list=ls())
load("../../out/Writeup6f/tmp.RData")

lrt_vec_all <- lrt_vec_all[!is.na(lrt_vec_all)]
lrt_vec_all <- lrt_vec_all[sort(names(lrt_vec_all))]
hist(lrt_vec_all, breaks = 50)

fn <- function(param_vec,
               N,
               N0,
               z0_vec,
               ub){
  df0 <- param_vec[1]
  theta <-  param_vec[2]
  if(theta > 0.99 | theta < 0.01) return(Inf)
  denom <- stats::pchisq(ub, df = df0)
  
  # compute each sample's likelihood
  sample_llvec <- sapply(z0_vec, function(z){
    stats::dchisq(z, df = df0, log = T)
  })
  
  # compute full log-likelihood, Equation 4.12
  loglik <- N0*log(theta) + (N-N0)*log(1-theta) + sum(sample_llvec) - N0*log(denom)
  -loglik
}

z_vec <- lrt_vec_all
N <- length(z_vec)
ub <- stats::quantile(z_vec, probs = 0.99)
idx <- which(z_vec <= ub)
N0 <- length(idx)
z0_vec <- z_vec[idx]

init_df <- mean(z0_vec)
init_theta <- 0.99* stats::pchisq(ub, df = init_df)

optim_res <- stats::optim(par = c(init_df, init_theta),
                          fn = fn,
                          method = "Nelder-Mead",
                          N = N,
                          N0 = N0,
                          z0_vec = z0_vec,
                          ub = ub)
optim_res$par
df_est <- optim_res$par[1]

hist_res <- hist(z_vec, breaks = 50, plot = F)
hist_res$counts <- hist_res$density

x_vec <- hist_res$breaks
midpoints <- sapply(2:length(x_vec), function(i){
  mean(x_vec[c(i,i-1)])
})
prob_vec <- sapply(2:length(x_vec), function(i){
  stats::pchisq(x_vec[i], df = df_est) - stats::pchisq(x_vec[i-1], df = df_est) 
})
plot(hist_res, col = "gray")
points(midpoints, prob_vec, pch = 16, col = 2)

####################

i <- which(names(result_list) == "NDRG1")
par(mfrow = c(1,2))
plot(result_list[[i]]$res_win$grenander_obj$x,
     result_list[[i]]$res_win$grenander_obj$pdf, main = "Win")
plot(result_list[[i]]$res_die$grenander_obj$x,
     result_list[[i]]$res_die$grenander_obj$pdf, main = "Die")

#####################

# what if we compute the center of mass
obj <- result_list[[i]]$res_win$grenander_obj
bin_cutoff <- seq(min(obj$x), max(obj$x), length.out = 1000)[-c(1,1000)]
obj <- .add_cutoffs_to_grenander(obj = obj,
                                 bin_cutoff = bin_cutoff)
area_vec <- cumsum(diff(obj$x)*obj$pdf[-length(obj$pdf)])
idx <- which.min(abs(area_vec - 0.5))
if(area_vec[idx] >= 0.5){
  # take the left
  midpoint <- obj$x[idx]
} else {
  # take the right
  midpoint <- obj$x[idx+1]
}
midpoint*obj$scaling_factor/933

obj <- result_list[[i]]$res_die$grenander_obj
bin_cutoff <- seq(min(obj$x), max(obj$x), length.out = 1000)[-c(1,1000)]
obj <- .add_cutoffs_to_grenander(obj = obj,
                                 bin_cutoff = bin_cutoff)
area_vec <- cumsum(diff(obj$x)*obj$pdf[-length(obj$pdf)])
idx <- which.min(abs(area_vec - 0.5))
if(area_vec[idx] >= 0.5){
  # take the left
  midpoint <- obj$x[idx]
} else {
  # take the right
  midpoint <- obj$x[idx+1]
}
midpoint*obj$scaling_factor/933

plot(result_list[[i]]$res_die$grenander_obj$x,
     log(result_list[[i]]$res_die$grenander_obj$pdf))
points(result_list[[i]]$res_win$grenander_obj$x,
     log(result_list[[i]]$res_win$grenander_obj$pdf), col = 2, cex = 0.5, pch = 16)

obj_win <- result_list[[i]]$res_win$grenander_obj
obj_die <- result_list[[i]]$res_die$grenander_obj
x_lim <- c(max(min(obj_win$x), min(obj_die$x)),
           min(max(obj_win$x), max(obj_die$x)))
x_vec <- seq(x_lim[1], x_lim[2], length.out = 1000)[-1]*obj_win$scaling_factor
pdf_win <- sapply(x_vec, function(x){evaluate_grenander(obj = obj_win, x = x)})
pdf_die <- sapply(x_vec, function(x){evaluate_grenander(obj = obj_die, x = x)})

plot(pdf_win, pdf_die, asp = T); lines(c(-1e6,1e6), c(-1e6,1e6), col = 2, lwd = 2, lty = 2)

# [[IDEAS ABOUT APPROXIMATING IT BY AN EXPONENTIAL IS BORKED]]
# x_vec <- seq(min(result_list[[i]]$res_win$grenander_obj$x), 
#              max(result_list[[i]]$res_win$grenander_obj$x),
#              length.out = 2000)[-1] * result_list[[i]]$res_win$grenander_obj$scaling_factor
# pdf_vec <- sapply(x_vec, function(x){
#   evaluate_grenander(obj = result_list[[i]]$res_win$grenander_obj,
#                      x = x)})
# log_vec <- log(pdf_vec)
# x_vec <- x_vec / 1000
# 
# fn <- function(lambda, log_vec, x_vec){
#   theoretical_log_vec <- log(lambda) - lambda*x_vec
#   sum((theoretical_log_vec - log_vec)^2)/length(x_vec)
# }
# 
# opt_res <- stats::optimise(f = fn, interval = c(1e-6,100),
#                            log_vec = log_vec, x_vec = x_vec,
#                            maximum = F)
# 
# plot(x_vec, log(pdf_vec))
# lambda_est <- opt_res$minimum
# theoretical_pdf <- sapply(x_vec, function(x){
#   lambda_est*exp(-lambda_est*x)
# })
# points(x_vec, log(theoretical_pdf), col = 2, pch = 16)
# 
# # theoretical_pdf <- sapply(2:length(x_vec), function(i){
# #   stats::pexp(x_vec[i], rate = lambda_est) - stats::pexp(x_vec[i-1], rate = lambda_est) 
# # })
# # x_midpoint <- sapply(2:length(x_vec), function(i){
# #   mean(x_vec[c(i,i-1)])
# # })
# # points(x_midpoint, theoretical_pdf, col = 2, pch = 16, cex = 0.5)
# 
