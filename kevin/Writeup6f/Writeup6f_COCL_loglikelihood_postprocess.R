rm(list=ls())
load("../../../../out/kevin/Writeup6f/COCL_loglikelihood_tmp.RData")

idx <- which(sapply(result_list, length)>0)
lrt_vec <- sapply(idx, function(i){
  ll_win <- result_list[[i]]$res_win$loglikelihood_val
  ll_die <- result_list[[i]]$res_die$loglikelihood_val
  ll_both <- result_list[[i]]$res_both$loglikelihood_val
  
  lrt <- -2 * (ll_both - (ll_win + ll_die))
})