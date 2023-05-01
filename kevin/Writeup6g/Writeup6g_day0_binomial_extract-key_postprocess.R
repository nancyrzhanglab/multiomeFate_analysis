rm(list=ls())

load("../../../../out/kevin/Writeup6g/day0_binomial_extract-key_tmp.RData")

pvalue_list <- lapply(result_list, function(lis){
  len <- length(lis)
  pvalue_mat <- sapply(1:len, function(i){
    print(i)
    if(length(lis[[i]]) == 0) return(rep(NA,2))
    
    bin_win <- lis[[i]][1,]
    bin_die <- lis[[i]][2,]
    
    if(all(bin_win[c("Mid", "Close")] == 0) | all(bin_die[c("Mid", "Close")] == 0)){
      val1 <- NA
    } else {
      x_mat <- rbind(bin_win[c("Mid", "Close")], bin_die[c("Mid", "Close")])
      res1 <- stats::prop.test(x = x_mat, alternative = "greater")
      val1 <- res1$p.value
    }
    
    if(all(bin_win[c("Far", "Close")] == 0) | all(bin_die[c("Far", "Close")] == 0)){
      val2 <- NA
    } else {
      x_mat <- rbind(bin_win[c("Far", "Close")], bin_die[c("Far", "Close")])
      res2 <- stats::prop.test(x = x_mat, alternative = "greater")
      val2 <- res2$p.value
    }
   
    c(val1, val2)
  })
  colnames(pvalue_mat) <- names(lis)
  
  pvalue_mat
})

for(i in 1:3){
  print("====")
  print(i)
  zz <- round(-log10(pvalue_list[[i]]),2)
  idx <- which(zz > 2, arr.ind = T)
  # print(zz[,sort(unique(idx[,2]))])
  
  gene_tmp <- sort(unique(colnames(zz[,idx[,2]])))
  print(result_list[[i]][gene_tmp])
}

