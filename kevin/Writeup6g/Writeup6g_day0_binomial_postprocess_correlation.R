rm(list=ls())

load("../../../../out/kevin/Writeup6g/day0_binomial_extract-all_CIS_tmp.RData")
result_list_bg <- result_list
load("../../../../out/kevin/Writeup6g/day0_binomial_extract-all_COCL2_tmp.RData")
result_list_bg[["COCL2"]] <- result_list[["COCL2"]]
load("../../../../out/kevin/Writeup6g/day0_binomial_extract-key_tmp.RData")

#########################

pvalue_func <- function(mat){
  bin_win <- mat[1,]
  bin_die <- mat[2,]
  
  if(all(bin_win[c("Mid", "Close")] == 0) | all(bin_die[c("Mid", "Close")] == 0)){
    val1 <- NA
  } else {
    x_mat <- rbind(bin_win[c("Mid", "Close")], bin_die[c("Mid", "Close")])
    res1 <- stats::prop.test(x = x_mat, alternative = "greater")
    val1 <- -log10(res1$p.value)
  }
  
  if(all(bin_win[c("Far", "Close")] == 0) | all(bin_die[c("Far", "Close")] == 0)){
    val2 <- NA
  } else {
    x_mat <- rbind(bin_win[c("Far", "Close")], bin_die[c("Far", "Close")])
    res2 <- stats::prop.test(x = x_mat, alternative = "greater")
    val2 <- -log10(res2$p.value)
  }
  
  c(val1, val2)
}

########################

len <- length(result_list[["CIS"]])
log10pval_mat_cis <- sapply(1:len, function(i){
  if(i %% floor(len/10) == 0) cat('*')
  if(length(result_list[["CIS"]][[i]]) == 0) return(rep(NA,2))
  pvalue_func(result_list[["CIS"]][[i]])
})
colnames(log10pval_mat_cis) <- names(result_list[["CIS"]])
rownames(log10pval_mat_cis) <- c("Mid", "Far")

log10pval_mat_cocl2 <- sapply(1:len, function(i){
  if(i %% floor(len/10) == 0) cat('*')
  if(length(result_list[["COCL2"]][[i]]) == 0) return(rep(NA,2))
  pvalue_func(result_list[["COCL2"]][[i]])
})
colnames(log10pval_mat_cocl2) <- names(result_list[["COCL2"]])
rownames(log10pval_mat_cocl2) <- c("Mid", "Far")

#####

result_list_bg[["CIS"]] <- result_list_bg[["CIS"]][which(sapply(result_list_bg[["CIS"]], length) > 0)]

len <- length(result_list_bg[["CIS"]])
log10pval_mat_cis_bg <- sapply(1:len, function(i){
  if(i %% floor(len/10) == 0) cat('*')
  if(length(result_list_bg[["CIS"]][[i]]) == 0) return(rep(NA,2))
  pvalue_func(result_list_bg[["CIS"]][[i]])
})
colnames(log10pval_mat_cis_bg) <- names(result_list_bg[["CIS"]])
rownames(log10pval_mat_cis_bg) <- c("Mid", "Far")


len <- length(result_list_bg[["COCL2"]])
log10pval_mat_cocl2_bg <- sapply(1:len, function(i){
  if(i %% floor(len/10) == 0) cat('*')
  if(length(result_list_bg[["COCL2"]][[i]]) == 0) return(rep(NA,2))
  pvalue_func(result_list_bg[["COCL2"]][[i]])
})
colnames(log10pval_mat_cocl2_bg) <- names(result_list_bg[["COCL2"]])
rownames(log10pval_mat_cocl2_bg) <- c("Mid", "Far")


################################

png(paste0("../../../../out/figures/Writeup6g/Writeup6g_day0_binomial-pvalue_correlation_CIS-COCL2.png"),
    height = 1200, width = 2500, res = 300, units = "px")
par(mfrow = c(1,2))
vec1 <- log10pval_mat_cis_bg[1,]; vec2 <- log10pval_mat_cocl2_bg[1,]
vec1 <- vec1[!is.na(vec1)]; vec2 <- vec2[!is.na(vec2)]
gene_vec <- sort(intersect(names(vec1), names(vec2)))
vec1 <- vec1[gene_vec]; vec2 <- vec2[gene_vec]

plot(vec1, vec2, xlab = "-Log10 p-value (CIS)", ylab = "-Log10 p-value (COCL2)",
     main = "Mid vs. close (CIS-COCL2)",
     pch = 16, col = rgb(0.5,0.5,0.5,0.5))
points(log10pval_mat_cis[1,],
       log10pval_mat_cocl2[1,],
       pch = 16, col = "firebrick", cex = 1.5)

##

vec1 <- log10pval_mat_cis_bg[2,]; vec2 <- log10pval_mat_cocl2_bg[2,]
vec1 <- vec1[!is.na(vec1)]; vec2 <- vec2[!is.na(vec2)]
gene_vec <- sort(intersect(names(vec1), names(vec2)))
vec1 <- vec1[gene_vec]; vec2 <- vec2[gene_vec]

plot(vec1, vec2, xlab = "-Log10 p-value (CIS)", ylab = "-Log10 p-value (COCL2)",
     main = "Mid vs. far (CIS-COCL2)",
     pch = 16, col = rgb(0.5,0.5,0.5,0.5))
points(log10pval_mat_cis[2,],
       log10pval_mat_cocl2[2,],
       pch = 16, col = "firebrick", cex = 1.5)
graphics.off()
