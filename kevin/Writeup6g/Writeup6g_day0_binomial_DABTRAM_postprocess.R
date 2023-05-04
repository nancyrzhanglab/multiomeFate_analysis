rm(list=ls())

treatment <- "DABTRAM"
load(paste0("../../../../out/kevin/Writeup6g/day0_binomial_extract-all_", treatment, "_tmp.RData"))
result_list_bg <- result_list[[treatment]]
load("../../../../out/kevin/Writeup6g/day0_binomial_extract-key_tmp.RData")

pvalue_func <- function(mat){
  bin_win <- mat[1,]
  bin_die <- mat[2,]
  
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
}

len <- length(result_list[["DABTRAM"]])
pvalue_mat <- sapply(1:len, function(i){
  if(i %% floor(len/10) == 0) cat('*')
  if(length(result_list[["DABTRAM"]][[i]]) == 0) return(rep(NA,2))
  pvalue_func(result_list[["DABTRAM"]][[i]])
})
colnames(pvalue_mat) <- names(result_list[["DABTRAM"]])
rownames(pvalue_mat) <- c("Mid", "Far")

###########
result_list_bg <- result_list_bg[which(sapply(result_list_bg, length) > 0)]

len <- length(result_list_bg)
pvalue_mat_bg <- sapply(1:len, function(i){
  if(i %% floor(len/10) == 0) cat('*')
  if(length(result_list_bg[[i]]) == 0) return(rep(NA,2))
  pvalue_func(result_list_bg[[i]])
})
colnames(pvalue_mat_bg) <- names(result_list_bg)
rownames(pvalue_mat_bg) <- c("Mid", "Far")

round(quantile(pvalue_mat_bg[1,],na.rm=T,probs=seq(0,1,length.out=11)),2)
round(quantile(pvalue_mat_bg[2,],na.rm=T,probs=seq(0,1,length.out=11)),2)

png(paste0("../../../../out/figures/Writeup6g/Writeup6g_day0_binomial-pvalue_", treatment, ".png"),
    height = 1200, width = 2500, res = 300, units = "px")
par(mfrow = c(1,2))
hist(pvalue_mat_bg[1,], breaks = 50, xlab = "P-value", main = "Mid vs. Close (DABTRAM)")
rug(pvalue_mat[1,], col = "firebrick", lwd = 2)

hist(pvalue_mat_bg[2,], breaks = 50, xlab = "P-value", main = "Far vs. Close (DABTRAM)")
rug(pvalue_mat[2,], col = "firebrick", lwd = 2)
graphics.off()

##############################

teststat_func <- function(mat){
  bin_win <- mat[1,]
  bin_die <- mat[2,]
  
  if(all(bin_win[c("Mid", "Close")] == 0) | all(bin_die[c("Mid", "Close")] == 0)){
    val1 <- NA
  } else {
    x_mat <- rbind(bin_win[c("Mid", "Close")], bin_die[c("Mid", "Close")])
    res1 <- stats::prop.test(x = x_mat, alternative = "greater")
    val1 <- res1$statistic
  }
  
  if(all(bin_win[c("Far", "Close")] == 0) | all(bin_die[c("Far", "Close")] == 0)){
    val2 <- NA
  } else {
    x_mat <- rbind(bin_win[c("Far", "Close")], bin_die[c("Far", "Close")])
    res2 <- stats::prop.test(x = x_mat, alternative = "greater")
    val2 <- res2$statistic
  }
  
  c(val1, val2)
}

len <- length(result_list[["DABTRAM"]])
teststat_mat <- sapply(1:len, function(i){
  if(i %% floor(len/10) == 0) cat('*')
  if(length(result_list[["DABTRAM"]][[i]]) == 0) return(rep(NA,2))
  teststat_func(result_list[["DABTRAM"]][[i]])
})
colnames(teststat_mat) <- names(result_list[["DABTRAM"]])
rownames(teststat_mat) <- c("Mid", "Far")

len <- length(result_list_bg)
teststat_mat_bg <- sapply(1:len, function(i){
  if(i %% floor(len/10) == 0) cat('*')
  if(length(result_list_bg[[i]]) == 0) return(rep(NA,2))
  teststat_func(result_list_bg[[i]])
})
colnames(teststat_mat_bg) <- names(result_list_bg)
rownames(teststat_mat_bg) <- c("Mid", "Far")


png(paste0("../../../../out/figures/Writeup6g/Writeup6g_day0_binomial-teststat_", treatment, ".png"),
    height = 1200, width = 2500, res = 300, units = "px")
par(mfrow = c(1,2))
hist(teststat_mat_bg[1,], breaks = 50, xlab = "Chi-sq binomial test stat.", main = "Mid vs. Close (DABTRAM)")
rug(teststat_mat[1,], col = "firebrick", lwd = 2)

hist(teststat_mat_bg[2,], breaks = 50, xlab = "Chi-sq binomial test stat.", main = "Far vs. Close (DABTRAM)")
rug(teststat_mat[2,], col = "firebrick", lwd = 2)
graphics.off()

####################

logpmat <- -log10(pvalue_mat)
idx <- which(logpmat >= 2, arr.ind = T)
round(logpmat[,sort(unique(idx[,2]))], 3)

gene_vec <- colnames(logpmat[,sort(unique(idx[,2]))])

##############

library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
treatment <- "DABTRAM"

surviving_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("day10_", treatment)] >= 20)]
dying_lineages <- rownames(tab_mat)[which(apply(tab_mat,1,max)<=1)]
winning_idx <- intersect(
  intersect(which(all_data$assigned_lineage %in% surviving_lineages),
            which(all_data$assigned_posterior >= 0.5)),
  which(all_data$dataset == "day0")
)
dying_idx <- intersect(
  intersect(which(all_data$assigned_lineage %in% dying_lineages),
            which(all_data$assigned_posterior >= 0.5)),
  which(all_data$dataset == "day0")
)
winning_cells <- colnames(all_data)[winning_idx]
dying_cells <- colnames(all_data)[dying_idx]

keep_vec <- rep(NA, ncol(all_data))
keep_vec[winning_idx] <- paste0("day0_winner_", treatment)
keep_vec[dying_idx] <- paste0("day0_loser_", treatment)
keep_vec <- as.factor(keep_vec)
table(keep_vec)
all_data2 <- all_data
all_data2$keep <- keep_vec
all_data2 <- subset(all_data2, keep %in% c(paste0("day0_winner_", treatment),
                                           paste0("day0_loser_", treatment)))

Seurat::DefaultAssay(all_data2) <- "ATAC"
Seurat::Idents(all_data2) <- "keep"

pdf(paste0("../../../../out/figures/Writeup6g/Writeup6g_day0-coverage_", treatment, "_5000bp_select-genes.pdf"),
    onefile = T, width = 9, height = 4.5)
for(gene in gene_vec){
  plot1 <- Signac::CoveragePlot(
    object = all_data2,
    region = gene,
    features = gene,
    extend.upstream = 5000,
    extend.downstream = 5000
  )
  
  print(plot1)
}

dev.off() 

result_list[["DABTRAM"]][gene_vec]
