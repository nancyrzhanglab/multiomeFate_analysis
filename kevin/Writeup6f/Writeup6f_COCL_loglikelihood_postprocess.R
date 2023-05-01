rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

treatment <- "COCL2"

# find the winning and losing cells
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
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
length(winning_idx); length(dying_idx)
winning_cells <- colnames(all_data)[winning_idx]
dying_cells <- colnames(all_data)[dying_idx]

##############################

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

############################

lrt_vec_all <- lrt_vec_all[!is.na(lrt_vec_all)]
lrt_vec_all <- lrt_vec_all[sort(names(lrt_vec_all))]

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
ub <- stats::quantile(z_vec, probs = 0.95)
idx <- which(z_vec <= ub)
N0 <- length(idx)
z0_vec <- z_vec[idx]

init_df <- mean(z0_vec)
init_theta <- 0.95 * stats::pchisq(ub, df = init_df)

optim_res <- stats::optim(par = c(init_df, init_theta),
                          fn = fn,
                          method = "Nelder-Mead",
                          N = N,
                          N0 = N0,
                          z0_vec = z0_vec,
                          ub = ub)
optim_res$par
df_est <- optim_res$par[1]

pvalue_vec <- sapply(lrt_vec_all, function(x){
  1-stats::pchisq(x, df = df_est)
})
idx <- which(pvalue_vec <= 0.05)

round(sort(lrt_vec_all[idx], decreasing = T))

############################
# now let's make some coverage tracks

keep_vec <- rep(NA, ncol(all_data))
keep_vec[colnames(all_data) %in% winning_cells] <- "day0_win"
keep_vec[colnames(all_data) %in% dying_cells] <- "day0_lose"
all_data$keep <- keep_vec
all_data <- subset(all_data, keep %in% c("day0_win", "day0_lose"))

Seurat::DefaultAssay(all_data) <- "ATAC"
Seurat::Idents(all_data) <- "keep"

gene_vec <- sort(names(idx))

pdf(paste0("../../../../out/figures/Writeup6f/Writeup6f_coverage_", treatment, "_1000bp_highLRT.pdf"),
    onefile = T,
    width = 9, height = 4.5)
for(gene in gene_vec){
  plot1 <- Signac::CoveragePlot(
    object = all_data,
    region = gene,
    features = gene,
    extend.upstream = 1000,
    extend.downstream = 1000
  )
  print(plot1)
}
dev.off() 

#######################################

for(gene in gene_vec){
  if(!all(is.na(result_list_total[[gene]]))){
    obj_win <- result_list_total[[gene]]$res_win$grenander_obj
    obj_die <- result_list_total[[gene]]$res_die$grenander_obj
  } else {
    obj_win <- result_list[[gene]]$res_win$grenander_obj
    obj_die <- result_list[[gene]]$res_die$grenander_obj
  }
  
  area_vec1 <- cumsum(diff(obj_win$x)*obj_win$pdf[-length(obj_win$pdf)])
  area_vec2 <- cumsum(diff(obj_win$x)*obj_win$pdf[-length(obj_win$pdf)])
  
  xmax1 <- obj_win$x[which.min(abs(area_vec1 - 0.9))]
  xmax2 <- obj_win$x[which.min(abs(area_vec2 - 0.9))]
  xmax <- max(xmax1, xmax2); xlim <- c(0, xmax)
  ylim <- c(0, max(c(obj_win$pdf, obj_die$pdf)))
  
  png(paste0("../../../../out/figures/Writeup6f/Writeup6f_COCL_entropy_decr-density_", gene, ".png"),
      height = 1500, width = 2500, units = "px", res = 300)
  par(mfrow = c(1,2), mar = c(4,4,4,0.5))
  plot_grenander(obj = obj_win,
                 xlim = xlim, ylim = ylim, 
                 main = paste0(gene, ": Win"))
  plot_grenander(obj = obj_die,
                 xlim = xlim, ylim = ylim, 
                 main = paste0(gene, ": Die"))
  graphics.off()
  
}
