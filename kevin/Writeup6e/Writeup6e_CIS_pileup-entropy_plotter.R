rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)

load("../../../../out/kevin/Writeup6e/Writeup6e_CIS_peak-entropy.RData")
load("../../../../out/kevin/Writeup6d/Writeup6d_coverage_pileup_CIS.RData")
source("../Writeup6d/coverage_extractor_singlecell.R")
source("../Writeup6d/coverage_extractor_singlecell-plotter.R")
source("../Writeup6d/coverage_pileup.R")
source("../Writeup6d/coverage_pileup_plotter.R")

.entropy <- function(x, tol = 1e-6){
  idx <- which(x <= tol)
  tmp <- sapply(x, function(i){i*log(i)})
  tmp[idx] <- 0
  -sum(tmp)
}

entropy_mat <- sapply(result_list, function(x){
  tmp <- c(.entropy(x$theta_winning), .entropy(x$theta_dying))
  names(tmp) <- c("win", "die")
  tmp
})

all(names(result_list) == names(winning_curve_list_master))

png(paste0("../../../../out/figures/Writeup6d/Writeup6d_", treatment, "_coverage-pileup_", i, ".png"),
    width = 1200, height = 3000, units = "px", res = 300)
par(mfrow = c(5,2))
for(j in 1:length(gene_vec_small)){
  ymax <- max(winning_curve_list[[j]]$pileup_vec, dying_curve_list[[j]]$pileup_vec)
  
  coverage_pileup_plotter(
    pileup_res = winning_curve_list[[j]],
    col_curve = 3,
    main = paste0(gene_vec_small[j], " Winning (", treatment, ")"),
    ylim = c(0,ymax)
  )
  
  coverage_pileup_plotter(
    pileup_res = dying_curve_list[[j]], 
    col_curve = 2,
    main = paste0(gene_vec_small[j], " Losing (", treatment, ")"),
    ylim = c(0,ymax)
  )
}

graphics.off()

