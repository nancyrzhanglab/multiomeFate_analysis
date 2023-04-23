rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)

treatment <- "COCL2"

load(paste0("../../../../out/kevin/Writeup6e/Writeup6e_", treatment, "_peak-entropy.RData"))
load(paste0("../../../../out/kevin/Writeup6d/Writeup6d_coverage_pileup_", treatment, ".RData"))
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
gene_vec <- names(result_list)

######################################
# TEMPORARY FIX:
load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")

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

num_frags <- sapply(gene_vec, function(gene){
  print(gene)
  cutmat_winning <- multiomeFate:::extract_cutmatrix(
    object = all_data,
    gene = gene,
    cells = winning_cells
  )
  cutmat_dying <- multiomeFate:::extract_cutmatrix(
    object = all_data,
    gene = gene,
    cells = dying_cells
  )
  
  c(winning_num_frags = length(cutmat_winning@x),
    dying_num_frags = length(cutmat_dying@x))
})

for(gene in gene_vec){
  result_list[[gene]]$num_frag_winning <- num_frags["winning_num_frags",gene]
  result_list[[gene]]$num_frag_dying <- num_frags["dying_num_frags",gene]
}

######################################

for(kk in 1:ceiling(length(gene_vec)/5)){
  print(kk)
  gene_vec_small <- gene_vec[((kk-1)*5+1):min((kk*5), length(gene_vec))]
  
  png(paste0("../../../../out/figures/Writeup6e/Writeup6e_", treatment, "_coverage-pileup_", kk, ".png"),
      width = 1800, height = 3000, units = "px", res = 300)
  par(mfrow = c(5,3), mar = c(4,4,6,0.5))
  for(j in 1:length(gene_vec_small)){
    par(mar = c(4,4,5,0.5))
    gene <- gene_vec_small[j]
    ymax <- max(winning_curve_list_master[[gene]]$pileup_vec, 
                dying_curve_list_master[[gene]]$pileup_vec)
    
    coverage_pileup_plotter(
      pileup_res = winning_curve_list_master[[gene]],
      col_curve = 3,
      main = paste0(gene, " Winning (", treatment, "),\nn=", result_list[[gene]]$num_frag_winning),
      ylim = c(0,ymax)
    )
    
    coverage_pileup_plotter(
      pileup_res = dying_curve_list_master[[gene]], 
      col_curve = 2,
      main = paste0(gene, " Losing (", treatment, "),\nn=", result_list[[gene]]$num_frag_dying),
      ylim = c(0,ymax)
    )
    
    if(entropy_mat["win",gene] > entropy_mat["die",gene]){
      starwin <- "*"; starlose <- ""
    } else {
      starwin <- ""; starlose <- "*"
    }
    par(mar = rep(0.5, 4))
    plot(NA, xlim = c(0,1), ylim = c(0,1), 
         xaxt='n', yaxt = "n", bty = "n", ann=FALSE, 
         xlab = "", ylab = "")
    text(0.5, 0.5, labels = paste0(
      "Win, Entropy: ", round(entropy_mat["win",gene], 2), starwin,
      "\n(", paste0(round(result_list[[gene]]$theta_winning*100), collapse = ", "),
      ")\n\nDie, Entropy: ", round(entropy_mat["die",gene], 2), starlose,
      "\n(", paste0(round(result_list[[gene]]$theta_dying*100), collapse = ", "),
      ")\n\n-Log10 pval = ",round(-log10(result_list[[gene]]$pval+1e-6),2)
    ))
  }
  
  graphics.off()
}

