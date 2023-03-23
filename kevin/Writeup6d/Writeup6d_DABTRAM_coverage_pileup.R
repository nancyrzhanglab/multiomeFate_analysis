rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")
source("coverage_extractor_singlecell.R")
source("coverage_extractor_singlecell-plotter.R")
source("coverage_pileup.R")
source("coverage_pileup_plotter.R")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

treatment <- "DABTRAM"

source("../Writeup6b/gene_list.R")
source("gene_list_csc.R")
gene_vec <- sort(unique(c(unlist(keygenes), keygenes_csc)))

# remove any genes not found
tmp <- sapply(gene_vec, function(gene){
  tmp <- Signac::LookupGeneCoords(
    object = all_data,
    gene = gene,
    assay = "ATAC"
  )
  # make sure gene exists
  if(all(is.null(tmp))) return(T) else return(F)
})
if(any(tmp)){
  gene_vec <- gene_vec[-which(tmp)]
}

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

############################
winning_curve_list_master <- numeric(0)
dying_curve_list_master <- numeric(0)

for(i in 1:ceiling(length(gene_vec)/5)){
  print(paste0("Batch number ", i, " out of ", ceiling(length(gene_vec)/5)))
  gene_vec_small <- gene_vec[((i-1)*5+1):min(i*5, length(gene_vec))]
 
  winning_curve_list <- vector("list", length(gene_vec_small))
  names(winning_curve_list) <- gene_vec_small
  dying_curve_list <- winning_curve_list
  
  for(j in 1:length(gene_vec_small)){
    print(paste0("Working on gene: ", gene_vec_small[j]))
    
    # first winning cells
    print("Working on winning cells")
    tmp <- coverage_pileup(
      object = all_data,
      gene = gene_vec_small[j],
      cells = winning_cells
    )
    winning_curve_list[[j]] <- compute_pileup_curve(
      pileup_mat = tmp$pileup_mat,
      peak_width_max = tmp$peak_width_max,
      peak_width_median = tmp$peak_width_median
    )
    
    # next dying cells
    print("Working on dying cells")
    tmp <- coverage_pileup(
      object = all_data,
      gene = gene_vec_small[j],
      cells = dying_cells
    )
    dying_curve_list[[j]] <- compute_pileup_curve(
      pileup_mat = tmp$pileup_mat,
      peak_width_max = tmp$peak_width_max,
      peak_width_median = tmp$peak_width_median
    )
  }
  
  # save outputs
  winning_curve_list_master <- c(winning_curve_list_master, winning_curve_list)
  dying_curve_list_master <- c(dying_curve_list_master, dying_curve_list)
  save(date_of_run, session_info,
       winning_curve_list_master, dying_curve_list_master,
       treatment,
       file = paste0("../../../../out/kevin/Writeup6d/Writeup6d_coverage_pileup_", treatment, ".RData"))
  
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
}



