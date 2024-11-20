rm(list=ls())
library(Seurat)
library(Signac)
library(multiomeFate)

plot_folder <- "~/project/Multiome_fate/git/multiomeFate_analysis_kevin/fig/kevin/Writeup15/"

#######

all_data <- multiomeFate:::data_loader(which_files = "fasttopics")

treatment_vec <- c("CIS", "COCL2", "DABTRAM")


.compute_cumulative_vec <- function(vec){
  vec <- sort(vec, decreasing = FALSE)
  n <- length(vec)
  cumsum_vec <- cumsum(vec)
  cumsum_vec <- cumsum_vec/max(cumsum_vec)
  
  list(x = seq(0, 1, length.out = n+1), y = c(0, cumsum_vec))
}

for(kk in 1:2){
  if(kk == 1){
    day_early <- "day10"
    day_later <- "week5"
    cell_threshold <- 4
  } else {
    day_early <- "day0"
    day_later <- "day10"
    cell_threshold <- 10
  }
  
  
  for(treatment in treatment_vec){
    day_early_full <- ifelse(day_early == "day0", day_early, paste0(day_early, "_", treatment))
    day_early_abrv <- ifelse(day_early == "day0", "d0", "d10")
    day_late_abrv <- ifelse(day_later == "day10", "d10", "w5")
    fate_variable <- paste0("fatepotential_", treatment, "_", day_early_abrv, "_", day_late_abrv)
    cell_idx <- which(all_data$dataset == day_early_full)
    cell_early <- all_data@meta.data[cell_idx, c("assigned_lineage", fate_variable)]
    
    # compute the global gini curve
    future_cell <- as.numeric(10^cell_early[,fate_variable])
    global_gini <- .compute_cumulative_vec(future_cell)
    
    # compute the local gini curve
    tab_vec <- table(cell_early[,"assigned_lineage"])
    lineage_names <- names(tab_vec)[tab_vec >= cell_threshold]
    lineage_gini <- lapply(lineage_names, function(lineage){
      cell_early_subset <- cell_early[cell_early$assigned_lineage == lineage,]
      future_cell_subset <- as.numeric(10^cell_early_subset[,fate_variable])
      .compute_cumulative_vec(future_cell_subset)
    })
    names(lineage_gini) <- lineage_names
    
    # reorder lineage sizes
    future_size_vec <- sapply(lineage_names, function(lineage){
      idx <- which(cell_early$assigned_lineage == lineage)
      sum(as.numeric(10^cell_early[idx,fate_variable]))
    })
    names(future_size_vec) <- lineage_names
    lineage_order <- order(future_size_vec, decreasing = FALSE)
    future_size_vec <- future_size_vec[lineage_order]
    lineage_gini <- lineage_gini[lineage_order]
    
    # create the scales
    n_colors <- 10
    color_vec <- viridisLite::viridis(n_colors+1)[-(n_colors+1)]
    break_color_vec <- seq(min(future_size_vec), max(future_size_vec), length = n_colors)
    
    n_size <- 5
    size_vec <- seq(0.5, 1.5, length.out = n_size)
    break_size_vec <- seq(log(min(future_size_vec)), log(max(future_size_vec)), length = n_size)
    
    n_alpha <- 5
    alpha_vec <- seq(0.2, 1, length.out = n_alpha)
    break_alpha_vec <- seq(min(future_size_vec), max(future_size_vec), length = n_alpha)
    
    # plot
    png(paste0(plot_folder, "Writeup15_predicted-gini_", treatment, "_", day_early_abrv, "-", day_late_abrv,".png"),
        height = 1200, width = 1200, units = "px", res = 300)
    plot(NA, 
         xlab = paste0("% of ", day_early, " cells"), 
         ylab = paste0("% of expansion at ", day_later, " (Predicted)"),
         xlim = c(0,1),
         ylim = c(0,1),
         asp = TRUE,
         main = paste0("Predicted Gini: ", treatment, " (", day_early_abrv, "-", day_late_abrv, ")"))
    lines(x = global_gini$x,
          y = global_gini$y,
          lwd = 3)
    for(i in 1:length(lineage_gini)){
      col = color_vec[which.min(abs(break_color_vec-future_size_vec[i]))]
      alpha = alpha_vec[which.min(abs(break_alpha_vec-future_size_vec[i]))]
      col = scales::alpha(col, alpha = alpha)
      lwd = size_vec[which.min(abs(break_size_vec-log(future_size_vec[i])))]
      lines(x = lineage_gini[[i]]$x,
            y = lineage_gini[[i]]$y,
            lwd = lwd,
            col = col)
    }
    graphics.off()
  }
  
}





