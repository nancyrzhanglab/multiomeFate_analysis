barcode_assignment <- function(seurat_obj,
                               fixing_multiplier = 1,
                               lineage_maximum_size = 10,
                               max_threshold = NULL,
                               threshold_strategy = c("Max_single", "Max_ratio", "Max_difference"),
                               verbose = T){
  stopifnot(inherits(seurat_obj, "Seurat"), "Lineage" %in% Seurat::Assays(seurat_obj))
  if(length(threshold_strategy) > 1) threshold_strategy <- threshold_strategy[1]
  
  counts <- Seurat::GetAssayData(seurat_obj, assay = 'Lineage')
  num_lin_orig <- nrow(counts)
  
  # Only care about lineages that have greater than one total count
  counts_filt <- counts[which(Matrix::rowSums(counts) > 1),]
  counts_filt <- as.matrix(counts_filt)
  n <- ncol(counts_filt)
  num_lin_filt <- nrow(counts_filt)
  percent_lin_filt <- num_lin_filt/num_lin_orig*100
  
  # Find how many barcodes each cell has that has greater than 0 counts
  barcodes_per_cell <- apply(counts_filt,2, function(x) sum(x>0))
  barcodes_per_cell_table <- table(barcodes_per_cell)
  
  if(is.null(max_threshold)) max_threshold <- max(barcodes_per_cell)
  threshold_df <- as.data.frame(t(sapply(1:max_threshold, function(threshold){
    if(verbose) print(paste0("Trying a lineage count threshold of ", threshold))
    num_exceed_vec <- sapply(1:n, function(i){
      length(which(counts_filt[,i] >= threshold))
    })
    
    c(zero = length(which(num_exceed_vec == 0)), 
      one = length(which(num_exceed_vec == 1)),
      two_or_more = length(which(num_exceed_vec > 1)))
  })))
  
  # Evaluate different ways of filtering ----
  # Find best way to set the threshold by trying a couple of different metrics
  strategy_vec <- rep(NA, 3)
  names(strategy_vec) <- c("Max_single", "Max_ratio", "Max_difference")
  threshold_vec <- as.numeric(rownames(threshold_df))
  strategy_vec[1] <- threshold_vec[which.max(threshold_df$one)]
  strategy_vec[2] <- threshold_vec[which.max(threshold_df$one/threshold_df$two_or_more)]
  strategy_vec[3] <- threshold_vec[which.max(threshold_df$one - threshold_df$two_or_more)]
  
  # filter based on decided threshold
  threshold_val <- strategy_vec[threshold_strategy]
  selected_lineage_list <- lapply(1:n, function(i){
    rownames(counts_filt)[which(counts_filt[,i] >= threshold_val)]
  })
  
  # See how many cells still have reads from how many barcodes
  num_lineages_table <- table(lengths(selected_lineage_list))
  
  # Build the indices of the cells that still have multiple barcodes - max_single
  multi_lineage_idx <- which(lengths(selected_lineage_list)>1)
  
  # Assign dominant barcode if the number of counts in the highest expressed barcodes is >fixing_multiplier that of second highest barcode
  # Replace the values in fam that had 2+ barcodes with with either the dominant barcode if found or say that its still multiple 
  for(i in multi_lineage_idx){ 
    topn_vec <- kit::topn(counts_filt[,i],2)
    if(counts_filt[topn_vec[1],i] > fixing_multiplier*counts_filt[topn_vec[2],i]){
      selected_lineage_list[[i]] <- rownames(counts_filt)[topn_vec[1]]
    }else{
      selected_lineage_list[[i]] <- "Still multiple"
    }
  }
  
  # Give cells without a barcode the label 'No barcode'
  selected_lineage_list[lengths(selected_lineage_list) == 0] <- "No barcode"
  
  # Format into metadata so that final lineage determinations can be added to the seurat object
  selected_lineage_vec <- unlist(selected_lineage_list)
  names(selected_lineage_vec) <- colnames(counts_filt)
  seurat_obj$final_lineage <- selected_lineage_vec
  
  # Find out the size of lineages in cells
  cells_per_lineage  <- table(seurat_obj$final_lineage[!seurat_obj$final_lineage %in% c("No barcode", "Still multiple")])
  
  # Set a maximum lineage size and see how many barcoded cells remain in the naive condition ----
  # See total number of lineages in each case 
  lineages_large <- names(cells_per_lineage[cells_per_lineage > lineage_maximum_size])
  lineages_lessequal <- names(cells_per_lineage[cells_per_lineage <= lineage_maximum_size])
  
  # Lets see how many naive cells were in these large lineages 
  num_cells_lineages_large <- length(which(seurat_obj$final_lineage %in% lineages_large))
  num_cells_lineages_lessequal <-  length(which(seurat_obj$final_lineage %in% lineages_lessequal))
  
  # Make cells in the large lineage into a "too large Naive"
  seurat_obj$final_lineage[seurat_obj$final_lineage %in% lineages_large] <- "Too large Naive"
  
  # Write an excel file of the number of the number of cells/lineages above/below each threshold
  output_list <-  list(Percent_filtered_lineages = percent_lin_filt,
                       Num_lineages_per_cell = num_lineages_table,
                       Num_large_lineages = length(lineages_large),
                       Num_cells_in_large_lineages = num_cells_lineages_large,
                       Num_lessEqual_lineages = length(lineages_lessequal),
                       Num_cells_in_lessEqual_lineages = num_cells_lineages_lessequal)
  final_vec <- c(
    Assigned_lineage = length(which(!seurat_obj$final_lineage %in% c("No barcode", "Still multiple", "Too large Naive"))),
    No_barcode = length(which(seurat_obj$final_lineage == "No barcode")),
    Still_multiple = length(which(seurat_obj$final_lineage == "Still multiple")),
    Too_large_Naive = length(which(seurat_obj$final_lineage == "Too large Naive"))
  )
  
  param_list <- list(fixing_multiplier = fixing_multiplier,
                     lineage_maximum_size = lineage_maximum_size,
                     max_threshold = max_threshold,
                     threshold_strategy = threshold_strategy)
  
  list(seurat_obj = seurat_obj,
       threshold_df = threshold_df,
       strategy_vec = strategy_vec,
       output_list = output_list,
       cells_in_lineagetypes = final_vec,
       param_list = param_list)
  
}

plot_barcode_threshold <- function(threshold_df,
                                   main = "Number of cells remaining after different\nlineage barcode count cutoffs"){
  p <- nrow(threshold_df)
  tmp_df <- data.frame(
    count_filt = rep(as.numeric(rownames(threshold_df)),3),
    num_cells = as.numeric(unlist(threshold_df)),
    barcode_num = c(rep('Zero lineage',p),rep('One lineage',p), rep('Multiple lineages',p))
  )
  
  plot1 <- ggplot2::ggplot(data = tmp_df, 
                           ggplot2::aes(x = count_filt, y = num_cells, group = barcode_num, color = barcode_num)) 
  plot1 <- plot1 + ggplot2::geom_point() + ggplot2::geom_line() + ggplot2::labs(title = main)
  plot1
}

plot_lineage_barcodecounts <- function(seurat_obj,
                                       main = "For each lineage, average barcode count\nper cell vs log2(number of cells)",
                                       alpha = 0.05,
                                       height_jitter = 0.2,
                                       width_jitter = 0.1){
  stopifnot("final_lineage" %in% colnames(seurat_obj@meta.data))
  
  counts <- Seurat::GetAssayData(seurat_obj, assay = 'Lineage')
  num_lin_orig <- nrow(counts)
  
  # Only care about lineages that have greater than one total count
  counts_filt <- counts[which(Matrix::rowSums(counts) > 1),]
  counts_filt <- as.matrix(counts_filt)
  
  # Make a scatterplot of the number of cells in each lineage vs ave number of counts ----
  # build a named vector of 0s of the length of the number of lineages 
  uniq_lineage <- unique(seurat_obj$final_lineage)
  uniq_lineage <- uniq_lineage[!uniq_lineage %in% c("No barcode", "Still multiple", "Too large Naive")]
  
  tmp_lineage_df <- data.frame(t(sapply(uniq_lineage, function(lineage_name){
    cell_idx <- which(seurat_obj$final_lineage == lineage_name)
    lineage_idx <- which(rownames(counts_filt) == lineage_name)
    
    num_cells <- length(cell_idx)
    total_count <- sum(counts_filt[lineage_idx, cell_idx])
    
    c(num_cells = num_cells, total_count = total_count)
  })))
  
  tmp_df <- data.frame(Lineage = uniq_lineage, 
                       Log2_Number_cells = log2(tmp_lineage_df$num_cells),
                       Total_counts = tmp_lineage_df$total_count,
                       Average_counts = tmp_lineage_df$total_count/tmp_lineage_df$num_cells)
  tmp_df$Log2_Number_cells <- tmp_df$Log2_Number_cells + stats::runif(length(tmp_df$Log2_Number_cells), min = -width_jitter, max = width_jitter)
  tmp_df$Average_counts <- tmp_df$Average_counts + stats::runif(length(tmp_df$Average_counts), min = -height_jitter, max = height_jitter)
  
  # Make scatterplot of average number of counts vs lineage size and save
  plot1 <- ggplot2::ggplot(tmp_df, 
                           ggplot2::aes(x = Log2_Number_cells, y = Average_counts))
  plot1 <- plot1 + ggplot2::geom_point(alpha = alpha) + ggplot2::labs(title = main) 
  
  plot1
}
