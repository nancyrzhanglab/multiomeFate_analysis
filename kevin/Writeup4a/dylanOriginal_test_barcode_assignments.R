# Script for testing different ways to assigning the final barcodes to each cell

# Initialize ----
rm(list = ls())
setwd('/Volumes/GoogleDrive/My Drive/Schaff_Shared/Cloud/Experiment_IDs/DLS005/Segmented_Analysis_Pipeline/preprocess_cDNA_BCs/test_barcode_assignments')
library(dplyr)
library(Seurat)
library(xlsx)
library(ggplot2)
library(RColorBrewer)
library(kit)

# Load in the merged dataset for all data ----
load('/Volumes/GoogleDrive/My Drive/Schaff_Shared/Cloud/Experiment_IDs/DLS005/Segmented_Analysis_Pipeline/preprocess_GEX/all_data_merged.RData')
reference <- read.delim('/Volumes/GoogleDrive/My Drive/Schaff_Shared/Super_cells/2022_01_11_DLS005_PrepBarcodesCellRanger/CellRanger_inputs/FeatureReference_filtered.csv', sep = ',')

# Build counts object  ----
counts <- GetAssayData(all_data, assay = 'lineage')
num_lin_orig <- nrow(counts)

# Remove any lineages that only have one total count across all cells ----
counts_filt <- counts[which(rowSums(counts) > 1),] # Only care about lineages that have greater than one total count
num_lin_filt <- nrow(counts_filt)
percent_lin_filt <- num_lin_filt/num_lin_orig*100

df_counts <- as.data.frame(counts_filt)

# Find how many barcodes each cell has that has greater than 0 counts
barcodes_per_cell <- apply(counts_filt,2, function(x) sum(x>0))
barcodes_per_cell_table <- table(barcodes_per_cell)
barcodes_per_cell_table

# Run loop to find how many cells per lineage in naive cells in each condition ----
condition <- c()
for (i in 1:length(all_data@assays$RNA@counts@Dimnames[2][[1]])){
  condition <- rbind(condition, strsplit(all_data@assays$RNA@counts@Dimnames[2][[1]][[i]], '_')[[1]][1])
}
all_data$OG_condition_naivesplit <- condition


# Pipe this to run as a for loop through each of the different conditions
plot_df_list <- list() 

for (k in unique(all_data$OG_condition_naivesplit)){
  print(k)
  plot_df <- data.frame()
  for (i in 1:50){ # 1:20
    temp_fam_list <- list()
    print(i)
    count <- 1
    for (j in grep(k,colnames(counts_filt))){
      temp_fam_list[[count]] <- which(counts_filt[,j] > i)
      count <- count+1
    }
    table <- table(lengths(temp_fam_list))
    temp <- data.frame(zero = as.numeric(table[1]), one = as.numeric(table[2]), two_or_more = length(temp_fam_list) - as.numeric(table[1]) - as.numeric(table[2]))
    plot_df <- rbind(plot_df,temp)
  }
  plot_df_list[[k]] <- plot_df
}

# Build a properly formatted plot list
plot_df_list2 <- list()
for (i in names(plot_df_list)){
  plot_df_list2[[i]] <- data.frame(
    count_filt = rep(1:50,3),
    num_cells = unlist(plot_df_list[[i]]),
    barcode_num = c(rep('zero',50),rep('one',50), rep('multiple',50))
  )
}

# Write the plots to a PDF
pdf('threshold_sweep.pdf')
for (i in names(plot_df_list2)){
  print(ggplot(data = plot_df_list2[[i]], aes(x = count_filt, y = num_cells, group = barcode_num, color = barcode_num)) + geom_point()+geom_line() + labs(title = i))
}
dev.off()

# Evaluate different ways of filtering ----

# Find best way to set the threshold by trying a couple of different metrics
cutoffs <-data.frame()
for(i in names(plot_df_list2)){
  zero <- plot_df_list2[[i]][plot_df_list2[[i]]$barcode_num == 'zero',]
  one <- plot_df_list2[[i]][plot_df_list2[[i]]$barcode_num == 'one',]
  multiple <- plot_df_list2[[i]][plot_df_list2[[i]]$barcode_num == 'multiple',]
  cutoffs <- rbind(cutoffs, data.frame(Max_single = which.max(one$num_cells), Max_ratio = which.max(one$num_cells/multiple$num_cells), Max_difference = which.max(one$num_cells - multiple$num_cells)))
}
rownames(cutoffs) <- names(plot_df_list2)

# Export the cutoffs so that it can be loaded in for final barcode assignment script
save(cutoffs, file = 'lineage_count_cutoffs.RData')

# filter based on the Max_single
filter_max_single <- cutoffs$Max_single
names(filter_max_single) <- names(plot_df_list)

# filter based on decided threshold
fam_max_single <- list()
for (k in names(filter_max_single)){
  for (i in grep(k,colnames(df_counts))){
    fam_max_single[[i]] <- rownames(df_counts)[which(df_counts[,i] > filter_max_single[k])]
  }
}
# See how many cells still have reads from how many barcodes
table_fam_max_single <- table(lengths(fam_max_single))

# filter based on the Max_difference
filter_max_difference <- cutoffs$Max_difference
names(filter_max_difference) <- names(plot_df_list)

# filter based on decided threshold
fam_max_difference <- list()
for (k in names(filter_max_difference)){
  for (i in grep(k,colnames(df_counts))){
    fam_max_difference[[i]] <- rownames(df_counts)[which(df_counts[,i] > filter_max_difference[k])]
  }
}
# See how many cells still have reads from how many barcodes
table_fam_max_difference <- table(lengths(fam_max_difference))

# Write these results to an xlsx file
write.xlsx(table_fam_max_single, 'filter_stats.xlsx', sheetName = 'Max_single')
write.xlsx(table_fam_max_difference, 'filter_stats.xlsx', sheetName = 'Max_difference', append = T)

# Build a boolean of the cells that still have multiple barcodes - max_single
boolean_multi_lin_single <- lengths(fam_max_single)>1

# Assign dominant barcode if the number of counts in the highest expressed barcodes is >triple that of second highest barcode
fixed_multi_single <- list()
for (i in (1:ncol(df_counts))){ 
  if ( boolean_multi_lin_single[i] == T) {
    if(df_counts[topn(as.numeric(df_counts[,i]),2),i][1] >= 3*df_counts[topn(as.numeric(df_counts[,i]),2),i][2]){
      fixed_multi_single[[i]] <- rownames(df_counts)[topn(as.numeric(df_counts[,i]),2)][1]
    }else{
      fixed_multi_single[[i]] <- "Still multiple"
    }
  } else{ fixed_multi_single[[i]] <- "Single"
  }
}

# Build a boolean of the cells that still have multiple barcodes - max_difference
boolean_multi_lin_difference <- lengths(fam_max_difference)>1

# Assign dominant barcode if the number of counts in the highest expressed barcodes is >triple that of second highest barcode
fixed_multi_difference <- list()
for (i in (1:ncol(df_counts))){ 
  if ( boolean_multi_lin_difference[i] == T) {
    if(df_counts[topn(as.numeric(df_counts[,i]),2),i][1] >= 3*df_counts[topn(as.numeric(df_counts[,i]),2),i][2]){
      fixed_multi_difference[[i]] <- rownames(df_counts)[topn(as.numeric(df_counts[,i]),2)][1]
    }else{
      fixed_multi_difference[[i]] <- "Still multiple"
    }
  } else{ fixed_multi_difference[[i]] <- "Single"
  }
}


# Replace the values in fam that had 2+ barcodes with with either the dominant barcode if found or say that its still multiple - max single
fam_multi_collapse_single <- fam_max_single
fam_multi_collapse_single[boolean_multi_lin_single] <- fixed_multi_single[boolean_multi_lin_single]

# Replace the values in fam that had 2+ barcodes with with either the dominant barcode if found or say that its still multiple - max difference
fam_multi_collapse_difference <- fam_max_difference
fam_multi_collapse_difference[boolean_multi_lin_difference] <- fixed_multi_difference[boolean_multi_lin_difference]

# Give cells without a barcode the label 'No barcode'
fam_multi_collapse_single[lengths(fam_multi_collapse_single) == 0] <- 'No Barcode'
fam_multi_collapse_difference[lengths(fam_multi_collapse_difference) == 0] <- 'No Barcode'

# Get just the naive cells
Idents(all_data) <- all_data$OG_condition
only_naive <- WhichCells(all_data, idents = all_data$OG_condition[all_data$OG_condition %in% 'naive'])
`%nin%` = Negate(`%in%`)

# Format into metadata so that final lineage determinations can be added to the seurat object
lin_meta_single <- t(as.data.frame(fam_multi_collapse_single))
colnames(lin_meta_single) <- 'final_lineage_single'
rownames(lin_meta_single) <- colnames(df_counts)
all_data$final_lineage_single <- lin_meta_single

lin_meta_difference <- t(as.data.frame(fam_multi_collapse_difference))
colnames(lin_meta_difference) <- 'final_lineage_difference'
rownames(lin_meta_difference) <- colnames(df_counts)
all_data$final_lineage_difference <- lin_meta_difference

# Find out the size of lineages in naive cells
naive_cells_per_lineage_single  <- as.data.frame(table(all_data$final_lineage_single[names(all_data$final_lineage_single) %in% only_naive & all_data$final_lineage_single %nin% c('No Barcode','Still multiple')]))
naive_cells_per_lineage_difference  <- as.data.frame(table(all_data$final_lineage_difference[names(all_data$final_lineage_difference) %in% only_naive & all_data$final_lineage_difference %nin% c('No Barcode','Still multiple')]))

# Write these results to an xlsx file
write.xlsx(naive_cells_per_lineage_single[order(naive_cells_per_lineage_single$Freq, decreasing = T),], 'naive_lin_sizes.xlsx', sheetName = 'Max_single')
write.xlsx(naive_cells_per_lineage_difference[order(naive_cells_per_lineage_difference$Freq, decreasing = T),], 'naive_lin_sizes.xlsx', sheetName = 'Max_difference', append = T)

# Set a maximum lineage size and see how many barcoded cells remain in the naive condition ----
lin_size_thresh <-10

# See total number of lineages in each case - max_single
naive_greater_single <- as.character(naive_cells_per_lineage_single$Var1[naive_cells_per_lineage_single$Freq > lin_size_thresh]) # 111 lineages
naive_lessEqual_single <- as.character(naive_cells_per_lineage_single$Var1[naive_cells_per_lineage_single$Freq <= lin_size_thresh]) # 952 lineages

# Lets see how many naive cells were in these large lineages - max_single
only_naive_largeLin_single <- names(all_data$final_lineage_single[all_data$OG_condition %in% 'naive' & all_data$final_lineage_single %in% naive_greater_single]) # 2315 cells
only_naive_smallLin_single <- names(all_data$final_lineage_single[all_data$OG_condition %in% 'naive' & all_data$final_lineage_single %in% naive_lessEqual_single]) # 2395 cells
# the small lineages make up slightly more cells than the large lineages

# Make cells in the large lineage into a "too large Naive" - max_single
all_data$final_lineage_single[all_data$final_lineage_single %in% naive_greater_single] <- "Too large Naive"

# See total number of lineages in each case - max_difference
naive_greater_difference <- as.character(naive_cells_per_lineage_difference$Var1[naive_cells_per_lineage_difference$Freq > lin_size_thresh]) # 95 lineages
naive_lessEqual_difference <- as.character(naive_cells_per_lineage_difference$Var1[naive_cells_per_lineage_difference$Freq <= lin_size_thresh]) # 937 lineages

# Lets see how many naive cells were in these large lineages - max_difference
only_naive_largeLin_difference <- names(all_data$final_lineage_difference[all_data$OG_condition %in% 'naive' & all_data$final_lineage_difference %in% naive_greater_difference]) # 2015 cells
only_naive_smallLin_difference <- names(all_data$final_lineage_difference[all_data$OG_condition %in% 'naive' & all_data$final_lineage_difference %in% naive_lessEqual_difference]) # 2359 cells
# the small lineages make up slightly more cells than the large lineages

# Make cells in the large lineage into a "too large Naive" - max_difference
all_data$final_lineage_difference[all_data$final_lineage_difference %in% naive_greater_difference] <- "Too large Naive"

# Write an excel file of the number of the number of cells/lineages above/below each threshold
output_df <-  data.frame(Lineages_greater10 = length(naive_greater_single),
                         Cells_lineages_greater10 = length(only_naive_largeLin_single),
                         Lineages_less10 = length(naive_lessEqual_single),
                         Cells_lineages_less10 = length(only_naive_smallLin_single))
output_df <- rbind(output_df,
                   c(length(naive_greater_difference),
                     length(only_naive_largeLin_difference),
                     length(naive_lessEqual_difference),
                     length(only_naive_smallLin_difference)))

rownames(output_df) <- c('Max_single','Max_difference')
write.xlsx(output_df, 'lineage_sizes_maxSingle_vs_maxDifference.xlsx')

# Make a scatterplot of the number of cells in each lineage vs ave number of counts ----

# build a named vector of 0s of the length of the number of lineages - Max_single
total_counts_max_single <- as.vector(matrix(0,1,nrow(naive_cells_per_lineage_single)))
names(total_counts_max_single) <- naive_cells_per_lineage_single$Var1

for (i in 1:length(only_naive)){
  if(lin_meta_single[[i]] %nin% c('Too large Naive','No Barcode','Still multiple')){
    total_counts_max_single[[lin_meta_single[[i]]]] <- total_counts_max_single[[lin_meta_single[[i]]]] + df_counts[lin_meta_single[[i]], colnames(df_counts)[i]]
  }
}

naive_ave_counts_max_single <- data.frame(Lineage = naive_cells_per_lineage_single$Var1, Number_cells = naive_cells_per_lineage_single$Freq,
                                          Total_counts = total_counts_max_single, Average_counts = total_counts_max_single/naive_cells_per_lineage_single$Freq)

# Make scatterplot of average number of counts vs lineage size and save
pdf('MaxSingle_aveCounts_vs_numCells.pdf')
ggplot(naive_ave_counts_max_single, aes(x = log(Number_cells), y = Average_counts)) + geom_point() + labs(title = 'Average counts per cell in lineage vs log(Number of cells in lineage)')
dev.off()

# build a named vector of 0s of the length of the number of lineages - Max_difference
total_counts_max_difference <- as.vector(matrix(0,1,nrow(naive_cells_per_lineage_difference)))
names(total_counts_max_difference) <- naive_cells_per_lineage_difference$Var1

for (i in 1:length(only_naive)){
  if(lin_meta_difference[[i]] %nin% c('Too large Naive','No Barcode','Still multiple')){
    total_counts_max_difference[[lin_meta_difference[[i]]]] <- total_counts_max_difference[[lin_meta_difference[[i]]]] + df_counts[lin_meta_difference[[i]], colnames(df_counts)[i]]
  }
}

naive_ave_counts_max_difference <- data.frame(Lineage = naive_cells_per_lineage_difference$Var1, Number_cells = naive_cells_per_lineage_difference$Freq,
                                              Total_counts = total_counts_max_difference, Average_counts = total_counts_max_difference/naive_cells_per_lineage_difference$Freq)

# Make scatterplot of average number of counts vs lineage size and save
pdf('Maxdifference_aveCounts_vs_numCells.pdf')
ggplot(naive_ave_counts_max_difference, aes(x = log(Number_cells), y = Average_counts)) + geom_point() + labs(title = 'Average counts per cell in lineage vs log(Number of cells in lineage)')
dev.off()

# Export the ave counts in naive cells to an excel spreadsheet
write.xlsx(naive_ave_counts_max_single, 'naive_ave_lin_sizes.xlsx', sheetName = 'Max_single')
write.xlsx(naive_ave_counts_max_difference, 'naive_ave_lin_sizes.xlsx', sheetName = 'Max_difference', append = T)
