# Assign the dominant lineage (if there) to each cell
# Will be using the count cutoff that gives the largest delta between the number of cells with a single barcode vs multiple barcodes

# Initialize ----
rm(list = ls())
gc()
setwd('/Volumes/GoogleDrive/My Drive/Schaff_Shared/Cloud/Experiment_IDs/DLS005/Segmented_Analysis_Pipeline/preprocess_cDNA_BCs/assign_dominant_barcodes')
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

# Load in cutoffs from the test_barcode_assignments.R script ----
load('/Volumes/GoogleDrive/My Drive/Schaff_Shared/Cloud/Experiment_IDs/DLS005/Segmented_Analysis_Pipeline/preprocess_cDNA_BCs/test_barcode_assignments/lineage_count_cutoffs.RData')

# filter based on the maximum difference between the number of cells with one lineage vs multiple
filter_max_difference <- cutoffs$Max_difference
names(filter_max_difference) <- rownames(cutoffs)

fam_max_difference <- list()
for (k in names(filter_max_difference)){
  for (i in grep(k,colnames(df_counts))){
    fam_max_difference[[i]] <- rownames(df_counts)[which(df_counts[,i] > filter_max_difference[k])]
  }
}

# Look through cells that still have multiple barcodes and see if any have one with much higher expression than the rest ----

# Build a boolean of the cells that still have multiple barcodes
boolean_multi_lin <- lengths(fam_max_difference)>1

# Assign dominant barcode if the number of counts in the highest expressed barcodes is >triple that of second highest barcode
fixed_multi<- list()
for (i in (1:ncol(df_counts))){ 
  if ( boolean_multi_lin[i] == T) {
    if(df_counts[topn(as.numeric(df_counts[,i]),2),i][1] >= 3*df_counts[topn(as.numeric(df_counts[,i]),2),i][2]){
      fixed_multi[[i]] <- rownames(df_counts)[topn(as.numeric(df_counts[,i]),2)][1]
    }else{
      fixed_multi[[i]] <- "Still multiple"
    }
  } else{ fixed_multi[[i]] <- "Single"
  }
}

# Replace the values in fam that had 2+ barcodes with with either the dominant barcode if found or say that its still multiple
fam_multi_collapse <- fam_max_difference
fam_multi_collapse[boolean_multi_lin] <- fixed_multi[boolean_multi_lin]

fam_multi_collapse[lengths(fam_multi_collapse) == 0] <- 'No Barcode'

# Write the fixed lineages into the all_data object metadata ----
lin_meta <- t(as.data.frame(fam_multi_collapse))
colnames(lin_meta) <- 'Lineage'
rownames(lin_meta) <- colnames(df_counts)
all_data$Lineage <- lin_meta

# Remove lineages that have too many cells in the naive cells ----

# setup new metadata list for removing the lineages that are way too large amongst naive cells
all_data$Lineage_noLarge <- all_data$Lineage

# Get just the naive cells
Idents(all_data) <- all_data$OG_condition
only_naive <- WhichCells(all_data, idents = all_data$OG_condition[all_data$OG_condition %in% 'naive'])
`%nin%` = Negate(`%in%`)

naive_cells_per_lineage <- as.data.frame(table(all_data$Lineage[names(all_data$Lineage) %in% only_naive & all_data$Lineage %nin% c('No Barcode','Still multiple')]))

# Set a maximum lineage size and see how many barcoded cells remain in the naive condition
lin_size_threshold <- 10

# See total number of lineages in each case
naive_greater_threshold <- as.character(naive_cells_per_lineage$Var1[naive_cells_per_lineage$Freq > lin_size_threshold])

# Make cells in the large lineage into a "too large Naive"
all_data$Lineage_noLarge[all_data$Lineage_noLarge %in% naive_greater_threshold] <- "Lineage too large in Naive cells"

# Save the updated all_data object ----
save(all_data, file = 'all_data_final_lineages.RData')

# Load in the naive and resistant only objects and add the lineage information ----

# Naive
load('/Volumes/GoogleDrive/My Drive/Schaff_Shared/Cloud/Experiment_IDs/DLS005/Segmented_Analysis_Pipeline/preprocess_GEX/all_naive_merged.RData')
all_naive$Lineage <- all_data$Lineage[grep('naive',names(all_data$Lineage))]
all_naive$Lineage_noLarge <- all_data$Lineage_noLarge[grep('naive',names(all_data$Lineage_noLarge))]
save(all_naive, file = 'all_naive_final_lineages.RData')
rm(all_naive)

# Resistant 
load('/Volumes/GoogleDrive/My Drive/Schaff_Shared/Cloud/Experiment_IDs/DLS005/Segmented_Analysis_Pipeline/preprocess_GEX/all_resistant_merged.RData')
all_resistant$Lineage <- all_data$Lineage[-grep('naive',names(all_data$Lineage))]
all_resistant$Lineage_noLarge <- all_data$Lineage_noLarge[-grep('naive',names(all_data$Lineage_noLarge))]
save(all_resistant, file = 'all_resistant_final_lineages.RData')
rm(all_resistant)