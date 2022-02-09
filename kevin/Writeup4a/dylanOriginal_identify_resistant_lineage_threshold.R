# Use single-cell data to set a threshold for resistant lineages in gDNA

# Initialize ----
rm(list = ls())
gc()
setwd('/Volumes/GoogleDrive/My Drive/Schaff_Shared/Cloud/Experiment_IDs/DLS005/Segmented_Analysis_Pipeline/preprocess_cDNA_BCs/identify_resistant_lineage_threshold')
library(dplyr)
library(Seurat)
library(xlsx)
library(ggplot2)
library(RColorBrewer)
library(kit)
`%nin%` = Negate(`%in%`)

# Load in the all_data object with the final lineages assigned ----
load('/Volumes/GoogleDrive/My Drive/Schaff_Shared/Cloud/Experiment_IDs/DLS005/Segmented_Analysis_Pipeline/preprocess_cDNA_BCs/assign_dominant_barcodes/all_data_final_lineages.RData')

# Load in the gDNA counts ----
load('/Volumes/GoogleDrive/My Drive/Schaff_Shared/Cloud/Experiment_IDs/DLS005/Segmented_Analysis_Pipeline/preprocess_gDNA_BCs/preprocessed_gDNA.RData')

# for each resistant condition, pull out the number of cells in each lineage ---- 
cocl2_lineages <- data.frame(table(all_data$Lineage_noLarge[all_data$OG_condition == 'cocl2']))
cocl2_lineages <- cocl2_lineages[cocl2_lineages$Var1 %nin% c('No Barcode', 'Still multiple', 'Lineage too large in Naive cells'),]
colnames(cocl2_lineages) <- c('Lineage','Num_Cells')
cocl2_df <- merge(cocl2_lineages,gDNA_anno$sc_output_counts_DLS005_CoCl2_gDNA_S35.txt, by = 'Lineage' )
cocl2_lineages$RPM_gDNA <- sort(as.numeric(gDNA_anno$sc_output_counts_DLS005_CoCl2_gDNA_S35.txt$RPM[cocl2_lineages$Lineage]),decreasing = T)


acid_lineages <- data.frame(table(all_data$Lineage_noLarge[all_data$OG_condition == 'acid']))
acid_lineages <- acid_lineages[acid_lineages$Var1 %nin% c('No Barcode', 'Still multiple', 'Lineage too large in Naive cells'),]
colnames(acid_lineages) <- c('Lineage','Num_Cells')
acid_df <- merge(acid_lineages,gDNA_anno$sc_output_counts_DLS005_Acid_gDNA_S36.txt, by = 'Lineage' )
acid_lineages$RPM_gDNA <- sort(as.numeric(gDNA_anno$sc_output_counts_DLS005_Acid_gDNA_S36.txt$RPM[acid_lineages$Lineage]),decreasing = T)


cis_lineages <- data.frame(table(all_data$Lineage_noLarge[all_data$OG_condition == 'cis']))
cis_lineages <- cis_lineages[cis_lineages$Var1 %nin% c('No Barcode', 'Still multiple', 'Lineage too large in Naive cells'),]
colnames(cis_lineages) <- c('Lineage','Num_Cells')
cis_df <- merge(cis_lineages,gDNA_anno$sc_output_counts_DLS005_Cis_gDNA_S40.txt, by = 'Lineage' )
cis_lineages$RPM_gDNA <- sort(as.numeric(gDNA_anno$sc_output_counts_DLS005_Cis_gDNA_S40.txt$RPM[cis_lineages$Lineage]),decreasing = T)


dox_lineages <- data.frame(table(all_data$Lineage_noLarge[all_data$OG_condition == 'dox']))
dox_lineages <- dox_lineages[dox_lineages$Var1 %nin% c('No Barcode', 'Still multiple', 'Lineage too large in Naive cells'),]
colnames(dox_lineages) <- c('Lineage','Num_Cells')
dox_df <- merge(dox_lineages,gDNA_anno$sc_output_counts_DLS005_Dox_gDNA_S39.txt, by = 'Lineage' )
dox_lineages$RPM_gDNA <- sort(as.numeric(gDNA_anno$sc_output_counts_DLS005_Dox_gDNA_S39.txt$RPM[dox_lineages$Lineage]),decreasing = T)


dab_lineages <- data.frame(table(all_data$Lineage_noLarge[all_data$OG_condition == 'dab']))
dab_lineages <- dab_lineages[dab_lineages$Var1 %nin% c('No Barcode', 'Still multiple', 'Lineage too large in Naive cells'),]
colnames(dab_lineages) <- c('Lineage','Num_Cells')
dab_df <- merge(dab_lineages,gDNA_anno$sc_output_counts_DLS005_Dab_gDNA_S37.txt, by = 'Lineage' )
dab_lineages$RPM_gDNA <- sort(as.numeric(gDNA_anno$sc_output_counts_DLS005_Dab_gDNA_S37.txt$RPM[dab_lineages$Lineage]),decreasing = T)


tram_lineages <- data.frame(table(all_data$Lineage_noLarge[all_data$OG_condition == 'tram']))
tram_lineages <- tram_lineages[tram_lineages$Var1 %nin% c('No Barcode', 'Still multiple', 'Lineage too large in Naive cells'),]
colnames(tram_lineages) <- c('Lineage','Num_Cells')
tram_df <- merge(tram_lineages,gDNA_anno$sc_output_counts_DLS005_Tram_gDNA_S38.txt, by = 'Lineage' )
tram_lineages$RPM_gDNA <- sort(as.numeric(gDNA_anno$sc_output_counts_DLS005_Tram_gDNA_S38.txt$RPM[tram_lineages$Lineage]),decreasing = T)

# Export xlsx file with minimum RPM cutoffs for resistant lineages ----
output_df <- data.frame( Condition = c('Acid','Cis','Cocl2','Dab','Dox','Tram'),
                         cutoff_1cell = c(mean(acid_df$RPM[acid_df$Num_Cells == 1]),
                           mean(cis_df$RPM[cis_df$Num_Cells == 1]),
                           mean(cocl2_df$RPM[cocl2_df$Num_Cells == 1]),
                           mean(dab_df$RPM[dab_df$Num_Cells == 1]),
                           mean(dox_df$RPM[dox_df$Num_Cells == 1]),
                           mean(tram_df$RPM[tram_df$Num_Cells == 1])),
                         cutoff_2cell = c(mean(acid_df$RPM[acid_df$Num_Cells == 2]),
                                          mean(cis_df$RPM[cis_df$Num_Cells == 2]),
                                          mean(cocl2_df$RPM[cocl2_df$Num_Cells == 2]),
                                          mean(dab_df$RPM[dab_df$Num_Cells == 2]),
                                          mean(dox_df$RPM[dox_df$Num_Cells == 2]),
                                          mean(tram_df$RPM[tram_df$Num_Cells == 2])))

write.xlsx(output_df, file = 'Resistant_lineage_RPM_cutoffs.xlsx')
