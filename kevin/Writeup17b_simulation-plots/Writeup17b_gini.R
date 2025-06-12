rm(list=ls())
library(Seurat)
library(multiomeFate)

# first: Plastic

out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup14/Writeup14_plastic-setting_"
load(paste0(out_folder, "simulation.RData"))

cell_lineage <- as.character(simulation_res$lineage_assignment)
uniq_lineage <- sort(unique(cell_lineage))
lineage_future_count <- simulation_res$lineage_future_size
dineq::gini.wtd(lineage_future_count) # 0.28
# from the documentation: 
# The Gini coefficient is a measure of inequality among values of a distribution. 
# The most used single measure for income inequality. The coefficient can theoretically 
# range between 0 and 1, with 1 being the highest possible inequality (for instance: 
# 1 person in a society has all income; the others none). But coefficients that are 
# negative or greater than 1 are also possible because of negative values in the 
# distribution. Compared to other measures of inequality, the Gini coefficient is 
# especially sensitive for changes in the middle of the distribution.

###############

rm(list=ls())
out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup14/Writeup14_priming-setting_"
load(paste0(out_folder, "simulation_v2.RData"))

cell_lineage <- as.character(simulation_res$lineage_assignment)
uniq_lineage <- sort(unique(cell_lineage))
lineage_future_count <- simulation_res$lineage_future_size
dineq::gini.wtd(lineage_future_count) # 0.2657663


