rm(list=ls())
library(Seurat)
library(Signac)
library(multiomeFate)
library(ggplot2)
library(Ckmeans.1d.dp)

out_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/out/Writeup14/Writeup14_priming-setting_"
plot_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/git/multiomeFate_analysis/fig/kevin/Writeup14/tmp_Writeup14_priming-setting_"
func_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/git/multiomeFate_analysis/kevin/Writeup14_simulation/"

source(paste0(func_folder, "func_seurat.R"))

load(paste0(out_folder, "simulation_v2.RData"))
load(paste0(out_folder, "all_data_fasttopics.RData"))
