rm(list = ls())

library(multiomeFate)
library(Seurat)
library(UCell)
library(tidyverse)
library(data.table)
library(ggplot2)
library(ggpubr)
library(gridExtra)

figure_dir <- "~/project/Multiome_fate/git/multiomeFate_analysis_kevin/fig/kevin/Writeup18/"

results_dir <- "~/project/Multiome_fate/out/emilia/task2_correlate_fate_potential_and_features_V2/"

all_data <- multiomeFate:::data_loader(which_files = c("fasttopics", "saver"))
