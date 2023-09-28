rm(list=ls())
library(Seurat)
library(Signac)

load("../../../../out/kevin/Writeup6m/Writeup6m_all-data.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

cis_fasttopics <- all_data[["fasttopic_CIS"]]
cocl2_fasttopics <- all_data[["fasttopic_COCL2"]]
dabtram_fasttopics <- all_data[["fasttopic_DABTRAM"]]

note <- paste0("From Writeup4e_timeAll_peakmerging_combining_part2.R")

save(date_of_run, session_info,
     cis_fasttopics,
     cocl2_fasttopics,
     dabtram_fasttopics,
     file = "../../../../out/kevin/Writeup6n/Writeup6n_topics.RData")
