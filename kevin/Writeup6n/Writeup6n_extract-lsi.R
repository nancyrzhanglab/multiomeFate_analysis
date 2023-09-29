rm(list=ls())
library(Seurat)
library(Signac)

load("../../../../out/kevin/Writeup6m/Writeup6m_all-data.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

lsi <- all_data[["lsi"]]

save(date_of_run, session_info,
     lsi,
     file = "../../../../out/kevin/Writeup6n/Writeup6n_lsi.RData")
