rm(list=ls())
load("Writeup6p_all-data_lightweight_noATAC.RData")

ls()
all_data 
# contains: assays (RNA, Saver -- the "denoised" RNA, and lineage. Note: The ATAC assay is not present for filesize reasons)
# dimension reductions: 
## -- from "raw" RNA: pca, umap (technically -- after log-normalized and rescaled)
## -- from "denoised" RNA (i.e., SAVER): saverpca, saverumap
## -- from "raw" ATAC: atac.umap  (technically -- after TF-IDF)
## -- The following dimension reductions are only of the relevant day0, day10, and week5 cells of the that treatment
## --- based on RNA, constructed via fastTopics: fasttopic_CIS, fasttopic_COCL2, fasttopic_DABTRAM
## --- based on ATAC, constructed via peackVI: peakVI_CIS, peakVI_COCL2, peakVI_DABTRAM

session_info
class(all_data) # object was created via Seurat 4.4.0

# lineage data
head(all_data@meta.data) ## the important columns is "assigned_lineage"
head(all_data$assigned_lineage)
tab_mat <- table(all_data$assigned_lineage, all_data$dataset) # table of how many lineages are in each condition
head(tab_mat)

# RNA data
all_data[["Saver"]]
head(Seurat::VariableFeatures(all_data[["Saver"]]))
mat <- all_data[["Saver"]]$data
dim(mat)
mat[1:5,1:5]
