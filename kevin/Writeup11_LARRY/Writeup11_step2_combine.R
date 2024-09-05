rm(list=ls())
library(Seurat)
library(multiomeFate)

out_folder <- "~/project/Multiome_fate/out/kevin/Writeup11/"
load(paste0(out_folder, "Writeup11_larry_seurat.RData"))
load(paste0(out_folder, "Writeup11_larry_fasttopics.RData"))

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

##########

# preparing
cell_names <- Seurat::Cells(seurat_obj)
K <- ncol(topic_res$L)
topic_mat <- matrix(NA, nrow = length(cell_names), ncol = K)
rownames(topic_mat) <- cell_names

for(i in 1:nrow(topic_res$L)){
  topic_mat[rownames(topic_res$L)[i],] <- topic_res$L[i,]
}
colnames(topic_mat) <- paste0("fastTopic_", 1:K)
colnames(topic_res$F) <- paste0("fastTopic_", 1:K)

seurat_obj[["fasttopic"]] <- Seurat::CreateDimReducObject(embeddings = topic_mat, 
                                                          loadings =  topic_res$F,
                                                          assay = "RNA")

seurat_obj <- subset(seurat_obj, state_info %in% c("Monocyte", "Neutrophil", "Undifferentiated"))

tmp <- cbind(as.character(seurat_obj$state_info), 
             as.character(seurat_obj$time_info))
seurat_obj$time_celltype <- apply(tmp, 1, function(x){
  paste0(x[1], "-", x[2])
})
seurat_obj$time_celltype <- factor(seurat_obj$time_celltype)

#######

table(seurat_obj$time_celltype)

save(seurat_obj, date_of_run, session_info,
     file = paste0(out_folder, "Writeup11_larry_seurat_fasttoics.RData"))

print("Done! :)")
