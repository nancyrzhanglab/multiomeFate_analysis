rm(list=ls()); set.seed(10)

library(Seurat); library(Signac); library(EnsDb.Hsapiens.v86); 
library(BSgenome.Hsapiens.UCSC.hg38); library(GenomeInfoDb)
library(dplyr); library(ggplot2); library(multiomicCCA)

mbrain <- readRDS("../../../../data/Embrian_mouse/obj_seurat.rds")

# port over labels
mbrain_jane <- readRDS("../../../../data/Embrian_mouse/data_tenx_labels_Jane.rds")
cell_id <- rownames(mbrain@meta.data)
idx <- sapply(1:length(cell_id), function(x){
  which(rownames(mbrain_jane@meta.data) == cell_id[x])[1]
})
mbrain@meta.data$celltype <- mbrain_jane$savercatLable[idx]

dim(mbrain[["SCT"]]@data[Seurat::VariableFeatures(mbrain, assay = "SCT"),])
dim(mbrain[["ATAC"]]@data)

#################

getRegion <- function(vicinity.size, TSS){
  upstream.size <- vicinity.size[1]
  downstream.size <- vicinity.size[2]
  if (upstream.size == 0) {
    upstream <- TSS
  } else if (upstream.size > 0) {
    upstream <- GenomicRanges::flank(TSS, width=upstream.size, start=T) # upstream excluding the TSS
  }
  if (downstream.size == 0) {
    downstream <- TSS
  } else if (downstream.size > 0) {
    downstream <- GenomicRanges::flank(TSS, width=downstream.size, start=F) # downstream excluding the TSS
  }
  vicinity <- GenomicRanges::punion(upstream, downstream, fill.gap=T)
  BiocGenerics::start(vicinity) <- pmax(0, BiocGenerics::start(vicinity))
  vicinity$gene_name <- TSS$gene_name
  return(vicinity)
}

##############################

seurat_obj <- mbrain
gene_names <- rownames(mbrain[["SCT"]]@data[Seurat::VariableFeatures(mbrain, assay = "SCT"),])
p2 <- length(gene_names)
peaks_gr <- GenomicRanges::granges(seurat_obj@assays$ATAC)
annotations <- seurat_obj@assays$ATAC@annotation

# first find all the relevant genes and peaks that we'll need
relevant_peaks <- lapply(1:p2, function(j){
  if(j %% floor(p2/10) == 0) cat('*')
  
  temp <- Signac::LookupGeneCoords(seurat_obj, gene_names[j])
  if(is.null(temp)) return(numeric(0))
  GenomeInfoDb::seqlevelsStyle(temp)="UCSC"  # just so that the seqlevelsStyle matches up with peaks_gr
  temp2 <- annotations[which(annotations$gene_name==gene_names[j] & annotations$type=="cds")] #cds stands for "coding sequence"
  if(length(temp2)==0) return(numeric(0))
  gene_strand <- BiocGenerics::strand(temp2)[1]
  TSS_position <- ifelse(gene_strand == "+", BiocGenerics::start(temp), BiocGenerics::end(temp))
  TSS <- GenomicRanges::GRanges(seqnames = GenomeInfoDb::seqnames(temp),
                                ranges = IRanges::IRanges(start = TSS_position, width = 1),
                                strand = BiocGenerics::strand(temp),
                                gene_name = gene_names[j])
  upstream_gr <- getRegion(c(50000, 0), TSS)
  downstream_gr <- getRegion(c(0, 50000), TSS)
  
  # populate the peaks counts array.
  peaks_upstream <- GenomicRanges::countOverlaps(peaks_gr, upstream_gr)>0
  peaks_downstream <- GenomicRanges::countOverlaps(peaks_gr, downstream_gr)>0
  
  c(sort(unique(c(which(peaks_upstream), which(peaks_downstream)))))
})

quantile(sapply(relevant_peaks, length))
length(which(sapply(relevant_peaks, length) == 0))
table(table(unlist(relevant_peaks)))

# reparameterize the data
gene_reparam_mat <- data.frame(org_idx = which(sapply(relevant_peaks, length) > 0))
gene_reparam_mat$new_idx <- 1:nrow(gene_reparam_mat)

peak_reparam_mat <- data.frame(org_idx = sort(unique(unlist(relevant_peaks))))
peak_reparam_mat$new_idx <- 1:nrow(peak_reparam_mat)

relevant_peaks2 <- relevant_peaks[gene_reparam_mat$org_idx]
for(i in 1:length(relevant_peaks2)){
  relevant_peaks2[[i]] <- peak_reparam_mat$new_idx[which(peak_reparam_mat$org_idx %in% relevant_peaks2[[i]])]
}

mat_y <- t(mbrain[["SCT"]]@data[Seurat::VariableFeatures(mbrain, assay = "SCT"),])[,gene_reparam_mat$org_idx]
mat_x <- t(mbrain[["ATAC"]]@data)[,peak_reparam_mat$org_idx]
dim(mat_x); dim(mat_y)
p1 <- ncol(mat_x); p2 <- ncol(mat_y)

# form the hash mapping
ht_map <- hash::hash()
for(i in 1:p2){
  ht_map[[as.character(i)]] <- relevant_peaks2[[i]]
}

# remove too-outliery cell-types
idx <- which(mbrain@meta.data$celltype %in% c("Oligodendrocyte", "Hindbrain glycinergic", "Midbrain glutamatergic",
                                              "Forebrain glutamatergic", "Radial glia", "Forebrain GABAergic", 
                                              "Neuroblast", "Mixed region GABAergic", 
                                              "Glioblast", "Cortical or hippocampal glutamatergic"))
mat_x <- mat_x[idx,]; mat_y <- mat_y[idx,]
celltype <- mbrain@meta.data$celltype[idx]

save(mat_x, mat_y, ht_map, celltype, relevant_peaks2, 
     file = "../../../../out/kevin/Writeup3b/10x_embryo.RData")

#####################

rm(list=ls()); gc(T)
load("../../../../out/kevin/Writeup3b/10x_embryo.RData")

vec_start <- which(celltype == "Radial glia")
list_end <- list(which(celltype == "Oligodendrocyte"),
                 which(celltype == "Forebrain GABAergic"),
                 which(celltype == "Cortical or hippocampal glutamatergic"))
(length(vec_start) + length(unlist(list_end)))/nrow(mat_x)

rank_x <- 30
rank_y <- 50
df_x <- data.frame(name = colnames(mat_x))
df_y <- data.frame(name = colnames(mat_y))
mat_x <- as.matrix(mat_x)
mat_y <- as.matrix(mat_y)

set.seed(10)
prep_obj <- multiomeFate::chromatin_potential_prepare(mat_x, mat_y, df_x, df_y, 
                                        vec_start, list_end,
                                        form_method = "average_weighted",
                                        est_method = "threshold_glmnet",
                                        rec_method = "distant_cor",
                                        ht_map = ht_map,
                                        options = list(nn_nn = 30, dim_nlatent_x = rank_x,
                                                       dim_nlatent_y = rank_y, 
                                                       est_num_iterations = 4,
                                                       rec_bool_pred_nn = T,
                                                       est_cv_choice = "lambda.min",
                                                       form_stepsize = 0.5,
                                                       form_min_weight = 0,
                                                       est_verbose = T,
                                                       rec_verbose = T))

set.seed(10)
res <- multiomeFate::chromatin_potential(prep_obj, verbose = T, bool_oracle = F)
