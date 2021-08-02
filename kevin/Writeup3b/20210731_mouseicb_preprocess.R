rm(list=ls())
myeloid <- readRDS("../../../../data/ICB_mouse/Object_seurat_2000_each.rds")

library(Seurat); library(Signac)
myeloid
membership_vec <- as.character(myeloid@meta.data$celltype)
membership_vec <- sapply(membership_vec, function(x){
  paste0(strsplit(x, split = "_")[[1]][1:2], collapse = "_")
})
myeloid@meta.data$celltype <- membership_vec
table(myeloid@meta.data$celltype)

########################

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

gene_names <- rownames(myeloid[["SCT"]]@data[Seurat::VariableFeatures(myeloid, assay = "SCT"),])
p2 <- length(gene_names)
peaks_gr <- GenomicRanges::granges(myeloid@assays$ATAC)
annotations <- myeloid@assays$ATAC@annotation

relevant_peaks <- lapply(1:p2, function(j){
  if(j %% floor(p2/10) == 0) cat('*')
  
  temp <- Signac::LookupGeneCoords(myeloid, gene_names[j])
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

# reparameterize the data
gene_reparam_mat <- data.frame(org_idx = which(sapply(relevant_peaks, length) > 0))
gene_reparam_mat$new_idx <- 1:nrow(gene_reparam_mat)

peak_reparam_mat <- data.frame(org_idx = sort(unique(unlist(relevant_peaks))))
peak_reparam_mat$new_idx <- 1:nrow(peak_reparam_mat)

relevant_peaks2 <- relevant_peaks[gene_reparam_mat$org_idx]
for(i in 1:length(relevant_peaks2)){
  relevant_peaks2[[i]] <- peak_reparam_mat$new_idx[which(peak_reparam_mat$org_idx %in% relevant_peaks2[[i]])]
}

cell_idx <- which(myeloid@meta.data$celltype %in% c("B16_dICB", "B16_NT"))
mat_y <- Matrix::t(myeloid[["SCT"]]@data[Seurat::VariableFeatures(myeloid, assay = "SCT"),])[cell_idx,gene_reparam_mat$org_idx]
mat_x <- Matrix::t(myeloid[["ATAC"]]@data)[cell_idx,peak_reparam_mat$org_idx]
dim(mat_x); dim(mat_y)
p1 <- ncol(mat_x); p2 <- ncol(mat_y)

##########

set.seed(10)
pc_x <- irlba::irlba(mat_x, nv = 50)
tmp2 <- multiomeFate:::.mult_mat_vec(pc_x$u[,-1], pc_x$d[-1])
set.seed(10)
mat_umap_x <- Seurat::RunUMAP(tmp2)@cell.embeddings

set.seed(10)
pc_y <- irlba::irlba(mat_y, nv = 50)
tmp2 <- multiomeFate:::.mult_mat_vec(pc_y$u, pc_y$d)
set.seed(10)
mat_umap_y <- Seurat::RunUMAP(tmp2)@cell.embeddings

tmp_combined <- cbind(multiomeFate:::.mult_mat_vec(pc_x$u[,-1], pc_x$d[-1]/pc_x$d[2]),
                      multiomeFate:::.mult_mat_vec(pc_y$u, pc_y$d/pc_y$d[1]))
set.seed(10)
mat_umap_xy <- Seurat::RunUMAP(tmp_combined)@cell.embeddings

save(myeloid, cell_idx, mat_x, mat_y, mat_umap_x, mat_umap_y, mat_umap_xy,
     relevant_peaks2, 
     file = "../../../../out/kevin/Writeup3b/mouseicb_fate_prep.RData")

########

png(file = "../../../../out/figures/Writeup3b/Writeup3b_mouseicb_umap_rna.png",
    height = 3000, width = 3000, res = 300, units = "px")
plot(mat_umap_y[,1], mat_umap_y[,2], asp = T, col = as.numeric(as.factor(myeloid@meta.data$celltype)),
     pch = 16)
graphics.off()

png(file = "../../../../out/figures/Writeup3b/Writeup3b_mouseicb_umap_atac.png",
    height = 3000, width = 3000, res = 300, units = "px")
plot(mat_umap_x[,1], mat_umap_x[,2], asp = T, col = as.numeric(as.factor(myeloid@meta.data$celltype)),
     pch = 16)
graphics.off()

png(file = "../../../../out/figures/Writeup3b/Writeup3b_mouseicb_umap_both.png",
    height = 3000, width = 3000, res = 300, units = "px")
plot(mat_umap_xy[,1], mat_umap_xy[,2], asp = T, col = as.numeric(as.factor(myeloid@meta.data$celltype)),
     pch = 16)
graphics.off()

#####################



