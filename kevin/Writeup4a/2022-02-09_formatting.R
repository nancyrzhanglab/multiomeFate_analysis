rm(list=ls())
library(Seurat)
barcode_csv <- read.csv("../../../BarcodeOutputs/2022_02/Cellranger_output/2022_02_Time0/outs/filtered_feature_bc_matrix/barcodes.tsv", header = F)
features_csv <- read.csv("../../../BarcodeOutputs/2022_02/Cellranger_output/2022_02_Time0/outs/filtered_feature_bc_matrix/features.tsv", 
                         sep = "\t", header = F)
gene_vec <- features_csv[features_csv[,3]=="Gene Expression",2]
lineage_vec <- features_csv[features_csv[,3]=="Custom",2]

duplicated_genes <- gene_vec[which(duplicated(gene_vec))]
for(gene_name in duplicated_genes){
  idx <- which(gene_vec == gene_name)
  for(i in idx){
    gene_vec[i] <- paste0(features_csv[i,2], "_", features_csv[i,1])
  } 
}
any(duplicated(gene_vec))
gene_vec[grep("ENSG0", gene_vec)]

mat <- Matrix::readMM("../../../BarcodeOutputs/2022_02/Cellranger_output/2022_02_Time0/outs/filtered_feature_bc_matrix/matrix.mtx")
mat1 <- mat[features_csv[,3]=="Gene Expression",]
mat2 <- mat[features_csv[,3]=="Custom",]
  
colnames(mat1) <- barcode_csv[,1]
colnames(mat2) <- barcode_csv[,1]
rownames(mat1) <- gene_vec
rownames(mat2) <- lineage_vec

t0_obj <- Seurat::CreateSeuratObject(counts = mat1, assay = "RNA")
t0_obj[["Lineage"]] <- Seurat::CreateAssayObject(counts = mat2)

t0_obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(object = t0_obj, pattern = "^MT-")
t0_obj[["percent.rb"]] <- Seurat::PercentageFeatureSet(object = t0_obj, pattern = "^RPS")
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
t0_obj <- Seurat::CellCycleScoring(t0_obj, 
                                   s.features = s.genes, 
                                   g2m.features = g2m.genes)

head(t0_obj@meta.data)

save(t0_obj, file = "../../../out/kevin/Writeup4a/2022-02-09_time0_formatted.RData")

################3

idx <- grep("^MT-", rownames(t0_obj[["RNA"]]@counts))
round(Matrix::rowMeans(t0_obj[["RNA"]]@counts[idx,]),2)
round(quantile(Matrix::rowMeans(t0_obj[["RNA"]]@counts)),2)

