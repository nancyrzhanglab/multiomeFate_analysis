rm(list=ls())
library(Seurat)
library(fastTopics)
library(ggplot2)
library(GGally)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

in_dir <- '~/Dropbox/Thesis/Lineage_trace/data/Watermelon/'
out_dir <- '~/Dropbox/Thesis/Lineage_trace/data/Watermelon/'

# ==============================================================================
# Read data
# ==============================================================================

seurat_object <- readRDS(paste0(in_dir, 'pc9_time_course.rds'))
metadat <- seurat_object@meta.data

# ==============================================================================
# Basic info
# ==============================================================================
metadat$lineage_barcode_present <- ifelse(is.na(metadat$lineage_barcode), 'No', 'Yes')
total_cells <- metadat %>% 
  group_by(time_point) %>% 
  summarise(n = n())

barcode_recovery <- metadat %>% 
  group_by(lineage_barcode_present, time_point) %>% 
  summarise(num_cells = n())

# Recovery
barcode_recovery$time_point <- as.factor(barcode_recovery$time_point)
ggplot(barcode_recovery) +
  geom_col(aes(x = time_point, y = num_cells, fill = lineage_barcode_present), position="stack", width = 0.5) +
  scale_fill_manual(values = c('#C9C9C9', '#01204E')) +
  theme_bw()

# Detection
all_barcodes <- unique(metadat$lineage_barcode)
barcode_detection <- metadat %>% 
  filter(!is.na(lineage_barcode)) %>% 
  group_by(time_point) %>% 
  summarise(n_lineage = n_distinct(lineage_barcode))
barcode_detection$total_barcodes <- length(all_barcodes)

barcode_detection$time_point <- as.factor(barcode_detection$time_point)
ggplot(barcode_detection) +
  geom_col(aes(x = time_point, y = n_lineage), fill = '#01204E', width = 0.5) +
  geom_hline(yintercept = length(all_barcodes)) +
  theme_bw()
  

# Lineage size
lineage_size <- metadat %>% 
  group_by(lineage_barcode, time_point) %>% 
  summarise(num_cells = n())
lineage_size <- merge(lineage_size, total_cells, by = c('time_point'))
lineage_size$lin_freq <- lineage_size$num_cells / lineage_size$n
lineage_size$log10_n <- log10(lineage_size$num_cells + 1)
lineage_size <- lineage_size[, c('lineage_barcode', 'time_point', 'log10_n')]
lineage_size <- lineage_size %>% drop_na()
lineage_size_w <- dcast(lineage_size, lineage_barcode ~ time_point)
lineage_size_w[is.na(lineage_size_w)] <- 0
rownames(lineage_size_w) <- lineage_size_w$lineage_barcode
lineage_size_w <- subset(lineage_size_w, select = -c(lineage_barcode))


ggpairs(lineage_size_w) + 
  theme_bw()


# ==============================================================================
# Basic info (D3, 7, 14)
# ==============================================================================
metadat <- metadat[metadat$lineage_barcode_present == 'Yes', ]
d0 <- metadat[metadat$time_point == 0, ]
d3 <- metadat[metadat$time_point == 3, ]
d7 <- metadat[metadat$time_point == 7, ]
d14 <- metadat[metadat$time_point == 14, ]
d14 <- d14[d14$lineage_barcode %in% d7$lineage_barcode, ]

intersect_d3_d7_d14 <- intersect(d3$lineage_barcode, d7$lineage_barcode)
intersect_d3_d7_d14 <- intersect(intersect_d3_d7_d14, d14$lineage_barcode)

nrow(d3[d3$lineage_barcode %in% intersect_d3_d7_d14, ])
nrow(d7[d7$lineage_barcode %in% intersect_d3_d7_d14, ])
nrow(d14[d14$lineage_barcode %in% intersect_d3_d7_d14, ])

# ==============================================================================
# Basic preprocessing
# ==============================================================================

Seurat::VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# day14 <- subset(seurat_object, idents = 14)
# seurat_object <- day14
seurat_object <- Seurat::NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_object <- Seurat::FindVariableFeatures(seurat_object,
                                              selection.method = "vst",
                                              nfeatures = 2000)
seurat_object <- Seurat::ScaleData(seurat_object)
seurat_object <- Seurat::RunPCA(seurat_object,
                                features = Seurat::VariableFeatures(object = seurat_object),
                                verbose = F)
seurat_object <- Seurat::RunUMAP(seurat_object,
                                 dims = 1:30)


# plot based on time point
Seurat::Idents(seurat_object) <- "time_point"
p1 <- Seurat::DimPlot(seurat_object,
                      reduction = "umap",
                      label = T, label.box = T)
p1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 8)) 

seurat_object_subset <- subset(seurat_object, 
                               select = lineage_barcode %in% all_barcodes)

Seurat::Idents(seurat_object) <- "sample_type"
p2 <- Seurat::DimPlot(seurat_object,
                      reduction = "umap",
                      label = T, label.box = T)
p2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 8)) 


# ==============================================================================
# Basic preprocessing
# ==============================================================================
mat <- Seurat::GetAssayData(object = seurat_object, 
                            assay = "RNA", 
                            layer = "counts")
mat <- mat[Seurat::VariableFeatures(seurat_object),]
sum_vec <- Matrix::rowSums(mat)
if(any(sum_vec <= 0.3)){
  mat <- mat[sum_vec > 0.3, ]
}
mat <- Matrix::t(mat)

K <- 30
set.seed(10)
topic_res <- fastTopics::fit_topic_model(mat, k = K)

#########

topic_mat <- matrix(NA, nrow = ncol(seurat_object), ncol = ncol(topic_res$L))
rownames(topic_mat) <- colnames(seurat_object)
topic_res$L <- topic_res$L[rownames(topic_res$L) %in% colnames(seurat_object),]

for(i in 1:nrow(topic_res$L)){
  topic_mat[rownames(topic_res$L)[i],] <- topic_res$L[i,]
}
colnames(topic_mat) <- paste0("fastTopic_", 1:ncol(topic_mat))
colnames(topic_res$F) <- paste0("fastTopic_", 1:ncol(topic_mat))

seurat_object[["fasttopic"]] <- Seurat::CreateDimReducObject(embeddings = topic_mat, 
                                                             loadings =  topic_res$F,
                                                             assay = "RNA",
                                                             key =  "fastTopic_")

print("Finished")
date_of_run <- Sys.time()
session_info <- devtools::session_info()
save(seurat_object, date_of_run, session_info,
     file = paste0(out_dir, "/PC9_time_course_fasttopics.RData"))


