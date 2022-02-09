rm(list=ls())
load("final_image.RData")
library(Seurat); library(devtools)

all_data
head(all_data@meta.data)
table(all_data$OG_condition)
table(all_data$OG_condition_naivesplit)
table(all_data$Size_Lin)

coding_notes <- "This is a file extracted out by Kevin after tinkering with the data provided by Dylan (which was provided on 2021-10-28)"
date_of_run <- Sys.time()
session_info <- devtools::session_info()

save(all_data, coding_notes, date_of_run, session_info, 
     file = "final_all_data.RData")

################

load("../../../../data/Sydney_stressors/2021-10-28/2021_10_28_Cleaned_Data/final_all_data.RData")

all_data2 <- subset(all_data, final_lineage %nin% c("No Barcode", "Still multiple", "Too large Naive") )
Seurat::DefaultAssay(all_data2) <- "RNA"
all_data2[["SCT"]] <- NULL
all_data2[["umap"]] <- NULL
all_data2[["lineage"]] <- NULL

lineage_vec <- sort(unique(all_data2$final_lineage))
n <- ncol(all_data2)
j_vec <- 1:n
i_vec <- sapply(1:n, function(i){
  which(lineage_vec == all_data2$final_lineage[i])
})
all(1:length(lineage_vec) %in% i_vec)
lin_mat <- Matrix::sparseMatrix(i = i_vec, j = j_vec, 
                                x = rep(1, length(i_vec)),
                                dims = c(length(lineage_vec), n))
rownames(lin_mat) <- lineage_vec
colnames(lin_mat) <- rownames(all_data2@meta.data)
all_data2[["lineage"]] <- Seurat::CreateAssayObject(counts = lin_mat)

all_data2 <- Seurat::SCTransform(all_data2)
all_data2 <- Seurat::RunPCA(all_data2, verbose = FALSE)
all_data2 <- Seurat::RunUMAP(all_data2, dims = 1:50, verbose = FALSE)

plot1 <- Seurat::DimPlot(all_data2, 
                         reduction = "umap", 
                         group.by = "OG_condition", 
                         label = TRUE, repel = TRUE, label.size = 2.5) + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1 <- plot1 + ggplot2::ggtitle("Colored by treatment")
ggplot2::ggsave(filename = "../../../../out/figures/Writeup3f/Writeup3f_syndey_umap.png",
                plot1, device = "png", width = 6, height = 5, units = "in")

####

cell_names <- rownames(all_data2@meta.data)
p0 <- Seurat::DimPlot(all_data2, reduction = "umap", 
                      cells.highlight = cell_names[which(all_data2$OG_condition=="naive")])+ ggplot2::theme(legend.position="none") + ggplot2::ggtitle("Naive")
p1 <- Seurat::DimPlot(all_data2, reduction = "umap", 
                      cells.highlight = cell_names[which(all_data2$OG_condition=="cocl2")])+ ggplot2::theme(legend.position="none") + ggplot2::ggtitle("Cocl2")
p2 <- Seurat::DimPlot(all_data2, reduction = "umap", 
                      cells.highlight = cell_names[which(all_data2$OG_condition=="cis")])+ ggplot2::theme(legend.position="none") + ggplot2::ggtitle("Cis")
p3 <- Seurat::DimPlot(all_data2, reduction = "umap", 
                      cells.highlight = cell_names[which(all_data2$OG_condition=="acid")])+ ggplot2::theme(legend.position="none") + ggplot2::ggtitle("Acid")
p4 <- Seurat::DimPlot(all_data2, reduction = "umap", 
                      cells.highlight = cell_names[which(all_data2$OG_condition=="dab")])+ ggplot2::theme(legend.position="none") + ggplot2::ggtitle("Dab")
p5 <- Seurat::DimPlot(all_data2, reduction = "umap", 
                      cells.highlight = cell_names[which(all_data2$OG_condition=="tram")])+ ggplot2::theme(legend.position="none") + ggplot2::ggtitle("Tram")
p <- cowplot::plot_grid(p0, p1, p2, p3, p4, p5)
cowplot::save_plot(filename = "../../../../out/figures/Writeup3f/Writeup3f_syndey_umap_separate.png", 
                   p, ncol = 3, nrow = 2, base_asp = 1.1, device = "png")

###############

source("../Writeup3d/funcs.R")
source("../Writeup3e/select_cells.R")
tabulate_mat <- .tabulate_lineages(all_data2, condition_var = "OG_condition")
tabulate_mat[1:15,]
quantile(tabulate_mat[,"naive"])
max_val <- sapply(1:nrow(tabulate_mat), function(i){max(tabulate_mat[i,-5])})
tabulate_mat <- tabulate_mat[order(max_val, decreasing = T),]
tabulate_mat[1:15,]
tabulate_mat[nrow(tabulate_mat):1,][1:15,]
tabulate_mat[order(tabulate_mat[,"tram"], decreasing = T),][1:15,]

zz <- rowSums(tabulate_mat[,-which(colnames(tabulate_mat) == "naive")])
quantile(zz)

save(all_data2, coding_notes, date_of_run, session_info, 
     file = "final_all_data2.RData")



