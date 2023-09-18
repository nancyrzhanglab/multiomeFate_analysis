rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(IRanges)

## see https://www.r-bloggers.com/2020/03/what-is-a-dgcmatrix-object-made-of-sparse-matrix-format-in-r/
# if you want to find the nonzero entries for a row, I suggest
# first transposing via Matrix::t()
.nonzero_col <- function(mat, col_idx, bool_value){
  stopifnot(inherits(mat, c("dgCMatrix", "lgCMatrix")), col_idx %% 1 == 0,
            col_idx > 0, col_idx <= ncol(mat))
  
  val1 <- mat@p[col_idx]
  val2 <- mat@p[col_idx+1]
  
  if(val1 == val2) return(numeric(0))
  if(bool_value){
    # return the value
    mat@x[(val1+1):val2]
  } else {
    # return the row index
    mat@i[(val1+1):val2]+1
  }
}

load("../../../../out/kevin/Writeup6m/Writeup6m_all-data.RData")
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
load("../../../../out/kevin/Writeup6l/Writeup6l_COCL2-split-by-day10_differential-peak.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

# simplifying the dataset
print("Simplifying the dataset")
Seurat::DefaultAssay(all_data) <- "ATAC"
all_data[["geneActivity"]] <- NULL
all_data[["pca"]] <- NULL
all_data[["umap"]] <- NULL
all_data[["saverpca"]] <- NULL
all_data[["fasttopic_CIS"]] <- NULL
all_data[["fasttopic_COCL2"]] <- NULL
all_data[["fasttopic_DABTRAM"]] <- NULL
all_data[["common_tcca"]] <- NULL
all_data[["distinct1_tcca"]] <- NULL
all_data[["distinct2_tcca"]] <- NULL
keep_vec <- rep(FALSE, ncol(all_data))
keep_vec[which(all_data$dataset == "day0")] <- TRUE
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == TRUE)

peak_mat <- all_data[["ATAC"]]@counts
peak_mat_t <- Matrix::t(peak_mat)

# set up the motif of interest
print("Extracting relevant motifs")
motif_focus <- "FOS"
motif_matrix <- Signac::GetMotifData(object = all_data, slot = "data")
motif_name_vec <- unlist(all_data[["ATAC"]]@motifs@motif.names)
motif_idx <- grep(motif_focus, motif_name_vec)

# extract the positive and negative peaks
print("Finding positive/negative peaks")
idx <- intersect(which(de_res[,"p_val"] <= 1e-2), 
                 which(de_res[,"avg_log2FC"] > 0))
positive_peaks <- rownames(de_res)[idx]

idx <- intersect(which(de_res[,"p_val"] <= 1e-2), 
                 which(de_res[,"avg_log2FC"] < 0))
negative_peaks <- rownames(de_res)[idx]

# grab the relevant cells
print("Partitioning cells")
treatment <- "COCL2"
lineage_names_win <- rownames(tab_mat)[which(tab_mat[,paste0("day10_", treatment)] >= 10)]
cell_names_win <- colnames(all_data)[which(all_data$assigned_lineage %in% lineage_names_win)]
lineage_names_lose <- rownames(tab_mat)[which(tab_mat[,paste0("day10_", treatment)] == 0)]
cell_names_lose <- colnames(all_data)[which(all_data$assigned_lineage %in% lineage_names_lose)]
ident_vec <- rep(NA, ncol(all_data))
names(ident_vec) <- colnames(all_data)
ident_vec[cell_names_win] <- "winner"
ident_vec[cell_names_lose] <- "loser"
table(ident_vec)
winner_idx <- which(ident_vec == "winner")
loser_idx <- which(ident_vec == "loser")

#################################

print("Plotting")
pdf(paste0("../../../../out/figures/Writeup6m/Writeup6m_day0_motif-peak_split-by-day10COCL2_scatterplot_", motif_focus, ".png"),
    onefile = T, width = 8, height = 5)

for(kk in motif_idx[1:2]){
  # grab the relevant peaks
  motif <- motif_name_vec[kk]
  peak_idx <- .nonzero_col(motif_matrix, 
                           col_idx = which(colnames(motif_matrix) == names(motif)),
                           bool_value = F)
  
  # compute the relevant percentage
  percentage_mat <- sapply(peak_idx, function(j){
    cell_idx <- .nonzero_col(mat = peak_mat_t,
                             col_idx = j,
                             bool_value = F)
    winner_percentage <- length(intersect(winner_idx, cell_idx))/length(winner_idx)
    loser_percentage <- length(intersect(loser_idx, cell_idx))/length(loser_idx)
    
    c(winner = winner_percentage, 
      loser = loser_percentage)
  })
  colnames(percentage_mat) <- colnames(peak_mat_t)[peak_idx]
  
  positive_idx <- which(colnames(percentage_mat) %in% positive_peaks)
  negative_idx <- which(colnames(percentage_mat) %in% negative_peaks)
  status_vec <- rep("none", ncol(percentage_mat))
  status_vec[positive_idx] <- "positive"
  status_vec[negative_idx] <- "negative"
  status_vec <- factor(status_vec)
  
  # now make the plot
  xmax <- max(percentage_mat)
  value_vec <- c("black", "dodgerblue3",  "red")
  names(value_vec) <- c("none", "positive", "negative")
  
  df <- data.frame(winner = percentage_mat["winner",],
                   loser = percentage_mat["loser",],
                   status = status_vec)
  
  df <- df[c(which(df[,"status"] == "none"),  which(df[,"status"] == "positive"), which(df[,"status"] == "negative")),]
  
  p1 <- ggplot2::ggplot(df, ggplot2::aes(x = winner, 
                                         y = loser))
  p1 <- p1 +  ggplot2::geom_point( ggplot2::aes(color = status, 
                                                size = status, 
                                                alpha = ifelse(status == "none", 0.5, 1)))  # Set color and size based on 'status'
  p1 <- p1 + ggplot2::scale_color_manual(values = c("none" = "black", "positive" = "dodgerblue3", "negative" = "red"))  # Customize point colors
  p1 <- p1 + ggplot2::scale_size_manual(values = c("none" = 2, "positive" = 4, "negative" = 4))  # Customize point sizes
  p1 <- p1 + ggplot2::guides(color = ggplot2::guide_legend(title = "Status"))   # Add a legend title
  p1 <- p1 + ggplot2::labs(size = "Status")   # Label for size aesthetic
  p1 <- p1 + ggplot2::xlim(0, xmax) + ggplot2::ylim(0, xmax)
  p1 <- p1 + ggplot2::geom_abline(intercept = 0, 
                                  slope = 1, 
                                  linewidth = 2,
                                  linetype = "dashed", 
                                  color = "green")  # Add a dashed green x=y line
  p1 <- p1 + ggplot2::labs(y= "Loser peak percentage", x = "Winner peak percentage") 
  p1 <- p1 + ggplot2::ggtitle(paste0(motif, ": Percentage among winner or loser cells, defined by Day10 COCL2\n",
                                     length(winner_idx), " winner cells, ", length(loser_idx), " loser cells\n",
                                     length(positive_peaks), " positive peaks, ", length(negative_peaks), " negative peaks, ",  length(peak_idx), " total peaks")) 
  
  withr::with_tempfile({
    png_file <- paste0("temp_plot_", kk, ".png")
    ggplot2::ggsave(png_file, plot = p1, type = "cairo-png", dpi = 300)
    
    # Insert the rasterized PNG into the PDF
    grid::grid.raster(png_file)
  })
  
}

dev.off()

print("Done! :)")



