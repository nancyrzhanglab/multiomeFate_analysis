rm(list=ls())
library(Seurat)

load("../../../../out/kevin/Writeup6p/Writeup6p_all-data_lightweight_noATAC.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

set.seed(10)
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
tab_mat <- tab_mat[,c("day0", "day10_COCL2", "week5_COCL2")]
stats::cor(log10(tab_mat+1))
n <- nrow(tab_mat); p <- ncol(tab_mat)
tab_mat <- tab_mat + matrix(stats::runif(n*p, max = 0.5), nrow = n, ncol = p)
tab_mat <- log10(tab_mat+1)
tab_mat <- apply(as.matrix.noquote(tab_mat),2,as.numeric)
tab_mat <- as.data.frame(tab_mat)

max_val <- max(tab_mat)

plot1 <- ggplot2::ggplot(tab_mat,  ggplot2::aes(x=day0, y=day10_COCL2)) + 
  ggplot2::geom_point(alpha = 0.3, shape = 16) + ggplot2::xlim(0, max_val) + ggplot2::ylim(0, max_val)
ggplot2::ggsave(filename = paste0("../../../../out/figures/kevin/Writeup6s/Writeup6s_ppt-slides_lineage-scatterplot_day0-day10cocl2.png"),
                plot1, device = "png", width = 3, height = 3, units = "in", dpi=300)

plot1 <- ggplot2::ggplot(tab_mat,  ggplot2::aes(x=week5_COCL2, y=day10_COCL2)) + 
  ggplot2::geom_point(alpha = 0.3, shape = 16) + ggplot2::xlim(0, max_val) + ggplot2::ylim(0, max_val)
ggplot2::ggsave(filename = paste0("../../../../out/figures/kevin/Writeup6s/Writeup6s_ppt-slides_lineage-scatterplot_day10cocl2-week5cocl2.png"),
                plot1, device = "png", width = 3, height = 3, units = "in", dpi=300)

#####################

set.seed(10)
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
tab_mat <- tab_mat[,c("day10_COCL2", "day10_DABTRAM", "week5_COCL2", "week5_DABTRAM")]
stats::cor(log10(tab_mat+1))
n <- nrow(tab_mat); p <- ncol(tab_mat)
tab_mat <- tab_mat + matrix(stats::runif(n*p, max = 0.5), nrow = n, ncol = p)
tab_mat <- log10(tab_mat+1)
tab_mat <- apply(as.matrix.noquote(tab_mat),2,as.numeric)
tab_mat <- as.data.frame(tab_mat)

plot1 <- ggplot2::ggplot(tab_mat,  ggplot2::aes(x=day10_COCL2, y=day10_DABTRAM)) + 
  ggplot2::geom_point(alpha = 0.3, shape = 16)
ggplot2::ggsave(filename = paste0("../../../../out/figures/kevin/Writeup6s/Writeup6s_ppt-slides_lineage-scatterplot_day10cocl2-day10dab.png"),
                plot1, device = "png", width = 3, height = 3, units = "in", dpi=300)

plot1 <- ggplot2::ggplot(tab_mat,  ggplot2::aes(x=week5_COCL2, y=week5_DABTRAM)) + 
  ggplot2::geom_point(alpha = 0.3, shape = 16)
ggplot2::ggsave(filename = paste0("../../../../out/figures/kevin/Writeup6s/Writeup6s_ppt-slides_lineage-scatterplot_week5cocl2-week5dab.png"),
                plot1, device = "png", width = 3, height = 3, units = "in", dpi=300)

##########################

load(paste0("../../../../out/kevin/Writeup6r/Writeup6r_COCL2_day10_lineage-imputation_postprocess.RData"))

cell_imputed_score2 <- log10(exp(cell_imputed_score))
lineage_vec <- all_data$assigned_lineage[names(cell_imputed_score2)]

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
tab_vec <- table(lineage_vec)
tab_vec <- tab_vec[tab_vec >= 2]
week5_size <- sapply(names(tab_vec), function(lineage){
  tab_mat[lineage, "week5_COCL2"]
})
lineage_names <- names(sort(week5_size, decreasing = T))[1:8]
idx <- which(lineage_vec %in% lineage_names)

df <- data.frame(lineage = lineage_vec[idx],
                 imputed_count = cell_imputed_score2[idx])
df2 <- data.frame(lineage = "All",
                  imputed_count = cell_imputed_score2)
df <- rbind(df, df2)

col_vec <- c(rep("#999999", length(lineage_names)), "#E69F00")
names(col_vec) <- c(lineage_names, "All")

p1 <- ggplot2::ggplot(df, ggplot2::aes(x=lineage, y=imputed_count))
p1 <- p1 + ggplot2::geom_violin(trim=T, scale = "width", ggplot2::aes(fill=lineage))
p1 <- p1 + ggplot2::scale_fill_manual(values = col_vec) 
# p1 <- p1 + ggplot2::geom_jitter(shape=16, position=ggplot2::position_jitter(0.2), alpha = 0.3, size = 0.5)
p1 <- p1 + Seurat::NoLegend()
p1 <- p1 + ggplot2::geom_boxplot(width=0.05)
p1 <- p1 + ggplot2::scale_x_discrete(limits = c(lineage_names, "All"),
                                     guide = ggplot2::guide_axis(angle = 45))
p1 <- p1 + ggplot2::ylab("Week5 growth potential")
# p1 <- p1 + ggplot2::stat_summary(fun=mean, geom="point", shape=16, size=3, color="red")
# p1 <- p1 + ggplot2::stat_summary(fun=max, geom="point", shape=10, size=5, color="blue")
ggplot2::ggsave(filename = "../../../../out/figures/kevin/Writeup6s/Writeup6s_ppt-slides_lineage-growthpotential-violinplot.png",
                p1, device = "png", width = 6, height = 3, units = "in")

###################################

treatment <- "COCL2"
keep_vec <- rep(FALSE, ncol(all_data))
idx <- which(all_data$dataset %in% c("day0", paste0("day10_", treatment), paste0("week5_", treatment)))
keep_vec[idx] <- TRUE
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == TRUE)

set.seed(10)
all_data <- Seurat::RunUMAP(all_data, 
                            reduction = paste0("fasttopic_", treatment),
                            dims = 1:30)


cols.highlight_vec <- c("darkblue", "#DE2D26", "#00CC33")
names(cols.highlight_vec) <- c("week5", "day10", "day0")

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
lineage_vec <- order(tab_mat[,"week5_COCL2"], decreasing = T)[1:10]

pdf("../../../../out/figures/kevin/Writeup6s/Writeup6s_COCL2-lineages_week5.pdf", 
    onefile = T, width = 5, height = 5)
for(i in 1:length(lineage_vec)){
  print(paste0(i, " out of ", length(lineage_vec)))
  
  lineage_name <- rownames(tab_mat)[lineage_vec[i]]
  cell_idx <- intersect(which(all_data$dataset == paste0("week5_", treatment)),
                        which(all_data$assigned_lineage == lineage_name))
  cell_names <- colnames(all_data)[cell_idx]
  
  cell_idx2 <- intersect(which(all_data$dataset == paste0("day10_", treatment)),
                         which(all_data$assigned_lineage == lineage_name))
  cell_names2 <- colnames(all_data)[cell_idx2]
  
  cell_idx3 <- intersect(which(all_data$dataset == "day0"),
                         which(all_data$assigned_lineage == lineage_name))
  cell_names3 <- colnames(all_data)[cell_idx3]
  
  p1 <- Seurat::DimPlot(all_data, reduction = "umap",
                        cells.highlight = list(week5 = cell_names, day10 = cell_names2, day0 = cell_names3),
                        cols.highlight = cols.highlight_vec)
  p1 <- p1 + ggplot2::ggtitle(paste0(treatment, ": ", lineage_name, " (", length(cell_names), " week5 cells,\nfrom ",
                                     length(cell_idx2), " day10 cells and ", length(cell_idx3), " day0 cells)")) +
    ggplot2::theme(plot.title=ggplot2::element_text(size=9)) + Seurat::NoLegend()
  print(p1)
}

dev.off()