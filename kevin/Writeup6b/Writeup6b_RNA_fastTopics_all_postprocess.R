rm(list=ls())
library(Seurat)
library(Signac)
library(fastTopics)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")
load("../../../../out/kevin/Writeup6b/Writeup6b_RNA_fasttopics_all.RData")

ft_all <- Seurat::CreateDimReducObject(embeddings = topic_res$L,
                                       loadings = topic_res$F,
                                       key = "fastTopicAll_")
all_data[["fasttopic_all"]] <- ft_all

for(i in 1:6){
  print(i)
  
  Seurat::Idents(all_data) <- "dataset"
  plot1 <- Seurat::VlnPlot(all_data,
                           features = paste0("fastTopicAll_", (((i-1)*9)+1):min(i*9, 50)),
                           ncol = 3,
                           pt.size = 0)
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6b/Writeup6b_fasttopics_all_", i, ".png"),
                  plot1, device = "png", width = 15, height = 15, units = "in")
}

cor_mat <- cor(sqrt(topic_res$F))


png(paste0("../../../../out/figures/Writeup6b/Writeup6b_RNA_ftTopic-all_correlation.png"),
    width = 2000, height = 2000, units = "px", res = 300)
par(mar=c(5,5,1,1))
fields::image.plot(cor_mat)
graphics.off()

#######################################

mat <- topic_res$F
idx <- which(rownames(mat) %in% sort(unique(unlist(keygenes))))
mat <- mat[idx,]
rowname_vec <- rownames(mat)
l1_vec <- apply(mat, 1, sum)
mat <- diag(1/l1_vec) %*% mat
rownames(mat) <- rowname_vec

hclust_res <- stats::hclust(stats::dist(mat))
mat <- mat[hclust_res$order,]

col_vec <- viridis::viridis(25)
break_vec <- seq(0, 1, length.out = 26)
ratio <- min(2000*(nrow(mat)-1)/(ncol(mat)-1), 6500)/(2000)
png(paste0("../../../../out/figures/Writeup6b/Writeup6b_RNA_ftTopic-all_keygenes.png"),
    height = 2500*ratio, width = 2500, res = 500, units = "px")
par(mar = c(0.1, 0.1, 0.1, 0.1))
image(tiltedCCA:::.rotate(mat), 
      asp = ratio, ylim = c(0,1),
      main = "", xlab = "", ylab = "",
      xaxt = "n", yaxt = "n", bty = "n",
      breaks = break_vec, col = col_vec)

# draw lines
x_indent_val <- 1/(2*(ncol(mat)-1))
x_vec <- seq(-x_indent_val, 1+x_indent_val, by = 2*x_indent_val)
for(x in x_vec[seq(51,6,by=-5)[-1]]){
  graphics::lines(y = c(0,1), x = rep(x, 2), lwd = 2, col = "white")
  graphics::lines(y = c(0,1), x = rep(x, 2), lwd = 1.5, lty = 2)
}

# label genes
y_indent_val <- 1/(2*(nrow(mat)-1))
y_vec <- seq(-y_indent_val, 1+y_indent_val, by = 2*y_indent_val)
for(y in y_vec[seq(6,length(y_vec), by = 5)]){
  graphics::lines(x = c(0,1), y = rep(y, 2), lwd = 2, col = "white")
}

x_vec_label <- rep(x_vec[c(3,5)], times = ceiling(nrow(mat)/2))[1:nrow(mat)]
y_vec_label <- seq(1, 0, by = -2*y_indent_val)
graphics::text(x = x_vec_label, y = y_vec_label, labels = rownames(mat),
               col = "white", cex = 0.5)

graphics.off()
