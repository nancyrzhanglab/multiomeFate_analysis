rm(list=ls())
# from https://rpubs.com/KDVdecisions/triadtutorial1
library(ggplot2)
library(ggtern)

load("../../out/Writeup6n/Writeup6n_COCL2_topic-simplex.RData")

# adjust just so points don't lie on the boundary
cocl2_d10 <- cocl2_d10 + 0.1
for(i in 1:nrow(cocl2_d10)){
  cocl2_d10[i,] <- cocl2_d10[i,]/sum(cocl2_d10[i,])
}

cocl2_d10_v2 <- cocl2_d10*100
cocl2_d10_v2 <- as.data.frame(cocl2_d10_v2)
colnames(cocl2_d10_v2) <- sapply(colnames(cocl2_d10_v2),
                                 function(x){
                                   paste0("T", strsplit(x, split = "_")[[1]][2])
                                 })
week5_size <- tab_mat[lineage_names,"week5_COCL2"]
lineage_names <- names(week5_size)[order(week5_size, decreasing = T)]

pdf("../../out/figures/Writeup6n/Writeup6n_COCL2_topic-simplex.pdf", 
    onefile = T, width = 8, height = 8)

for(lineage_name in lineage_names){
  cell_idx <- names(lineage_vec)[which(lineage_vec == lineage_name)]
  cocl2_d10_v3 <- cocl2_d10_v2[cell_idx,]
  
  p1 <- ggtern::ggtern(cocl2_d10_v3, ggplot2::aes(x = T13,
                                                  y = T28,
                                                  z = T9))
  p1 <- p1 + ggtern::stat_density_tern(aes(fill=..level.., alpha=..level..), geom='polygon') 
  p1 <- p1 + ggplot2::scale_fill_gradient2(high = "blue")
  p1 <- p1 + ggplot2::geom_point()
  p1 <- p1 + ggtern::theme_showarrows() 
  p1 <- p1 + ggplot2::ggtitle(paste0("COCL2: D10 size = ", tab_mat[lineage_name, "day10_COCL2"],
                            ", W5 size = ", tab_mat[lineage_name, "week5_COCL2"]))
  p1 <- p1 + ggplot2::guides(color = "none", fill = "none", alpha = "none")
  print(p1)
}

dev.off()