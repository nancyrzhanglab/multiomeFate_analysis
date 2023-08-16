rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)


load("../../../../out/kevin/Writeup6j/Writeup6j_COCL2-day0_extracted.RData")
load("../../../../out/kevin/Writeup6j/Writeup6j_COCL2_day0_lineage-imputation_stepdown-LOOCV_concise-postprocessed.RData")

tmp <- rep(NA, ncol(all_data))
names(tmp) <- colnames(all_data)
cell_names <- intersect(colnames(all_data), names(cell_imputed_count))
tmp[cell_names] <- cell_imputed_count[cell_names]
all_data$imputed_count <- tmp

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

########################

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
tab_mat <- tab_mat[order(tab_mat[,"day10_COCL2"], decreasing = T),]
lineage_idx <- intersect(which(tab_mat[,"day0"] >= 3),
                         which(tab_mat[,"day10_COCL2"] >= 20))
large_lineages <- rownames(tab_mat)[lineage_idx]
length(large_lineages)

all_data3 <- subset(all_data, 
                    dataset == "day0" & assigned_lineage %in% large_lineages)

large_lineages2 <- large_lineages
lineage_vec <- all_data3$assigned_lineage
# replace the lineage names
for(lineage in large_lineages){
  idx <- which(lineage_vec == lineage)
  lineage_name <- paste0(lineage, "_D10:", tab_mat[lineage,"day10_COCL2"])
  lineage_vec[idx] <- lineage_name
  
  idx <- which(large_lineages == lineage)
  large_lineages2[idx] <- lineage_name
}

###############

df <- data.frame(lineage = lineage_vec,
                 imputed_count = all_data3$imputed_count)
df$lineage <- as.factor(df$lineage)

set.seed(10)
anova_res <- stats::oneway.test(imputed_count ~ lineage, data = df)

total_std <- sum((df$imputed_count - mean(df$imputed_count))^2)
within_lineage_std <- sum(sapply(levels(df$lineage), function(lineage_name){
  idx <- which(df$lineage == lineage_name)
  sum((df$imputed_count[idx] - mean(df$imputed_count[idx]))^2)
}))
across_lineage_std <- sum(sapply(levels(df$lineage), function(lineage_name){
  idx <- which(df$lineage == lineage_name)
  mean_val <- mean(df$imputed_count[idx])
  length(idx) * (mean_val - mean(df$imputed_count))^2 
}))
across_lineage_std + within_lineage_std - total_std
lineage_effect <- round(across_lineage_std/total_std*100,1)

###############

df <- data.frame(lineage = lineage_vec,
                 imputed_count = all_data3$imputed_count)
df2 <- data.frame(lineage = "All",
                  imputed_count = all_data$imputed_count[!is.na(all_data$imputed_count)])
df <- rbind(df, df2)

col_vec <- c(rep("#999999", length(large_lineages2)), "#E69F00")
names(col_vec) <- c(large_lineages2, "All")

p1 <- ggplot2::ggplot(df, ggplot2::aes(x=lineage, y=imputed_count))
p1 <- p1 + ggplot2::geom_violin(trim=T, scale = "width", ggplot2::aes(fill=lineage))
p1 <- p1 + ggplot2::scale_fill_manual(values = col_vec) 
p1 <- p1 + ggplot2::geom_jitter(shape=16, position=ggplot2::position_jitter(0.2), alpha = 0.3, size = 0.5)
p1 <- p1 + Seurat::NoLegend()
p1 <- p1 + ggplot2::geom_boxplot(width=0.05)
p1 <- p1 + ggplot2::scale_x_discrete(limits = c(large_lineages2, "All"),
                                     guide = ggplot2::guide_axis(angle = 45))
p1 <- p1 + ggplot2::ylab("Log10(Day10 prediction)")
p1 <- p1 + ggplot2::stat_summary(fun=mean, geom="point", shape=16, size=3, color="red")
p1 <- p1 + ggplot2::stat_summary(fun=max, geom="point", shape=16, size=5, color="blue")
p1 <- p1 + ggplot2::ggtitle(paste0("ANOVA -Log10(pvalue)=", round(-log10(anova_res$p.value), 2), ", Lineage effect = ", lineage_effect, "%"))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6k/Writeup6k_COCL2-day0_imputation_stepdown-LOOCV-violinplot.png"),
                p1, device = "png", width = 13, height = 6, units = "in")


