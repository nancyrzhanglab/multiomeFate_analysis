rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6j/Writeup6j_DABTRAM_day0_lineage-imputation_step-up_step2.RData")
load("../../../../out/kevin/Writeup6j/Writeup6j_DABTRAM-day0_extracted.RData")
all_data2$tier_vec <- all_data2$keep
load("../../../../out/kevin/Writeup6j/Writeup6j_DABTRAM_day0_lineage-imputation_stepwise-up_training.RData")

# extract the average test and train
len_vec <- sapply(coefficient_list_list, function(lis){
  length(lis$var_next_iteration)
})
len <- length(coefficient_list_list)
train_vec <- sapply(training_list, function(x){
  x$fit$objective_val
})
test_vec <- colMeans(loocv_mat)

idx <- which.min(test_vec)
lineage_res <- training_list[[idx]]$fit
var_names <- names(lineage_res$coefficient_vec)

fasttopic_mat <- all_data2[["fasttopic_DABTRAM"]]@cell.embeddings
lsi_mat <- all_data2[["lsi"]]@cell.embeddings

cell_features <- cbind(1, 
                       scale(fasttopic_mat), 
                       scale(lsi_mat))
colnames(cell_features)[1] <- "Intercept"
cell_features <- cell_features[,var_names]
cell_lineage <- all_data2$assigned_lineage
cell_lineage <- cell_lineage[!is.na(cell_lineage)]
cell_features <- cell_features[names(cell_lineage),]
uniq_lineage <- sort(unique(cell_lineage))
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)

cell_imputed_count <- as.numeric(exp(cell_features %*% lineage_res$coefficient_vec))
names(cell_imputed_count) <- rownames(cell_features)
lineage_imputed_count <- sapply(uniq_lineage, function(lineage){
  sum(cell_imputed_count[which(cell_lineage == lineage)])
})

imputed_vec <- rep(NA, ncol(all_data))
names(imputed_vec) <- colnames(all_data)
imputed_vec[names(cell_imputed_count)] <- log10(cell_imputed_count)
all_data$imputed_count <- imputed_vec

#################

tab_mat <- tab_mat[order(tab_mat[,"day10_DABTRAM"], decreasing = T),]
lineage_idx <- intersect(which(tab_mat[,"day0"] >= 3),
                         which(tab_mat[,"day10_DABTRAM"] >= 10))
large_lineages <- rownames(tab_mat)[lineage_idx]
all_data3 <- subset(all_data, 
                    dataset == "day0" & assigned_lineage %in% large_lineages)

large_lineages2 <- large_lineages
lineage_vec <- all_data3$assigned_lineage
# replace the lineage names
for(lineage in large_lineages){
  idx <- which(lineage_vec == lineage)
  lineage_name <- paste0(lineage, "_D10:", tab_mat[lineage,"day10_DABTRAM"])
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
p1 <- p1 + ggplot2::ggtitle(paste0("ANOVA -Log10(pvalue)=", round(-log10(anova_res$p.value), 2)))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6k/Writeup6k_DABTRAM-day0_imputation_stepup-violinplot.png"),
                p1, device = "png", width = 13, height = 6, units = "in")


