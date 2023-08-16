# THIS VERSION IS OUTDATED. PLEASE USE Writeup6k_COCL2_day10_DE_stepdown_at-day0.R

rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6j/Writeup6j_COCL2_day10_lineage-imputation_stepdown_step2.RData")
load("../../../../out/kevin/Writeup6i/Writeup6i_COCL2-day10_extracted.RData")
all_data$tier_vec <- all_data$keep

# extract the average test and train
len_vec <- sapply(coefficient_list_list, function(lis){
  length(lis$var_next_iteration)
})
len <- length(coefficient_list_list)
train_vec <- sapply(coefficient_list_list, function(x){
  x$fit$objective_val
})
test_vec <- colMeans(loocv_mat)

idx <- which.min(test_vec)
lineage_res <- coefficient_list_list[[idx]]$fit
var_names <- names(lineage_res$coefficient_vec)

fasttopic_mat <- all_data[["fasttopic_COCL2"]]@cell.embeddings
lsi_mat <- all_data[["lsi"]]@cell.embeddings

cell_features <- cbind(1, 
                       scale(fasttopic_mat), 
                       scale(lsi_mat))
colnames(cell_features)[1] <- "Intercept"
cell_features <- cell_features[,var_names]
cell_lineage <- all_data$assigned_lineage
cell_lineage <- cell_lineage[!is.na(cell_lineage)]
cell_features <- cell_features[names(cell_lineage),]
uniq_lineage <- sort(unique(cell_lineage))
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
lineage_current_count <- tab_mat[uniq_lineage,"day10_COCL2"]
lineage_future_count <- tab_mat[uniq_lineage,"week5_COCL2"]

cell_imputed_count <- as.numeric(exp(cell_features %*% lineage_res$coefficient_vec))
names(cell_imputed_count) <- rownames(cell_features)
lineage_imputed_count <- sapply(uniq_lineage, function(lineage){
  sum(cell_imputed_count[which(cell_lineage == lineage)])
})

imputed_vec <- rep(NA, ncol(all_data))
names(imputed_vec) <- colnames(all_data)
imputed_vec[names(cell_imputed_count)] <- log10(cell_imputed_count)
all_data$imputed_count <- imputed_vec

############################

lineage_score <- sapply(1:nrow(tab_mat), function(i){
  if(i %% floor(nrow(tab_mat)/10) == 0) cat('*')
  lineage_name <- rownames(tab_mat)[i]
  cell_names <- colnames(all_data)[which(all_data$assigned_lineage == lineage_name)]
  if(length(cell_names) == 0) return(NA)
  mean(all_data$imputed_count[cell_names], na.rm = T)
})
names(lineage_score) <- rownames(tab_mat)
quantile(lineage_score, na.rm = T, probs = seq(0,1,length.out=11))
length(which(lineage_score > 0))

lineage_high <- names(lineage_score)[which(lineage_score >= stats::quantile(lineage_score, probs = 0.75, na.rm = T))]
lineage_low <- names(lineage_score)[which(lineage_score <= stats::quantile(lineage_score, probs = 0.25, na.rm = T))]

level_vec <- rep(NA, ncol(all_data))
level_vec[intersect(
  which(all_data$assigned_lineage %in% lineage_high),
  which(all_data$dataset == "day0")
)] <- "high"
level_vec[intersect(
  which(all_data$assigned_lineage %in% lineage_low),
  which(all_data$dataset == "day0")
)] <- "low"
table(level_vec)
all_data$level <- level_vec
Seurat::Idents(all_data) <- "level"

pvalue_list <- lapply(all_data[["Saver"]]@var.features, function(gene){
  x_vec <- all_data[["Saver"]]@scale.data[gene, which(all_data$level == "high")]
  y_vec <- all_data[["Saver"]]@scale.data[gene, which(all_data$level == "low")]
  
  if(diff(range(x_vec)) <= 1e-6 || diff(range(y_vec)) <= 1e-6 ) return(NA)
  
  test_res <- stats::wilcox.test(x = x_vec, y = y_vec)
  list(stat = mean(x_vec) - mean(y_vec), pvalue = test_res$p.value)
})
names(pvalue_list) <- all_data[["Saver"]]@var.features
pvalue_list <- pvalue_list[sapply(pvalue_list, function(x){!all(is.na(x))})]

pvalue_vec <- sapply(pvalue_list, function(x){x$pvalue})
pvalue_vec <- stats::p.adjust(pvalue_vec, method = "BH")
quantile(pvalue_vec)
length(which(pvalue_vec <= 0.05))

quantile(sapply(pvalue_list, function(x){x$stat}))
quantile(sapply(pvalue_list, function(x){-log10(x$pvalue)}))

##############

source("../Writeup6b/gene_list.R")
keygene_vec <- sort(unique(c(keygenes$jackpot, keygenes$COCL2)))
keygene_vec <- sort(intersect(keygene_vec, names(pvalue_list)))

diff_vec <- sapply(pvalue_list, function(x){x$stat})
logpval_vec <- sapply(pvalue_list, function(x){-log10(x$pvalue)})

labeling_vec <- rep("1", length(diff_vec))
names(labeling_vec) <- names(diff_vec)
labeling_vec[names(logpval_vec)[which(logpval_vec >= 2)]] <- "3"
labeling_vec[intersect(keygene_vec, names(labeling_vec))] <- "2"

df <- data.frame(difference = diff_vec,
                 log10pval = logpval_vec,
                 name = names(logpval_vec),
                 labeling = labeling_vec)
# put all the labeling == TRUE on bottom
df <- df[c(which(df[,"labeling"] == "1"), which(df[,"labeling"] %in% c("2", "3"))),]

p1 <- ggplot2::ggplot(df, ggplot2::aes(x = difference, 
                                       y = log10pval))
p1 <- p1 + ggplot2::geom_point(ggplot2::aes(color = labeling))
p1 <- p1 + ggplot2::scale_colour_manual(values=c("black", "red", "blue"))
p1 <- p1 + ggrepel::geom_text_repel(data = subset(df, labeling %in% c("2", "3")),
                                    ggplot2::aes(label = name, color = labeling),
                                    box.padding = ggplot2::unit(0.5, 'lines'),
                                    point.padding = ggplot2::unit(1.6, 'lines'),
                                    max.overlaps = 50)
p1 <- p1 + ggplot2::ggtitle(paste0("COCL2 Day10 DE based on growth potential\n(Top 25% to Bottom 25%)")) +
  ggplot2::xlab("Mean difference in SAVER gene exp.: (Best-Worst)") + ggplot2::ylab("T test p-value (-Log10)")
p1 <- p1 + Seurat::NoLegend()

ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6k/Writeup6k_COCL2-day10_stepdown_DE_lineage-at-day0_volcano.png"),
                p1, device = "png", width = 10, height = 10, units = "in")

###################

# make a violin plot
gene_vec <- sort(names(logpval_vec)[which(logpval_vec >= 2)])
cell_idx <- which(all_data$level %in% c("high", "low"))

pdf("../../../../out/figures/Writeup6k/Writeup6k_COCL2-day10_stepdown_DE_lineage-at-day0_violin_highest-DE.pdf", 
    onefile = T, width = 5, height = 5)

for(gene in gene_vec){
  
  df <- data.frame(rna = all_data[["Saver"]]@scale.data[gene,cell_idx],
                   tier = all_data$level[cell_idx])
  col_vec <- c("lightgray", "dodgerblue4")
  names(col_vec) <- c("low", "high")
  
  p1 <- ggplot2::ggplot(df, ggplot2::aes(x=tier, y=rna)) +
    ggplot2::geom_violin(trim=FALSE, ggplot2::aes(fill=tier)) +
    ggplot2::scale_fill_manual(values=col_vec) +
    ggplot2::geom_jitter(shape=16, position=ggplot2::position_jitter(0.2), alpha = 0.3, size = 0.5) + 
    Seurat::NoLegend() + 
    ggplot2::geom_boxplot(width=0.05) +
    ggplot2::ggtitle(paste0("Saver: ", gene, ": -Log10pvalue=", round(-log10(pvalue_list[[gene]]$pvalue), 2)))
  print(p1)
}

dev.off()
