rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6k/Writeup6k_DABTRAM_day10_lineage-imputation_stepdown_postprocessed.RData")

############################

# do the chromvar analysis 
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)

# let's define winners and losers
# score each lineage by its mean day10 score
lineage_score <- sapply(1:nrow(tab_mat), function(i){
  if(i %% floor(nrow(tab_mat)/10) == 0) cat('*')
  lineage_name <- rownames(tab_mat)[i]
  cell_names <- colnames(all_data)[which(all_data$assigned_lineage == lineage_name)]
  if(length(cell_names) == 0) return(NA)
  # quantile(all_data$imputed_count[cell_names], probs = 0.75, na.rm = T)
  length(which(all_data$imputed_count[cell_names] >= 0))
})
names(lineage_score) <- rownames(tab_mat)
quantile(lineage_score, na.rm = T, probs = seq(0,1,length.out=11))
length(which(lineage_score > 0))

lineage_score2 <- sapply(1:nrow(tab_mat), function(i){
  if(i %% floor(nrow(tab_mat)/10) == 0) cat('*')
  lineage_name <- rownames(tab_mat)[i]
  cell_names <- colnames(all_data)[which(all_data$assigned_lineage == lineage_name)]
  if(length(cell_names) == 0) return(NA)
  mean(all_data$imputed_count[cell_names], na.rm = T)
})
names(lineage_score2) <- rownames(tab_mat)
quantile(lineage_score2, na.rm = T, probs = seq(0,1,length.out=11))
length(which(lineage_score2 > 0))

cell_winning_idx <- intersect(which(all_data$assigned_lineage %in% names(lineage_score2)[intersect(
  which(lineage_score > 0),
  which(lineage_score2 >= -1)
)]),
which(all_data$dataset == "day0"))
cell_winning_names <- colnames(all_data)[cell_winning_idx]
cell_losing_idx <- intersect(which(all_data$assigned_lineage %in% names(lineage_score2)[which(lineage_score2 <= -1)]),
                             which(all_data$dataset == "day0"))
cell_losing_names <- colnames(all_data)[cell_losing_idx]
length(cell_winning_names); length(cell_losing_names)

############################

level_vec <- rep(NA, ncol(all_data))
names(level_vec) <- colnames(all_data)
level_vec[cell_winning_names] <- "high"
level_vec[cell_losing_names] <- "low"
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
quantile(sapply(pvalue_list, function(x){x$stat}))
quantile(sapply(pvalue_list, function(x){-log10(x$pvalue)}))

pvalue_vec <- sapply(pvalue_list, function(x){x$pvalue})
names(pvalue_vec) <- names(pvalue_list)
pvalue_vec <- stats::p.adjust(pvalue_vec, method = "BH")

# make a volcano plot now
source("../Writeup6b/gene_list.R")
keygene_vec <- sort(unique(c(keygenes$jackpot, keygenes$DABTRAM)))
keygene_vec <- sort(intersect(keygene_vec, names(pvalue_list)))

diff_vec <- sapply(pvalue_list, function(x){x$stat})
logpval_vec <- sapply(pvalue_list, function(x){-log10(x$pvalue)})

labeling_vec <- rep("1", length(diff_vec))
names(labeling_vec) <- names(diff_vec)
labeling_vec[names(logpval_vec)[which(logpval_vec >= 5)]] <- "3"
labeling_vec[intersect(keygene_vec, names(labeling_vec))] <- "2"
table(labeling_vec)

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
p1 <- p1 + ggplot2::geom_hline(yintercept=min(logpval_vec[which(pvalue_vec <= 0.05)]), linetype="dashed", 
                               color = "red", linewidth=2)
p1 <- p1 + ggplot2::ggtitle(paste0("DABTRAM Day0 DE based on lineage's mean growth potential at Day10\n(Top 25% to Bottom 25%)")) +
  ggplot2::xlab("Mean difference in SAVER gene exp.: (Best-Worst)") + ggplot2::ylab("T test p-value (-Log10)")
p1 <- p1 + Seurat::NoLegend()

ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6k/Writeup6k_DABTRAM-day10_stepdown_at-day0_DE_volcano.png"),
                p1, device = "png", width = 10, height = 10, units = "in")

###################

# make a violin plot
gene_vec <- sort(names(logpval_vec)[which(logpval_vec >= 300)])
cell_idx <- which(all_data$level %in% c("high", "low"))

pdf("../../../../out/figures/Writeup6k/Writeup6k_DABTRAM-day10_stepdown_DE_violin_highest-DE.pdf", 
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
