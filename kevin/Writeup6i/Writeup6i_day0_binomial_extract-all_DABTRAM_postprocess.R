rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6i/day0_binomial_extract-DABTRAM_tmp.RData")
source("fisher_exact_test.R")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

tmp <- sapply(result_list, function(x){all(!is.null(x))})
result_list <- result_list[which(tmp)]

pvalue_func <- function(mat){
  bin_win <- mat[1,]
  bin_die <- mat[2,]
  
  if(all(bin_win[c("Mid", "Close")] == 0) | all(bin_die[c("Mid", "Close")] == 0)){
    val1 <- NA
  } else {
    count_group1 <- bin_win[c("Close", "Mid")]
    count_group2 <- bin_die[c("Close", "Mid")]
    res1 <- exact_fisher_pvalue(count_group1 = count_group1,
                                count_group2 = count_group2)
    val1 <- res1$log10prob_val
  }
  
  if(all(bin_win[c("Far", "Close")] == 0) | all(bin_die[c("Far", "Close")] == 0)){
    val2 <- NA
  } else {
    count_group1 <- bin_win[c("Close", "Far")]
    count_group2 <- bin_die[c("Close", "Far")]
    res2 <- exact_fisher_pvalue(count_group1 = count_group1,
                                count_group2 = count_group2)
    val2 <- res2$log10prob_val
  }
  
  c(val1, val2)
}

len <- length(result_list)
pvalue_mat <- sapply(1:len, function(i){
  if(i %% floor(len/10) == 0) cat('*')
  if(length(result_list[[i]]) == 0) return(rep(NA,2))
  pvalue_func(result_list[[i]])
})
colnames(pvalue_mat) <- names(result_list)
rownames(pvalue_mat) <- c("Mid", "Far")

vec1 <- pvalue_mat[1,]
vec1 <- vec1[!is.na(vec1)]
vec1 <- 10^(-vec1)

vec2 <- pvalue_mat[2,]
vec2 <- vec2[!is.na(vec2)]
vec2 <- 10^(-vec2)

png(paste0("../../../../out/figures/Writeup6i/Writeup6i_day0_binomial-pvalue_DABTRAM.png"),
    height = 1200, width = 2500, res = 300, units = "px")
par(mfrow = c(1,2))
hist(vec1, breaks = 15, xlab = "P-value", main = "Mid vs. Close (DABTRAM)")

hist(vec2, breaks = 15, xlab = "P-value", main = "Far vs. Close (DABTRAM)")
graphics.off()

vec1b <- stats::p.adjust(vec1, method = "BH")
vec1b[vec1b <= .05]

vec2b <- stats::p.adjust(vec2, method = "BH")
vec2b[vec2b <= .05]

###############################

entropy_func <- function(vec){
  if(all(vec == 0)) return(entropy_func(c(1,1)))
  vec <- vec/sum(vec)
  log_vec <- log(vec)
  log_vec[vec == 0] <- 0
  -sum(vec*log_vec)
}

len <- length(result_list)
mid_entropy_mat <- sapply(1:len, function(i){
  if(i %% floor(len/10) == 0) cat('*')
  if(length(result_list[[i]]) == 0) return(rep(NA,2))
  
  mat <- result_list[[i]]
  bin_win <- mat[1,]
  bin_die <- mat[2,]
  
  count_group1 <- bin_win[c("Close", "Mid")]
  count_group2 <- bin_die[c("Close", "Mid")]
  
  vec <- c(entropy_func(count_group1),
           entropy_func(count_group2))
  names(vec) <- c("Win", "Lose")
  vec
})
colnames(mid_entropy_mat) <- names(result_list)

far_entropy_mat <- sapply(1:len, function(i){
  if(i %% floor(len/10) == 0) cat('*')
  if(length(result_list[[i]]) == 0) return(rep(NA,2))
  
  mat <- result_list[[i]]
  bin_win <- mat[1,]
  bin_die <- mat[2,]
  
  count_group1 <- bin_win[c("Close", "Far")]
  count_group2 <- bin_die[c("Close", "Far")]
  
  vec <- c(entropy_func(count_group1),
           entropy_func(count_group2))
  names(vec) <- c("Win", "Lose")
  vec
})
colnames(far_entropy_mat) <- names(result_list)

##################################

source("../Writeup6b/gene_list.R")
labeling_vec <- rep("1", ncol(mid_entropy_mat))
names(labeling_vec) <- colnames(mid_entropy_mat)
labeling_vec[intersect(
  names(keygene_vec),
  sort(unlist(keygenes))
)] <- "2"
labeling_vec[names(vec1b)[which(vec1b <= .05)]] <- "3"

df <- data.frame(diff_entropy = mid_entropy_mat[1,]-mid_entropy_mat[2,],
                 log10pval = pvalue_mat[1,],
                 name = colnames(mid_entropy_mat),
                 labeling = labeling_vec)
df <- df[which(!is.na(df$log10pval)),]
# put all the labeling == TRUE on bottom
df <- df[c(which(df[,"labeling"] == "1"), which(df[,"labeling"] %in% c("2", "3"))),]

p1 <- ggplot2::ggplot(df, ggplot2::aes(x = diff_entropy, 
                                       y = log10pval))
p1 <- p1 + ggplot2::geom_point(ggplot2::aes(color = labeling))
p1 <- p1 + ggplot2::scale_colour_manual(values=c("black", "darkgoldenrod1", "red"))
p1 <- p1 + ggrepel::geom_text_repel(data = subset(df, labeling %in% c("2", "3")),
                                    ggplot2::aes(label = name, color = labeling),
                                    box.padding = ggplot2::unit(0.5, 'lines'),
                                    point.padding = ggplot2::unit(1.6, 'lines'),
                                    max.overlaps = 50)
p1 <- p1 + ggplot2::ggtitle(paste0("DABTRAM Day0 (Mid vs close)")) +
  ggplot2::xlab("Entropy difference: (Win-Lose)") + ggplot2::ylab("Fisher exact pvalue (Log10)")
p1 <- p1 + Seurat::NoLegend()

ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6i/Writeup6i_day0_binomial_DABTRAM_volcano_mid-vs-close.png"),
                p1, device = "png", width = 10, height = 10, units = "in")

###############

labeling_vec <- rep("1", ncol(mid_entropy_mat))
names(labeling_vec) <- colnames(mid_entropy_mat)
labeling_vec[intersect(
  names(keygene_vec),
  sort(unlist(keygenes))
)] <- "2"
labeling_vec[names(vec2b)[which(vec2b <= .05)]] <- "3"

pvalue_vec <- pvalue_mat[2,]
idx <- which(is.infinite(pvalue_vec))
pvalue_vec[idx] <- max(pvalue_vec[-idx], na.rm=T)

df <- data.frame(diff_entropy = far_entropy_mat[1,]-far_entropy_mat[2,],
                 log10pval = pvalue_vec,
                 name = colnames(mid_entropy_mat),
                 labeling = labeling_vec)
df <- df[which(!is.na(df$log10pval)),]
# put all the labeling == TRUE on bottom
df <- df[c(which(df[,"labeling"] == "1"), which(df[,"labeling"] %in% c("2", "3"))),]

p1 <- ggplot2::ggplot(df, ggplot2::aes(x = diff_entropy, 
                                       y = log10pval))
p1 <- p1 + ggplot2::geom_point(ggplot2::aes(color = labeling))
p1 <- p1 + ggplot2::scale_colour_manual(values=c("black", "darkgoldenrod1", "red"))
p1 <- p1 + ggrepel::geom_text_repel(data = subset(df, labeling %in% c("2", "3")),
                                    ggplot2::aes(label = name, color = labeling),
                                    box.padding = ggplot2::unit(0.5, 'lines'),
                                    point.padding = ggplot2::unit(1.6, 'lines'),
                                    max.overlaps = 50)
p1 <- p1 + ggplot2::ggtitle(paste0("DABTRAM Day0 (Far vs close)")) +
  ggplot2::xlab("Entropy difference: (Win-Lose)") + ggplot2::ylab("Fisher exact pvalue (Log10)")
p1 <- p1 + Seurat::NoLegend()

ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6i/Writeup6i_day0_binomial_DABTRAM_volcano_far-vs-close.png"),
                p1, device = "png", width = 10, height = 10, units = "in")

##########


load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")

treatment <- "DABTRAM"

Seurat::DefaultAssay(all_data) <- "RNA"
all_data[["geneActivity"]] <- NULL
all_data[["fasttopic_CIS"]] <- NULL
all_data[["fasttopic_COCL2"]] <- NULL
all_data[["common_tcca"]] <- NULL
all_data[["distinct1_tcca"]] <- NULL
all_data[["distinct2_tcca"]] <- NULL
all_data[["activity.umap"]] <- NULL

keep_vec <- rep(NA, ncol(all_data))
idx1 <- which(all_data$dataset %in% c("day0", 
                                      paste0("day10_", treatment), 
                                      paste0("week5_", treatment)))
idx2 <- which(all_data$assigned_posterior >= 0.5)
idx3 <- which(!is.na(all_data$assigned_lineage))
keep_vec[intersect(intersect(idx1, idx2), idx3)] <- TRUE
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == TRUE)

gene_vec_all <- sort(intersect(gene_vec_all,
                               rownames(all_data[["Saver"]]@scale.data)))

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)

surviving_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("day10_", treatment)] >= 20)]
dying_lineages <- rownames(tab_mat)[intersect(
  which(tab_mat[,paste0("day10_", treatment)] <= 2),
  which(tab_mat[,paste0("week5_", treatment)] <= 2)
)]
winning_idx <- intersect(
  which(all_data$assigned_lineage %in% surviving_lineages),
  which(all_data$dataset == "day0")
)
dying_idx <- intersect(
  which(all_data$assigned_lineage %in% dying_lineages),
  which(all_data$dataset == "day0")
)
winning_cells <- colnames(all_data)[winning_idx]
dying_cells <- colnames(all_data)[dying_idx]

keep_vec <- rep(NA, ncol(all_data))
keep_vec[winning_idx] <- paste0("day0_winner_", treatment)
keep_vec[dying_idx] <- paste0("day0_loser_", treatment)
keep_vec <- as.factor(keep_vec)
table(keep_vec)
all_data2 <- all_data
all_data2$keep <- keep_vec
all_data2 <- subset(all_data2, keep %in% c(paste0("day0_winner_", treatment),
                                           paste0("day0_loser_", treatment)))

Seurat::DefaultAssay(all_data2) <- "ATAC"
Seurat::Idents(all_data2) <- "keep"

gene_vec <- sort(c("ACBD3", "NEIL3", "IPCEF1", "EMC2", "NRCAM", 
                   "RNF144A", "LINC00973", "LAMTOR2"))

pdf(paste0("../../../../out/figures/Writeup6i/Writeup6i_binomial_coverage_", treatment, "_select-genes.pdf"),
    onefile = T, width = 9, height = 4.5)
for(gene in gene_vec){
  plot1 <- Signac::CoveragePlot(
    object = all_data2,
    region = gene,
    features = gene,
    extend.upstream = 5000,
    extend.downstream = 5000
  )
  
  print(plot1)
}

dev.off() 
