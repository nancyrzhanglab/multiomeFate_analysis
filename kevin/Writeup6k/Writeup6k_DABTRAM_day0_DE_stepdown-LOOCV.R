rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6k/Writeup6k_DABTRAM_day0_lineage-imputation_stepdown-LOOCV_postprocessed.RData")

############################

keep_vec <- !is.na(all_data$imputed_count)
all_data$keep <- keep_vec
all_data2 <- subset(all_data, keep == TRUE)
quantile(all_data2$imputed_count)

level_vec <- rep(NA, ncol(all_data2))
level_vec[which(all_data2$imputed_count >= quantile(all_data2$imputed_count, probs = 0.75))] <- "high"
level_vec[which(all_data2$imputed_count <= quantile(all_data2$imputed_count, probs = 0.25))] <- "low"
all_data2$level <- level_vec
Seurat::Idents(all_data2) <- "level"

pvalue_list <- lapply(all_data2[["Saver"]]@var.features, function(gene){
  x_vec <- all_data2[["Saver"]]@scale.data[gene, which(all_data2$level == "high")]
  y_vec <- all_data2[["Saver"]]@scale.data[gene, which(all_data2$level == "low")]
  
  if(diff(range(x_vec)) <= 1e-6 || diff(range(y_vec)) <= 1e-6 ) return(NA)
  
  test_res <- stats::wilcox.test(x = x_vec, y = y_vec)
  list(stat = mean(x_vec) - mean(y_vec), pvalue = test_res$p.value)
})
names(pvalue_list) <- all_data2[["Saver"]]@var.features
pvalue_list <- pvalue_list[sapply(pvalue_list, function(x){!all(is.na(x))})]
quantile(sapply(pvalue_list, function(x){x$stat}))
quantile(sapply(pvalue_list, function(x){-log10(x$pvalue)}))

pvalue_vec <- sapply(pvalue_list, function(x){x$pvalue})
names(pvalue_vec) <- names(pvalue_list)
pvalue_vec <- stats::p.adjust(pvalue_vec, method = "BH")
quantile(pvalue_vec)

# make a volcano plot now

source("../Writeup6b/gene_list.R")
keygene_vec <- sort(unique(c(keygenes$jackpot, keygenes$DABTRAM)))
keygene_vec <- sort(intersect(keygene_vec, names(pvalue_list)))

diff_vec <- sapply(pvalue_list, function(x){x$stat})
logpval_vec <- sapply(pvalue_list, function(x){-log10(x$pvalue)})
logpval_vec[is.infinite(logpval_vec)] <- max(logpval_vec[!is.infinite(logpval_vec)])
quantile(logpval_vec)

labeling_vec <- rep("1", length(diff_vec))
names(labeling_vec) <- names(diff_vec)
labeling_vec[names(logpval_vec)[which(logpval_vec >= 60)]] <- "3"
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
if(any(pvalue_vec <= 0.05)) {
  p1 <- p1 + ggplot2::geom_hline(yintercept=min(logpval_vec[which(pvalue_vec <= 0.05)]), linetype="dashed", 
                                 color = "red", linewidth=2)
}
p1 <- p1 + ggplot2::ggtitle(paste0("DABTRAM Day0 DE based on growth potential\n(Top 25% to Bottom 25%)")) +
  ggplot2::xlab("Mean difference in SAVER gene exp.: (Best-Worst)") + ggplot2::ylab("T test p-value (-Log10)")
p1 <- p1 + Seurat::NoLegend()

ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6k/Writeup6k_DABTRAM-day0_stepdown-LOOCV_DE_volcano.png"),
                p1, device = "png", width = 10, height = 10, units = "in")
