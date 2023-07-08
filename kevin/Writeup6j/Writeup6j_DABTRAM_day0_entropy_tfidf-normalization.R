rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

treatment <- "DABTRAM"
load(paste0("../../../../out/kevin/Writeup6i/day0_cutmatrix_extract-", treatment, "_tmp.RData"))

source("entropy_normalization.R")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

########

cutmat_list <- cutmat_list[which(sapply(cutmat_list, length) > 0)]

# idx <- which(names(cutmat_list) == "CD44")

len <- length(cutmat_list)
preprocessed_list <- vector("list", len)
names(preprocessed_list) <- names(cutmat_list)
for(idx in 1:len){
  print(idx)
  
  if(idx %% floor(len/10) == 0) {
    save(date_of_run, session_info,
         preprocessed_list,
         file = "../../../../out/kevin/Writeup6j/Writeup6j_DABTRAM_day0_entropy_tfidf-normalization.RData")
  }
  
  cutmat <- rbind(cutmat_list[[idx]]$cutmat_dying,
                  cutmat_list[[idx]]$cutmat_winning)
  peak_mat <- cutmat_list[[idx]]$peak_mat
  stopifnot(is.matrix(peak_mat))
  
  peak_prior <- cutmat_list[[idx]]$peak_prior
  thres1 <- 0.1; thres2 <- 0.05
  if(any(peak_prior >= thres1)){
    peak_mat <- peak_mat[which(peak_prior >= thres1),,drop = F]
  } else if(any(peak_prior >= thres2)){
    peak_mat <- peak_mat[which(peak_prior >= thres2),,drop = F]
  }
  
  peak_width <- stats::median(apply(peak_mat, 1, diff))
  data_mat <- .normalize_fragments(cutmat = cutmat,
                                   peak_mat = peak_mat,
                                   normalize_by_unique_cells = T)
  
  preprocessed_list[[idx]] <- list(
    idx_dying = which(data_mat$cell_name %in% rownames(cutmat_list[[idx]]$cutmat_dying)),
    idx_winning = which(data_mat$cell_name %in% rownames(cutmat_list[[idx]]$cutmat_winning)),
    data_mat = data_mat,
    peak_width = peak_width
  )
}

save(date_of_run, session_info,
     preprocessed_list,
     file = "../../../../out/kevin/Writeup6j/Writeup6j_DABTRAM_day0_entropy_tfidf-normalization.RData")

# break_vec <- sort(c(0, peak_width, 2*peak_width, 5*peak_width, max(data_mat$dist_to_peak)))
# table(cut(data_mat$dist_to_peak, breaks = break_vec),
#       cut(data_mat$tfidf, breaks = quantile(data_mat$tfidf)))
# 
# table(cut(data_mat$nearby_frags, breaks = quantile(data_mat$nearby_frags)),
#       cut(data_mat$tfidf, breaks = quantile(data_mat$tfidf)))


lrt_vec <- sapply(1:len, function(idx){
  print(idx)
  
  data_mat <- preprocessed_list[[idx]]$data_mat
  idx_dying <- preprocessed_list[[idx]]$idx_dying
  idx_winning <- preprocessed_list[[idx]]$idx_winning
  peak_width <- preprocessed_list[[idx]]$peak_width
  
  if(length(idx_dying) == 0 | length(idx_winning) == 0) return(NA)
  
  grenander_win <- multiomeFate:::estimate_grenander(
    values = data_mat$dist_to_peak[idx_winning],
    weights = data_mat$tfidf[idx_winning],
    # weights = rep(1, length(idx_winning)),
    scaling_factor = peak_width
  )
  loglik_win <- sum(sapply(data_mat$dist_to_peak[idx_winning], function(x){
    multiomeFate:::evaluate_grenander(
      obj = grenander_win,
      x = x,
      bool_log = T
    )
  }))
  
  grenander_die <- multiomeFate:::estimate_grenander(
    values = data_mat$dist_to_peak[idx_dying],
    weights = data_mat$tfidf[idx_dying],
    # weights = rep(1, length(idx_dying)),
    scaling_factor = peak_width
  )
  loglik_die <- sum(sapply(data_mat$dist_to_peak[idx_dying], function(x){
    multiomeFate:::evaluate_grenander(
      obj = grenander_die,
      x = x,
      bool_log = T
    )
  }))
  
  grenander_all <- multiomeFate:::estimate_grenander(
    values = data_mat$dist_to_peak,
    weights = data_mat$tfidf,
    # weights = rep(1, length(data_mat$dist_to_peak)),
    scaling_factor = peak_width
  )
  loglik_all <- sum(sapply(data_mat$dist_to_peak, function(x){
    multiomeFate:::evaluate_grenander(
      obj = grenander_all,
      x = x,
      bool_log = T
    )
  }))
  
  -2*(loglik_all - (loglik_win+loglik_die))
})
names(lrt_vec) <- names(cutmat_list) 

save(date_of_run, session_info,
     preprocessed_list, lrt_vec,
     file = "../../../../out/kevin/Writeup6j/Writeup6j_DABTRAM_day0_entropy_tfidf-normalization.RData")

#####################

length(which(lrt_vec <= 0))
quantile(lrt_vec, na.rm = T)
mean(lrt_vec, na.rm = T)

lrt_vec2 <- lrt_vec[intersect(which(lrt_vec >= 0), which(!is.na(lrt_vec)))]

df <- mean(lrt_vec2[lrt_vec2 <= stats::quantile(lrt_vec2, probs = 0.95)])

pval_vec <- stats::pchisq(lrt_vec2, df = df, lower.tail = F, log.p = F)
names(pval_vec) <- names(lrt_vec2)
quantile(pval_vec)
quantile(-log10(pval_vec), probs = seq(0.9, 1, length.out = 11))
pval_adj_vec <- stats::p.adjust(pval_vec, method = "BH")
quantile(pval_adj_vec)
length(which(pval_adj_vec <= 0.05))
pval_adj_vec[which(pval_adj_vec <= 0.05)]

x_vec <- seq(0, max(lrt_vec2), length.out = 1000)
y_vec <- stats::dchisq(x_vec, df = df)
png("../../../../out/figures/Writeup6j/Writeup6j_DABTRAM-day0_entropy_tfidf-normalization_lrt-teststat.png", 
    width = 3000, height = 1500, units = "px", res = 300)
par(mfrow = c(1,2))
hist(lrt_vec2, breaks = 50, col = "gray", 
     xlab = "LRT test stat", 
     main = "DABTRAM Day0 Grenander LRT")
y_vec <- y_vec*(max(res$counts)/max(y_vec))
lines(x = x_vec, y = y_vec, lwd = 2, lty = 2, col = "red")

hist(pval_vec, breaks = 50, col = "gray",
     main = "P-values based on empirical null",
     xlab = "P-value")
graphics.off()

#############################

preprocessed_list2 <- preprocessed_list[names(lrt_vec2)]
len <- length(preprocessed_list2)

.compute_halfcdf <- function(obj){
  idx <- max(which(obj$x <= 1))
  if(idx == length(obj$x)) return(1)
  
  x_new <- obj$x
  pdf_new <- obj$pdf
  
  x_new <- c(x_new[1:idx], 1, x_new[(idx+1):length(x_new)])
  pdf_new <- c(pdf_new[1:idx], pdf_new[idx:length(pdf_new)])
  
  cdf_vec <- cumsum(diff(x_new)*pdf_new[-length(pdf_new)])
  
  tab <- cbind(x_new[-1], cdf_vec)
  tab[which.min(abs(tab[,1] - 1)), 2]
}

entropy_vec <- sapply(1:len, function(idx){
  if(idx %% floor(len/10) == 0) cat('*')
  
  data_mat <- preprocessed_list2[[idx]]$data_mat
  idx_dying <- preprocessed_list2[[idx]]$idx_dying
  idx_winning <- preprocessed_list2[[idx]]$idx_winning
  peak_width <- preprocessed_list2[[idx]]$peak_width
  
  if(length(idx_dying) == 0 | length(idx_winning) == 0) return(NA)
  
  grenander_win <- multiomeFate:::estimate_grenander(
    values = data_mat$dist_to_peak[idx_winning],
    weights = data_mat$tfidf[idx_winning],
    scaling_factor = peak_width
  )
  
  grenander_die <- multiomeFate:::estimate_grenander(
    values = data_mat$dist_to_peak[idx_dying],
    weights = data_mat$tfidf[idx_dying],
    scaling_factor = peak_width
  )
  
  entropy_win <- 1-.compute_halfcdf(grenander_win)
  entropy_die <- 1-.compute_halfcdf(grenander_die)
  
  (entropy_win - entropy_die)/max(entropy_die, 1e-3)
})
names(entropy_vec) <- names(preprocessed_list2)

entropy_vec[names(pval_adj_vec)[which(pval_adj_vec <= 0.05)]]
quantile(entropy_vec)
quantile(entropy_vec, probs = seq(0.9, 1, length.out = 5))


###############

source("../Writeup6b/gene_list.R")
keygene_vec <- sort(unique(c(keygenes$jackpot, keygenes$DABTRAM)))

labeling_vec <- rep("1", length(entropy_vec))
names(labeling_vec) <- names(entropy_vec)
labeling_vec[names(pval_adj_vec)[which(pval_adj_vec <= .05)]] <- "3"
labeling_vec[intersect(keygene_vec, names(labeling_vec))] <- "2"

df <- data.frame(diff_entropy = entropy_vec,
                 log10pval = -log10(pval_vec),
                 name = names(entropy_vec),
                 labeling = labeling_vec)
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
p1 <- p1 + ggplot2::ggtitle(paste0("DABTRAM Day0 (TFIDF-normalization)")) +
  ggplot2::xlab("Entropy difference: (Win-Lose)") + ggplot2::ylab("LRT empirical null p-value (Log10)")
p1 <- p1 + Seurat::NoLegend()

ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6j/Writeup6j_DABTRAM-day0_entropy_tfidf-normalization_volcano.png"),
                p1, device = "png", width = 10, height = 10, units = "in")

###################

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

gene_vec <- names(pval_adj_vec)[intersect(
  which(pval_adj_vec <= .05),
  which(entropy_vec > 0)
)]

pdf(paste0("../../../../out/figures/Writeup6j/Writeup6j_DABTRAM-day0_entropy_tfidf-normalization_coverage.pdf"),
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

