rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(IRanges)
library(GGally)
library(ggplot2)

load("../../../../out/kevin/Writeup6l/Writeup6l_chromVar_rna-chromvar_lightweight.RData")
tab_mat_full <- table(all_data$assigned_lineage, all_data$dataset)

load("../../../../out/kevin/Writeup6l/Writeup6l_day0_chromvar_lightweight_noATAC.RData")

tab_vec <- table(all_data$assigned_lineage)
tab_vec <- tab_vec[tab_vec > 0]
lineage_names <- names(tab_vec)
tab_mat_full <- tab_mat_full[lineage_names,]

# compute which cells belong to which lineage
lineage_idx_list <- lapply(lineage_names, function(lineage){
  intersect(which(all_data$assigned_lineage == lineage),
            which(all_data$dataset == "day0"))
})
names(lineage_idx_list) <- lineage_names

# for each lineage, compute the mean scores
name_vec <- names(data.use) 
mean_mat <- sapply(1:length(name_vec), function(i){
  if(i %% floor(length(name_vec)/10) == 0) cat('*')
  sapply(lineage_idx_list, function(idx_vec){
    mean(all_data[["chromvar"]]@data[i,idx_vec], na.rm = T)
  })
})
rownames(mean_mat) <- lineage_names
colnames(mean_mat) <- name_vec

all(rownames(mean_mat) == rownames(tab_mat_full))
day0 <- log10(tab_mat_full[,"day0"]+1)
day10_cocl2 <- log10(tab_mat_full[,"day10_COCL2"]+1)

cor_day0 <- sapply(1:ncol(mean_mat), function(j){
  x <- day0
  y <- mean_mat[,j]
  rm_idx <- sort(unique(c(which(is.na(x)), which(is.na(y)))))
  if(length(rm_idx) > 0){x <- x[-rm_idx]; y <- y[-rm_idx]}
  stats::cor(x, y)
})
names(cor_day0) <- name_vec

cor_day10 <- sapply(1:ncol(mean_mat), function(j){
  x <- day10_cocl2
  y <- mean_mat[,j]
  rm_idx <- sort(unique(c(which(is.na(x)), which(is.na(y)))))
  if(length(rm_idx) > 0){x <- x[-rm_idx]; y <- y[-rm_idx]}
  stats::cor(x, y)
})
names(cor_day10) <- name_vec

quantile(cor_day0); quantile(cor_day10)

motif_idx <- sort(unique(c(grep("JUN", name_vec),
                           grep("FOS", name_vec),
                           grep("NFE2", name_vec),
                           grep("TEAD", name_vec))))
labeling_vec <- rep(0, length(cor_day0))
labeling_vec[motif_idx] <- 1
labeling_vec[order(cor_day0, decreasing = T)[1:20]] <- 2
labeling_vec[order(cor_day10, decreasing = T)[1:20]] <- 2

df <- data.frame(cor_day0 = cor_day0,
                 cor_day10 = cor_day10,
                 name = name_vec,
                 labeling = labeling_vec)
df[,"labeling"] <- as.factor(df[,"labeling"])
# put all the labeling == TRUE on bottom
df <- df[c(which(df[,"labeling"] == 0), which(df[,"labeling"] == 1), which(df[,"labeling"] == 2)),]

p1 <- ggplot2::ggplot(df, ggplot2::aes(x = cor_day0, 
                                       y = cor_day10))
p1 <- p1 + ggplot2::geom_point(ggplot2::aes(color = labeling))
p1 <- p1 + ggplot2::scale_colour_manual(values=c("black", "green",  "red"))
p1 <- p1 + ggrepel::geom_text_repel(data = subset(df, labeling == "2"),
                                    ggplot2::aes(label = name, color = labeling),
                                    box.padding = ggplot2::unit(0.5, 'lines'),
                                    point.padding = ggplot2::unit(1.6, 'lines'),
                                    max.overlaps = 50)
p1 <- p1 + ggplot2::ggtitle(paste0("COCL2")) +
  ggplot2::xlab("Correlation b/w chromVar and Day0 size") + 
  ggplot2::ylab("Correlation b/w chromVar and Day10 size")
p1 <- p1 + Seurat::NoLegend()

ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6l/Writeup6l_COCL2_chromvar_day0-size_day10-size_all-motifs.png"),
                p1, device = "png", width = 10, height = 10, units = "in")

################################################

# for each lineage, compute the max scores
name_vec <- names(data.use) 
max_mat <- sapply(1:length(name_vec), function(i){
  if(i %% floor(length(name_vec)/10) == 0) cat('*')
  sapply(lineage_idx_list, function(idx_vec){
    max(all_data[["chromvar"]]@data[i,idx_vec], na.rm = T)
  })
})
rownames(max_mat) <- lineage_names
colnames(max_mat) <- name_vec

all(rownames(max_mat) == rownames(tab_mat_full))
day0 <- log10(tab_mat_full[,"day0"]+1)
day10_cocl2 <- log10(tab_mat_full[,"day10_COCL2"]+1)

cor_day0 <- sapply(1:ncol(max_mat), function(j){
  x <- day0
  y <- max_mat[,j]
  rm_idx <- sort(unique(c(which(is.na(x)), which(is.na(y)))))
  if(length(rm_idx) > 0){x <- x[-rm_idx]; y <- y[-rm_idx]}
  stats::cor(x, y)
})
names(cor_day0) <- name_vec

cor_day10 <- sapply(1:ncol(max_mat), function(j){
  x <- day10_cocl2
  y <- max_mat[,j]
  rm_idx <- sort(unique(c(which(is.na(x)), which(is.na(y)))))
  if(length(rm_idx) > 0){x <- x[-rm_idx]; y <- y[-rm_idx]}
  stats::cor(x, y)
})
names(cor_day10) <- name_vec

quantile(cor_day0); quantile(cor_day10)

motif_idx <- sort(unique(c(grep("JUN", name_vec),
                           grep("FOS", name_vec),
                           grep("NFE2", name_vec),
                           grep("TEAD", name_vec))))
labeling_vec <- rep(0, length(cor_day0))
labeling_vec[order(cor_day0, decreasing = T)[1:20]] <- 1
labeling_vec[order(cor_day10, decreasing = T)[1:20]] <- 1
labeling_vec[motif_idx] <- 2

df <- data.frame(cor_day0 = cor_day0,
                 cor_day10 = cor_day10,
                 name = name_vec,
                 labeling = labeling_vec)
df[,"labeling"] <- as.factor(df[,"labeling"])
# put all the labeling == TRUE on bottom
df <- df[c(which(df[,"labeling"] == 0), which(df[,"labeling"] == 1), which(df[,"labeling"] == 2)),]

p1 <- ggplot2::ggplot(df, ggplot2::aes(x = cor_day0, 
                                       y = cor_day10))
p1 <- p1 + ggplot2::geom_point(ggplot2::aes(color = labeling))
p1 <- p1 + ggplot2::scale_colour_manual(values=c("black", "red", "green"))
p1 <- p1 + ggrepel::geom_text_repel(data = subset(df, labeling == "1"),
                                    ggplot2::aes(label = name, color = labeling),
                                    box.padding = ggplot2::unit(0.5, 'lines'),
                                    point.padding = ggplot2::unit(1.6, 'lines'),
                                    max.overlaps = 50)
p1 <- p1 + ggplot2::ggtitle(paste0("COCL2")) +
  ggplot2::xlab("Correlation b/w chromVar (max) and Day0 size") + 
  ggplot2::ylab("Correlation b/w chromVar (max) and Day10 size")
p1 <- p1 + Seurat::NoLegend()

ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6l/Writeup6l_COCL2_chromvar-max_day0-size_day10-size_all-motifs.png"),
                p1, device = "png", width = 10, height = 10, units = "in")


