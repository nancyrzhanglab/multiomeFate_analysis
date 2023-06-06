rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6h/Writeup6h_DABTRAM-day10_pseudotime.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

cell_features <- cbind(1, scale(all_data$pseudotime[colnames(all_data2)]))
colnames(cell_features) <- c("Intercept", "Pseudotime")
cell_lineage <- all_data2$assigned_lineage
cell_lineage <- cell_lineage[!is.na(cell_lineage)]
cell_features <- cell_features[names(cell_lineage),]
uniq_lineage <- sort(unique(cell_lineage))
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
lineage_current_count <- tab_mat[uniq_lineage,"day10_DABTRAM"]
lineage_future_count <- tab_mat[uniq_lineage,"week5_DABTRAM"]

tmp <- quantile(abs(cell_features[,2]), probs = 0.95)
coef_val <- 2*log(max(lineage_future_count/lineage_current_count))/tmp
coefficient_initial <- c(0, coef_val)
names(coefficient_initial) <- colnames(cell_features)

set.seed(10)
res1 <- multiomeFate:::lineage_imputation(cell_features = cell_features,
                                          cell_lineage = cell_lineage,
                                          coefficient_initial = coefficient_initial,
                                          lineage_future_count = lineage_future_count,
                                          random_initializations = 10,
                                          verbose = 2)

#####

# now try using the somewhat correlated features of the topic model
topic_mat <- all_data[["fasttopic_DABTRAM"]]@cell.embeddings[names(cell_lineage),]
pseudotime_vec <- all_data$pseudotime[names(cell_lineage)]
cor_vec <- stats::cor(topic_mat, pseudotime_vec)
quantile(abs(cor_vec))

cell_features <- cbind(1, scale(topic_mat[,which(abs(cor_vec) >= 0.5)]), scale(pseudotime_vec))
p <- ncol(cell_features)
colnames(cell_features)[1] <- "Intercept"
colnames(cell_features)[p] <- "Pseudotime"

tmp <- quantile(abs(cell_features[,-1]), probs = 0.95)
coef_val <- 2*log(max(lineage_future_count/lineage_current_count))/((p-1)*tmp)
coefficient_initial <- c(0, rep(coef_val, p-1))
names(coefficient_initial) <- colnames(cell_features)

set.seed(10)
res2 <- multiomeFate:::lineage_imputation(cell_features = cell_features,
                                          cell_lineage = cell_lineage,
                                          coefficient_initial = coefficient_initial,
                                          lineage_future_count = lineage_future_count,
                                          random_initializations = 10,
                                          verbose = 2)

cell_imputed_count <- as.numeric(exp(cell_features %*% res2$fit$coefficient_vec))
names(cell_imputed_count) <- rownames(cell_features)
quantile(round(cell_imputed_count), probs = seq(0.9,1,length.out=11))
round(res2$fit$coefficient_vec, 3)

lineage_imputed_count <- sapply(uniq_lineage, function(lineage){
  sum(cell_imputed_count[which(cell_lineage == lineage)])
})
stats::cor(lineage_imputed_count, lineage_future_count)
stats::cor(log1p(lineage_imputed_count), log1p(lineage_future_count))

imputed_vec <- rep(NA, ncol(all_data))
names(imputed_vec) <- colnames(all_data)
imputed_vec[names(cell_imputed_count)] <- log(cell_imputed_count)
all_data$imputed_count <- imputed_vec

p1 <- Seurat::FeaturePlot(all_data, reduction = "umap", 
                          features = "imputed_count")
p1 <- p1 + ggplot2::ggtitle("DABTRAM Day10 imputed counts (Log-scale)")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6i/Writeup6i_DABTRAM-day10_imputed-lineage_umap.png"),
                p1, device = "png", width = 7, height = 5, units = "in")

png("../../../../out/figures/Writeup6i/Writeup6i_DABTRAM-day10_lineage-level_counts.png",
    height = 1500, width = 1500, res = 300, units = "px")
plot(log10(lineage_future_count+1), 
     log10(lineage_imputed_count+1),
     xlab = "Observed lineage count", ylab = "Predicted lineage count",
     main = paste0("DABTRAM Day10 imputed counts\n(Log-scale), Cor: ", 
                   round(stats::cor(log10(lineage_imputed_count+1), 
                                    log10(lineage_future_count+1)),2)),
     pch = 16, asp = T, col = rgb(0.5,0.5,0.5,0.5))
graphics.off()

png("../../../../out/figures/Writeup6i/Writeup6i_DABTRAM-day10_pseudotime-vs-imputed.png",
    height = 1500, width = 1500, res = 300, units = "px")
plot(cell_features[,"Pseudotime"], 
     log(cell_imputed_count),
     xlab = "RNA Pseudotime", ylab = "Log imputated count",
     main = paste0("DABTRAM Day10 imputed counts\n(Log-scale), Cor: ", 
                   round(stats::cor(cell_features[,"Pseudotime"], 
                                    log(cell_imputed_count)),2)),
     pch = 16, col = rgb(0.5,0.5,0.5,0.5))
graphics.off()

########

all(names(lineage_imputed_count) == names(lineage_future_count))

lineage_imputed_count2 <- log10(lineage_imputed_count+1)
lineage_future_count2 <- log10(lineage_future_count+1)

labeling_vec <- rep(FALSE, length(lineage_imputed_count2))
labeling_vec[which(lineage_imputed_count2 >= 1.5)] <- TRUE
labeling_vec[which(lineage_future_count2 >= 1.5)] <- TRUE
table(labeling_vec)

df <- data.frame(lineage_imputed_count = lineage_imputed_count2,
                 lineage_future_count = lineage_future_count2,
                 name = names(lineage_imputed_count2),
                 labeling = labeling_vec)
# put all the labeling == TRUE on bottom
df <- df[c(which(!df[,"labeling"]), which(df[,"labeling"])),]

p1 <- ggplot2::ggplot(df, ggplot2::aes(x = lineage_future_count, y = lineage_imputed_count))
p1 <- p1 + ggplot2::geom_point(ggplot2::aes(color = labeling))
p1 <- p1 + ggplot2::scale_colour_manual(values=c("black", "red"))
p1 <- p1 + ggrepel::geom_text_repel(data = subset(df, labeling == TRUE),
                                    ggplot2::aes(label = name, color = labeling),
                                    box.padding = ggplot2::unit(0.5, 'lines'),
                                    point.padding = ggplot2::unit(1.6, 'lines'),
                                    max.overlaps = 50)
p1 <- p1 + ggplot2::ggtitle("DABTRAM Day10 imputed counts (Log-scale)") +
  ggplot2::xlab("Observed lineage count (Log10)") + ggplot2::ylab("Predicted lineage count (Log10)")
p1 <- p1 + Seurat::NoLegend()

ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6i/Writeup6i_DABTRAM-day10_lineage-level_counts_labeled.png"),
                p1, device = "png", width = 10, height = 10, units = "in")

###################################

# https://github.com/linnykos/permanent_notes/blob/master/negative-binomial/overdispersion.pdf
var_vec <- (lineage_imputed_count-lineage_future_count)^2
overdisperion_stat <- sum(var_vec/lineage_imputed_count)
df <- (nrow(cell_features) - ncol(cell_features))
overdisperion_stat <- overdisperion_stat/df

overdisperion_stat

################################

uniq_lineages <- sort(unique(names(lineage_future_count)))
cell_lineage_idx_list <- lapply(uniq_lineages, function(lineage){
  which(cell_lineage == lineage)
})
names(cell_lineage_idx_list) <- uniq_lineages

# plot lines and see what the objective value is like
coefficient_vec <- res2$fit$coefficient_vec

perturb_list <- vector("list", length = ncol(cell_features))
adj_vec <- seq(-2,2,length.out=1000)

for(j in 1:ncol(cell_features)){
  print(j)
  obj_vec <- sapply(adj_vec, function(adj){
    tmp <- coefficient_vec
    tmp[j] <- tmp[j]+adj
    
    multiomeFate:::.lineage_objective(cell_features = cell_features,
                                      cell_lineage = cell_lineage,
                                      cell_lineage_idx_list = cell_lineage_idx_list,
                                      coefficient_vec = tmp,
                                      lineage_future_count = lineage_future_count)
  })
  
  perturb_list[[j]] <- list(obj_vec = obj_vec,
                            x_vec = coefficient_vec[j]+adj_vec)
}

sapply(perturb_list, function(lis){
  quantile(lis$obj_vec)
})

pdf("../../../../out/figures/Writeup6i/Writeup6i_DABTRAM-day10_imputation_objective.pdf", 
    onefile = T, height = 5, width = 10)
for(j in 1:ncol(cell_features)){
  ylim <- stats::quantile(perturb_list[[j]]$obj_vec, probs = c(0, 0.75))
  len <- length(perturb_list[[j]]$x_vec)
  idx <- round(len/2)
  plot(x = perturb_list[[j]]$x_vec,
       y = perturb_list[[j]]$obj_vec,
       ylim = ylim,
       pch = 16, cex = 0.5,
       xlab = colnames(cell_features)[j],
       ylab = "Objective function")
  
  points(x = perturb_list[[j]]$x_vec[idx],
         y = perturb_list[[j]]$obj_vec[idx],
         pch = 16, cex = 2, col = "red")
}
graphics.off()