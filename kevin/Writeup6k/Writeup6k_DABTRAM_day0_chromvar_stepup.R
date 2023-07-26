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

all_data_original <- all_data
load("../../../../out/kevin/Writeup6b/Writeup6b_chromVar.RData")

# do the chromvar analysis 
tab_mat <- table(all_data_original$assigned_lineage, all_data_original$dataset)

#################

cell_winning_idx <- intersect(which(all_data_original$imputed_count > quantile(all_data_original$imputed_count, probs = 0.75, na.rm = T)),
                              which(all_data_original$dataset == "day0"))
cell_winning_names <- colnames(all_data_original)[cell_winning_idx]
cell_losing_idx <- intersect(which(all_data_original$imputed_count < quantile(all_data_original$imputed_count, probs = 0.25, na.rm = T)),
                             which(all_data_original$dataset == "day0"))
cell_losing_names <- colnames(all_data_original)[cell_losing_idx]
length(cell_winning_names); length(cell_losing_names)
keep_vec <- rep(NA, ncol(all_data))
names(keep_vec) <- colnames(all_data)
keep_vec[cell_winning_names] <- "winning_day0"
keep_vec[cell_losing_names] <- "losing_day0"

all_data$ident <- keep_vec
Seurat::Idents(all_data) <- "ident"

Seurat::DefaultAssay(all_data) <- "chromvar"
set.seed(10)
de_res <- Seurat::FindMarkers(
  object = all_data,
  ident.1 = "winning_day0",
  ident.2 = "losing_day0",
  only.pos = FALSE,
  min.pct = 0,
  logfc.threshold = 0,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  verbose = F
)

motifs <- rownames(de_res)
data.use <- Signac::GetMotifData(object = all_data,
                                 assay = "ATAC",
                                 slot = "pwm")
data.use <- data.use[motifs]
names(data.use) <- Signac::GetMotifData(
  object = all_data,
  assay = "ATAC",
  slot = "motif.names"
)[motifs]
tmp <- de_res; rownames(tmp) <- names(data.use); tmp[1:50,]

names(data.use) <- sapply(1:length(data.use), function(i){
  paste0(names(data.use)[i], ", ", rownames(de_res)[i], 
         "\n-Log10AdjPval=", round(-log10(de_res[i,"p_val_adj"]), 2),
         ", value=", round(de_res[i,"avg_diff"], 2))
})
data.use <- data.use[1:28]

plot1 <- ggseqlogo::ggseqlogo(data = data.use, ncol = 4)
plot1 <- plot1 + ggplot2::theme_bw()
width <- 15; height <- 10

ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6k/Writeup6k_DABTRAM_day0_chromvar_stepup.png"),
                plot1, device = "png", width = width, height = height, units = "in")
