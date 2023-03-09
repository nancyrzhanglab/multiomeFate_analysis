rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(hexbin)
library(RColorBrewer)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")
load("../../../../out/kevin/Writeup6c/Writeup6c_peak_counter_5000_onpeak.RData")
countmat_peak <- countmat_nopeak
load("../../../../out/kevin/Writeup6c/Writeup6c_peak_counter_5000.RData")

countmat_peak <- Matrix::Matrix(countmat_peak, sparse = T)
countmat_nopeak <- Matrix::Matrix(countmat_nopeak, sparse = T)

dim(countmat_peak)
dim(countmat_nopeak)
length(countmat_peak@x)/prod(dim(countmat_peak))
length(countmat_nopeak@x)/prod(dim(countmat_nopeak))

quantile(countmat_peak@x, probs = seq(0,1,length.out=11))
quantile(countmat_nopeak@x, probs = seq(0,1,length.out=11))

####################################

.nonzero_col <- function(mat, col_idx, bool_value){
  stopifnot(inherits(mat, "dgCMatrix"), col_idx %% 1 == 0,
            col_idx > 0, col_idx <= ncol(mat))
  
  val1 <- mat@p[col_idx]
  val2 <- mat@p[col_idx+1]
  
  if(val1 == val2) return(numeric(0))
  if(bool_value){
    # return the value
    mat@x[(val1+1):val2]
  } else {
    # return the row index
    mat@i[(val1+1):val2]+1
  }
}

.mult_mat_vec <- function(mat, vec){
  stopifnot(inherits(mat, c("matrix", "dgCMatrix")), 
            !is.matrix(vec), length(vec) == ncol(mat))
  
  if(inherits(mat, "dgCMatrix")) {
    mat %*% Matrix::Diagonal(x = vec)
  } else {
    mat * rep(vec, rep(nrow(mat), length(vec)))
  }
}

####################################

p <- ncol(countmat_peak); n <- nrow(countmat_peak)
round(quantile(sapply(1:p, function(j){length(.nonzero_col(countmat_peak, j, F))}), 
               probs = seq(0,1,length.out=11)))
round(quantile(sapply(1:p, function(j){length(.nonzero_col(countmat_nopeak, j, F))}), 
               probs = seq(0,1,length.out=11)))

ratio_mat <- (countmat_nopeak)/(countmat_peak+1)
value_vec <- sapply(1:p, function(j){
  num_count <- length(.nonzero_col(countmat_nopeak, j, F))
  log(n/num_count)
})

norm_mat <- .mult_mat_vec(ratio_mat, value_vec)
length(norm_mat@x)/prod(dim(norm_mat))
norm_mat@x <- log1p(norm_mat@x)

####################################################

all_data[["offpeak"]] <- Seurat::CreateAssayObject(counts = Matrix::t(countmat_nopeak))

lin_mat <- table(all_data$assigned_lineage, all_data$dataset)
DABTRAM_lineage <- rownames(lin_mat)[which(lin_mat[,"day10_DABTRAM"] >= 20)]
day0_survive_idx <- intersect(which(all_data$assigned_lineage %in% DABTRAM_lineage),
                              which(all_data$dataset == "day0"))
day0_die_idx <- setdiff(intersect(which(!is.na(all_data$assigned_lineage)),
                                  which(all_data$dataset == "day0")),
                        day0_survive_idx)
day0_unknown_idx <- setdiff(which(all_data$dataset == "day0"),
                            c(day0_survive_idx, day0_die_idx))

tmp_vec <- rep(NA, ncol(all_data))
tmp_vec[day0_survive_idx] <- "day0_survive"
tmp_vec[day0_die_idx] <- "day0_die"
tmp_vec[day0_unknown_idx] <- "day0_unknown"
all_data$day0_DABTRAM_status <- tmp_vec

cell_idx <- sort(unique(c(day0_survive_idx, day0_die_idx, day0_unknown_idx)))
cell_idx2 <- sort(unique(c(day0_survive_idx, day0_die_idx)))
lineage_vec <- all_data$assigned_lineage[cell_idx2]
fitness_vec_onlyday0 <- sapply(1:length(cell_idx2), function(i){
  log1p(lin_mat[lineage_vec[i],"day10_DABTRAM"])
})
fitness_vec <- rep(NA, ncol(all_data))
fitness_vec[cell_idx2] <- fitness_vec_onlyday0
all_data$day0_DABTRAM_fitness <- fitness_vec

n <- ncol(all_data)

cell_idx <- sort(unique(c(day0_survive_idx, day0_die_idx, day0_unknown_idx)))
survive_idx <- which(cell_idx %in% day0_survive_idx)
die_idx <- which(cell_idx %in% day0_die_idx)
unknown_idx <- which(cell_idx %in% day0_unknown_idx)

fitness_vec2 <- all_data$day0_DABTRAM_fitness[cell_idx2]
survive_idx2 <- which(cell_idx2 %in% day0_survive_idx)
die_idx2 <- which(cell_idx2 %in% day0_die_idx)

length(survive_idx) == length(survive_idx2)
length(die_idx) == length(die_idx2)

####################################################

norm_mat_t <- Matrix::t(norm_mat)
active_gene_vec <- sapply(cell_idx, function(i){
  sum(
    .nonzero_col(norm_mat_t, col_idx = i, bool_value = T)
  )
})
peak_plotting_func(active_gene_vec = active_gene_vec, 
                   fitness_vec2 = fitness_vec2, 
                   survive_idx = survive_idx, die_idx = die_idx, unknown_idx = unknown_idx,
                   survive_idx2 = survive_idx2, die_idx2 = die_idx2,
                   filename = "../../../../out/figures/Writeup6c/Writeup6c_DABTRAM_peakcounter_plotter_custnorm1-offpeak.png",
                   condition_name = "DABTRAM",
                   main4 = "cust-norm of offpeaks", 
                   xlab = "Custom norm.1 of offpeaks", 
                   ylab4 = "Custom norm.1 of offpeaks")

#######################################################

countmat_peak_dense <- as.matrix(countmat_peak)
countmat_nopeak_dense <- as.matrix(countmat_nopeak)

# Make the plot
my_colors <- grDevices::colorRampPalette(rev(brewer.pal(11,'Spectral')))
range_vec <- c(-5, quantile(c(countmat_peak@x, countmat_nopeak@x), probs = 0.95)+1)
range_vec[2] <- max(range_vec[2], 50)

png("../../../../out/figures/Writeup6c/Writeup6c_DABTRAM_peakcounter_hexbin_allgenes_surviving.png",
    height = 1500, width = 1500, units = "px", res = 300)
x_vec1 <- as.numeric(countmat_peak_dense[day0_survive_idx,])
y_vec1 <- as.numeric(countmat_nopeak_dense[day0_survive_idx,])
idx <- unique(c(which(x_vec1 >= range_vec[2]), which(y_vec1 >= range_vec[2])))
x_vec1 <- x_vec1[-idx]; y_vec1 <- y_vec1[-idx]
bin <- hexbin(x_vec1, y_vec1, 
              xbins = 40,
              xbnds = range_vec, 
              ybnds = range_vec)
plot(bin, main = "Surviving cells w/ lineage", colramp = my_colors, 
     legend = F,
     xlab = "On-peak count",
     ylab = "Off-peak count",
     trans = function(x){log(x+1)}, inv = function(x){exp(x)-1}) 
graphics.off()

png("../../../../out/figures/Writeup6c/Writeup6c_DABTRAM_peakcounter_hexbin_allgenes_dying.png",
    height = 1500, width = 1500, units = "px", res = 300)
x_vec1 <- as.numeric(countmat_peak_dense[day0_die_idx,])
y_vec1 <- as.numeric(countmat_nopeak_dense[day0_die_idx,])
idx <- unique(c(which(x_vec1 >= range_vec[2]), which(y_vec1 >= range_vec[2])))
x_vec1 <- x_vec1[-idx]; y_vec1 <- y_vec1[-idx]
bin <- hexbin(x_vec1, y_vec1, 
              xbins = 40,
              xbnds = range_vec, 
              ybnds = range_vec)
plot(bin, main = "Dying cells w/ lineage", colramp = my_colors, 
     legend = F,
     xlab = "On-peak count",
     ylab = "Off-peak count",
     trans = function(x){log(x+1)}, inv = function(x){exp(x)-1}) 
graphics.off()

png("../../../../out/figures/Writeup6c/Writeup6c_DABTRAM_peakcounter_hexbin_allgenes_unknown.png",
    height = 1500, width = 1500, units = "px", res = 300)
x_vec1 <- as.numeric(countmat_peak_dense[c(day0_die_idx, day0_unknown_idx),])
y_vec1 <- as.numeric(countmat_nopeak_dense[c(day0_die_idx, day0_unknown_idx),])
idx <- unique(c(which(x_vec1 >= range_vec[2]), which(y_vec1 >= range_vec[2])))
x_vec1 <- x_vec1[-idx]; y_vec1 <- y_vec1[-idx]
bin <- hexbin(x_vec1, y_vec1, 
              xbins = 40,
              xbnds = range_vec, 
              ybnds = range_vec)
plot(bin, main = "Dying cells w/+w/o lineage", colramp = my_colors, 
     legend = F,
     xlab = "On-peak count",
     ylab = "Off-peak count",
     trans = function(x){log(x+1)}, inv = function(x){exp(x)-1}) 
graphics.off()


