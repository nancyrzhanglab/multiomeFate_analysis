rm(list=ls())

reference <- read.csv("../../../../data/Sydney_stressors/Sydney_stressors_2021-09-24/starcode_barcodes/FeatureReference.csv")

treatment_files <- list.files(path = "../../../../data/Sydney_stressors/Sydney_stressors_2021-09-24/starcode_barcodes", 
                               pattern="output.*gDNA*", full.names = T)
treatment_list <- lapply(treatment_files, function(file){
  tab <- read.csv(file, sep = "\t", header = F)
  colnames(tab) <- c("sample", "reads")
  tab
})
names(treatment_list) <- c("Acid", "Cis", "CoCl2", "Dab", "Tram")

naive_files <- list.files(path = "../../../../data/Sydney_stressors/Sydney_stressors_2021-09-24/starcode_barcodes", 
                              pattern="*Naive*", full.names = T)
naive_list <- lapply(naive_files, function(file){
  tab <- read.csv(file, sep = "\t", header = F)
  colnames(tab) <- c("sample", "reads")
  tab
})
names(naive_list) <- c("Sample1", "Sample2", "Sample3")

startseq <- "GCTGTACAAGTAGGAT"

###################################

reference_original <- reference

# pair down reference to be only the samples that start with AT to make the loop run quicker
reference <- reference[(unlist(lapply(reference$sequence, function(x){ grepl('N',x, fixed = T)}))==0),]
reference <- reference[(unlist(lapply(reference$sequence, function(x){ grepl('AAAA',x, fixed = T)}))==0),]
reference <- reference[(unlist(lapply(reference$sequence, function(x){ grepl('TTTT',x, fixed = T)}))==0),]
reference <- reference[(unlist(lapply(reference$sequence, function(x){ grepl('CCCC',x, fixed = T)}))==0),]
reference <- reference[(unlist(lapply(reference$sequence, function(x){ grepl('GGGG',x, fixed = T)}))==0),] 

#############
treatment_list_org <- treatment_list
naive_list_org <- naive_list

# format treatments
for(i in 1:length(treatment_list)){
  start_vec <- sapply(treatment_list[[i]]$sample, function(x){
    substr(x, 1, 16)
  })
  treatment_list[[i]] <- treatment_list[[i]][which(start_vec == startseq),]
  
  sequence_vec <- sapply(treatment_list[[i]]$sample, function(x){
    substr(x, 17, 86)
  })
  treatment_list[[i]]$sample <- sequence_vec
  
  idx <- order(treatment_list[[i]]$sample)
  treatment_list[[i]] <- treatment_list[[i]][idx,]
}
sapply(treatment_list, nrow)
sapply(treatment_list, function(x){sum(x$reads)})

# format naives
for(i in 1:length(naive_list)){
  start_vec <- sapply(naive_list[[i]]$sample, function(x){
    substr(x, 1, 16)
  })
  naive_list[[i]] <- naive_list[[i]][which(start_vec == startseq),]
  
  sequence_vec <- sapply(naive_list[[i]]$sample, function(x){
    substr(x, 17, 86)
  })
  naive_list[[i]]$sample <- sequence_vec
  
  idx <- order(naive_list[[i]]$sample)
  naive_list[[i]] <- naive_list[[i]][idx,]
}
sapply(naive_list, nrow)
sapply(naive_list, function(x){sum(x$reads)})


#########

# check all the barcodes are the same
all(naive_list[[1]]$sample == naive_list[[2]]$sample)
all(naive_list[[1]]$sample == naive_list[[3]]$sample)
for(i in 1:length(treatment_list)){
  print(all(naive_list[[1]]$sample == treatment_list[[i]]$sample))
}

# keep only sequences in the reference
idx <- which(naive_list[[1]]$sample %in% reference$sequence)
for(i in 1:length(treatment_list)){
  treatment_list[[i]] <- treatment_list[[i]][idx,]
}
for(i in 1:length(naive_list)){
  naive_list[[i]] <- naive_list[[i]][idx,]
}

# remove the naive_mat rows with counts of 0, and subsequently all the barcodes in treatment_list
tmp <- sapply(naive_list, function(x){x$reads})
read_sum <- matrixStats::rowSums2(tmp)
idx <- which(read_sum > 0)
for(i in 1:length(treatment_list)){
  treatment_list[[i]] <- treatment_list[[i]][idx,]
}
for(i in 1:length(naive_list)){
  naive_list[[i]] <- naive_list[[i]][idx,]
}
sapply(treatment_list, function(x){sum(x$reads)})
sapply(naive_list, function(x){sum(x$reads)})

# normalize the counts
for(i in 1:length(treatment_list)){
  sum_val <- sum(treatment_list[[i]]$reads)
  treatment_list[[i]]$reads <- treatment_list[[i]]$reads/sum_val*1e6
}
for(i in 1:length(naive_list)){
  sum_val <- sum(naive_list[[i]]$reads)
  naive_list[[i]]$reads <- naive_list[[i]]$reads/sum_val*1e6
}
quantile(treatment_list[[1]]$reads)
naive_mat <- naive_list[[1]]
tmp <- sapply(naive_list, function(x){x$reads})
read_avg <- matrixStats::rowMeans2(tmp)
naive_mat$reads <- read_avg

############################

expansion_list <- lapply(1:length(treatment_list), function(i){
  vec <- log((treatment_list[[i]]$reads)/naive_mat$reads)
  vec <- vec[which(vec >= 0)]
  print(quantile(vec))
  EnvStats::epareto(vec)
})

############################

quantile(naive_mat$reads)
quantile(log(naive_mat$reads))
for(i in 1:length(treatment_list)){
  print(quantile(log((treatment_list[[i]]$reads + .1)/naive_mat$reads)))
}

#######################

png("../../../../out/figures/Writeup3e/10272021_sydney_barcode_naive.png",
    height = 1200, width = 1600, res = 300, units = "px")
hist(log(naive_mat$reads), col = "gray", xlab = "Log normalized reads", main = "Naive barcode reads\n(Averaged over 3 samples, and depth-normalized)")
graphics.off()

for(i in 1:length(treatment_list)){
  png(paste0("../../../../out/figures/Writeup3e/10272021_sydney_barcode_hist_naive-to-", names(treatment_list)[i], ".png"),
      height = 1200, width = 1600, res = 300, units = "px")
  len <- length(which(treatment_list[[i]]$reads > 0))
  hist(log((treatment_list[[i]]$reads + .1)/naive_mat$reads), 
       col = "gray", xlab = "Log ratio of normalized reads", 
       main = paste0("Log(", names(treatment_list)[i], "/naive) reads\n(", round(len/nrow(treatment_list[[i]])*100, 2),"% of treatment reads > 0)"))
  graphics.off()
}

for(i in 1:length(treatment_list)){
  png(paste0("../../../../out/figures/Writeup3e/10272021_sydney_barcode_hist_naive-to-", names(treatment_list)[i], "_onlyexpand.png"),
      height = 1200, width = 1600, res = 300, units = "px")
  vec <- log((treatment_list[[i]]$reads)/naive_mat$reads)
  vec <- vec[which(vec >= 0)]
  len <- length(vec)
  hist(vec, 
       col = "gray", xlab = "Log ratio of normliazed reads", 
       main = paste0("Log(", names(treatment_list)[i], "/naive) reads\n(Only the ", len, " barcodes that expanded)"))
  graphics.off()
}

for(i in 1:length(treatment_list)){
  png(paste0("../../../../out/figures/Writeup3e/10272021_sydney_barcode_hist_naive-to-", names(treatment_list)[i], "_onlyexpand_withpareto.png"),
      height = 1200, width = 1600, res = 300, units = "px")
  vec <- log((treatment_list[[i]]$reads)/naive_mat$reads)
  vec <- vec[which(vec >= 0)]
  len <- length(vec)
  hist(vec, 
       col = "gray", xlab = "Log ratio of normliazed reads", 
       main = paste0("Log(", names(treatment_list)[i], "/naive) reads, ",
                     len, " barcodes expanded)\n(Pareto, Location: ", 
       round(expansion_list[[i]]$parameters["location"], 4), 
       ", Shape: ", round(expansion_list[[i]]$parameters["shape"], 4), ")"))
  x_vec <- seq(0, ceiling(max(vec)), length.out = 1000)[-1]
  y_vec <- EnvStats::dpareto(x_vec, 
                             location = expansion_list[[i]]$parameters["location"], 
                             shape = expansion_list[[i]]$parameters["shape"])
  y_binned <- sapply(1:ceiling(max(vec)), function(x){
    idx <- intersect(which(x_vec <= x), which(x_vec >= x-1))
    sum(y_vec[idx])
  })
  bin_vec <- sapply(1:ceiling(max(vec)), function(x){
    length(intersect(which(vec <= x), which(vec >= x-1)))
  })
  multiplier1 <- mean(bin_vec/y_binned)
  multiplier2 <- mean(sapply(1:ceiling(max(vec)), function(x){
    idx <- intersect(which(x_vec <= x), which(x_vec >= x-1))
    mean(sum(y_vec[idx])/y_vec[idx])
  }))
  y_vec <- multiplier1*multiplier2*y_vec
  lines(x_vec, y_vec, col = i, lwd = 2)
  graphics.off()
}

#####################################3

png(paste0("../../../../out/figures/Writeup3e/10272021_sydney_barcode_hist_naive-to-all_onlyexpand_withpareto.png"),
    height = 1200, width = 1600, res = 300, units = "px")
plot(NA, xlim = c(0, 12), ylim = c(0, 400), xlab = "Log ratio of normalized reads",
     ylab = "Frequency")
for(i in 1:length(treatment_list)){
  vec <- log((treatment_list[[i]]$reads)/naive_mat$reads)
  vec <- vec[which(vec >= 0)]
  len <- length(vec)
  
  x_vec <- seq(0, 12, length.out = 1000)[-1]
  y_vec <- EnvStats::dpareto(x_vec, 
                             location = expansion_list[[i]]$parameters["location"], 
                             shape = expansion_list[[i]]$parameters["shape"])
  y_binned <- sapply(1:ceiling(max(vec)), function(x){
    idx <- intersect(which(x_vec <= x), which(x_vec >= x-1))
    sum(y_vec[idx])
  })
  bin_vec <- sapply(1:ceiling(max(vec)), function(x){
    length(intersect(which(vec <= x), which(vec >= x-1)))
  })
  multiplier1 <- mean(bin_vec/y_binned)
  multiplier2 <- mean(sapply(1:ceiling(max(vec)), function(x){
    idx <- intersect(which(x_vec <= x), which(x_vec >= x-1))
    mean(sum(y_vec[idx])/y_vec[idx])
  }))
  y_vec <- multiplier1*multiplier2*y_vec
  lines(x_vec, y_vec, col = i, lwd = 2)
}

legend("topright", names(treatment_list), col=1:5, lwd=2, bty="n");
graphics.off()


#####################################3

for(i in 1:length(treatment_list)){
  png(paste0("../../../../out/figures/Writeup3e/10272021_sydney_barcode_scatter_naive-to-", names(treatment_list)[i], ".png"),
      height = 1200, width = 1200, res = 300, units = "px")
  idx <- which(treatment_list[[i]]$reads > 0)
  cor_val <- stats::cor( log(naive_mat$reads[idx]), log(treatment_list[[i]]$reads[idx]))
  plot(x = log(naive_mat$reads[idx]), y = log(treatment_list[[i]]$reads[idx]),
       pch = 16, asp = T, xlab = "Log naive reads", 
       ylab = paste0("Log ", names(treatment_list)[i], " reads"),
       col = rgb(0.5, 0.5, 0.5, 0.2),
       main = paste0("Log reads of positive ", names(treatment_list)[i], " vs. naive\nCorrelation: ",
                     round(cor_val,2)))
  lines(c(-1e5, 1e5), c(-1e5, 1e5), col = "red", lty = 2)
  graphics.off()
}

idx_list <- lapply(treatment_list, function(mat){
  which(mat$reads > 0)
})
idx_vec <- unique(unlist(idx_list))
tmp <- sapply(treatment_list, function(mat){
  mat$reads[idx_vec]
})
tmp <- cbind(tmp, naive_mat$reads[idx_vec])
colnames(tmp) <- c(names(treatment_list), "Naive")
tmp <- log(tmp + .1)
tmp <- as.data.frame(tmp)
plot1 <- GGally::ggpairs(tmp, 
                         lower = list(continuous = GGally::wrap("points", alpha = 0.1)))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup3e/10272021_sydney_barcode_scatter_pairs.png"),
                plot1, device = "png", width = 7, height = 7, units = "in")


png(paste0("../../../../out/figures/Writeup3e/10272021_sydney_barcode_hist_acid.png"),
    height = 1200, width = 1600, res = 300, units = "px")
idx <- which(treatment_list_org[[1]]$reads > 1)
hist(log(treatment_list_org[[1]]$reads[idx]), 
     col = "gray", xlab = "Log of reads", 
     main = paste0("Log acid reads (Not normalized,\nonly raw count larger than 1)"))
graphics.off()

#################

idx1 <- which(treatment_list[[1]]$reads <= 1e-6)
idx2 <- which(treatment_list[[2]]$reads <= 1e-6)
length(intersect(idx1, idx2))