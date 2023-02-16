rm(list=ls())
load("../../dbox_MultiomeFate/ShafferLab_Data/2022_12_31_coverageplot_track.RData")

class(track_list)
names(track_list)
length(track_list)

class(track_list[["ANXA1"]])
names(track_list[["ANXA1"]])

colnames(track_list[["ANXA1"]][["coverage_sum"]])
dim(track_list[["ANXA1"]][["coverage_sum"]])
dim(track_list[["ANXA1"]][["coverage_count"]])
length(track_list[["ANXA1"]][["total_vec"]])

track_list[["ANXA1"]][["coverage_sum"]][1:5,1:5]
track_list[["ANXA1"]][["total_vec"]]

dim(track_list[["AXL"]][["coverage_sum"]])
dim(track_list[["AXL"]][["coverage_count"]])
length(track_list[["AXL"]][["total_vec"]])
track_list[["AXL"]][["total_vec"]]

k <- "EGR3"
zz <- track_list[[k]]$coverage_sum
plot(as.numeric(rownames(zz)), zz[,8], type = "l", 
     main = colnames(zz)[8],
     xlab = "Base pair", ylab = "Normalized signal")
lines(as.numeric(rownames(zz)), zz[,1], col = 2)

zz <- track_list[[k]]$coverage_count
par(mfrow = c(1,2))
plot(as.numeric(rownames(zz)), zz[,1], col = 2, type = "l", 
     main = colnames(zz)[1],
     xlab = "Base pair", ylab = "Raw count")
plot(as.numeric(rownames(zz)), zz[,8], type = "l", 
     main = colnames(zz)[8],
     xlab = "Base pair", ylab = "Raw count")

