rm(list=ls())
load("../../out/Writeup6b/Writeup6b_coverageplot_track.RData")

k <- 12
zz <- track_list[[k]]$coverage_sum
plot(as.numeric(rownames(zz)), zz[,8], type = "l", main = colnames(zz)[8])
lines(as.numeric(rownames(zz)), zz[,1], col = 2)

zz <- track_list[[k]]$coverage_count
par(mfrow = c(1,2))
plot(as.numeric(rownames(zz)), zz[,1], col = 2, type = "l", main = colnames(zz)[8])
plot(as.numeric(rownames(zz)), zz[,8], type = "l", main = colnames(zz)[8])

