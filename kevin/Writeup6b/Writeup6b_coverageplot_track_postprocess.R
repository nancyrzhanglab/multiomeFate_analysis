rm(list=ls())
load("../../out/Writeup6b/Writeup6b_coverageplot_track.RData")

zz <- track_list[[8]]$coverage_mean
plot(as.numeric(rownames(zz)), zz[,8], type = "l", main = colnames(zz)[8])
lines(as.numeric(rownames(zz)), 8*zz[,1], col = 2)

zz <- track_list[[8]]$coverage_mean
plot(as.numeric(rownames(zz)), zz[,8]/max(zz[,8]), type = "l", main = colnames(zz)[8])
lines(as.numeric(rownames(zz)), zz[,1]/max(zz[,1]), col = 2)

zz <- track_list[[8]]$coverage_count
plot(as.numeric(rownames(zz)), zz[,8]/track_list[[2]]$total_vec[8],
     type = "l", main = colnames(zz)[8])
lines(as.numeric(rownames(zz)), zz[,1]/track_list[[2]]$total_vec[1], col = 2)
