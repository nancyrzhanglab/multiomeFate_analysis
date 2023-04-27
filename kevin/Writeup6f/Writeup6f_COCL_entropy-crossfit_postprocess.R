rm(list=ls())
load("../../../../out/kevin/Writeup6f/COCL_entropy-crossfit_tmp.RData")

zz <- sapply(result_list, function(x){
  if(is.list(x)) x$pvalue else NA
})
zz <- zz[!is.na(zz)]
zz

i <- which(names(result_list) == "MYL6")
par(mfrow = c(1,2))
plot(result_list[[i]]$grenander_die1$x, 
     result_list[[i]]$grenander_die1$pdf, main = "Loser1", xlab = "x", ylab = "pdf")
plot(result_list[[i]]$grenander_win1$x, 
     result_list[[i]]$grenander_win1$pdf, main = "Winner1", xlab = "x", ylab = "pdf")

par(mfrow = c(1,2))
plot(result_list[[i]]$grenander_die2$x, 
     result_list[[i]]$grenander_die2$pdf, main = "Loser2", xlab = "x", ylab = "pdf")
plot(result_list[[i]]$grenander_win2$x, 
     result_list[[i]]$grenander_win2$pdf, main = "Winner2", xlab = "x", ylab = "pdf")

result_list[[i]]$teststat
