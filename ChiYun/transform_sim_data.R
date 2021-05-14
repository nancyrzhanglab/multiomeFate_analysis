## transform the simulated data
# min value set to 0 from the simulateAtacRna.r script
dim(Y)
Y2=Y
for(ii in 1:nrow(Y)){
  if(min(Y[ii,])<0){
    Y2[ii,]=Y2[ii,]-min(Y[ii,])
  }
}
summary(as.numeric(Y2))
dim(Y2)

dat2=list()
dat2$df_x=data.frame(name=paste0("peak", 1:nrow(Xflat)), location=rep(((1:nrow(Y))*5000-2500), rep(5, nrow(Y))), stringsAsFactors = F)
dat2$df_y=data.frame(name=paste0("gene", 1:nrow(Y)), location=((1:nrow(Y))*5000-2500),baseline=0, stringsAsFactors = F)
dat2$obs_x=t(Xflat)
dat2$obs_y=t(Y2)

dat2$true_x=t(Xflat)
dat2$true_y=t(Y2)



dat2$df_info$time=trueTime/1000
dat2$df_info$counter=trueTime/1000
dat2$df_info$branch=trueBranch


