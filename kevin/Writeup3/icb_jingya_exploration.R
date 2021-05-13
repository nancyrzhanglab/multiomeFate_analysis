rm(list=ls())

gmat <- readRDS("../../../../data/ICB_mouse/gmat_RNA_counts_ALL_dim35.rds")
dim(gmat)
gmat[1:5,1:5]
length(gmat@x)
length(gmat@x)/prod(dim(gmat))
quantile(gmat@x)

metadata <- readRDS("../../../../data/ICB_mouse/metadat.rds")
head(metadata)

amat <- readRDS("../../../../data/ICB_mouse/ALL_processed_FM.rds")
class(amat)
dim(amat)
amat[1:5,1:5]
length(amat@x)
length(amat@x)/prod(dim(amat))
quantile(amat@x)

