# see https://cran.r-project.org/web/packages/RcppAnnoy/RcppAnnoy.pdf

rm(list=ls())
p <- 10; n <- 100
a <- new(RcppAnnoy::AnnoyEuclidean, p)
dat <- sapply(1:n, function(x){runif(p)})
for (i in 1:n) a$addItem(i-1, dat[,i]) # Annoy uses zero indexing
a$getNItems()
sum(abs(dat[,1] - a$getItemsVector(0)))
a$build(50) # you need to build the tree

a$getNNsByItem(0, 5)
a$getNNsByItemList(0, 5, -1, FALSE) # the same

v <- runif(p)
a$getNNsByVector(v, 5)
