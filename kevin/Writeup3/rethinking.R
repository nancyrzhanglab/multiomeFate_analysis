rm(list=ls())
len <- 100
x <- sin(seq(0, 0.75*pi, length.out = len))
y <- sapply(2:len, function(i){
  3*x[i]+50*(x[i]-x[i-1])
})
png("../../out/fig/briefnote1/petal.png", height = 1000, width = 2000, res = 300, units = "px")
par(mfrow = c(1,2), mar = c(4,4,4,0.5))
plot(x[-1], 3*x[-1], asp = T, pch = 16,
     col = colorRampPalette(c("red", "blue"))(len-1),
     xlab = "RNA", ylab = "ATAC", main = "Only linear")
plot(x[-1], y-min(y), asp = T, pch = 16,
     col = colorRampPalette(c("red", "blue"))(len-1),
     xlab = "RNA", ylab = "ATAC", main = "With previous ATAC")
graphics.off()

###################

rm(list=ls())
len <- 100; seq_vec <- seq(-5, 5, length.out = len)
mat_sigmoid <- sapply(1:len, function(i){
  1/(1+exp(-(seq_vec - seq_vec[i])))
})

generate_graph <- function(prob_mat, ...){
  mat <- matrix(stats::rbinom(prod(dim(prob_mat)), size = 1, 
                              prob = prob_mat), 
                nrow = nrow(prob_mat), ncol = ncol(prob_mat))
  print(sum(mat))
  k_vec <- 2:10; counter <- 1
  while(TRUE){
    rann_res <- RANN::nn2(mat, k = k_vec[counter])
    n <- nrow(rann_res$nn.idx)
    adj_mat <- sapply(1:n, function(i){
      vec <- rep(0, n)
      vec[rann_res$nn.idx[i,]] <- 1
      vec
    })
    
    g <- igraph::graph_from_adjacency_matrix(adj_mat, mode = "undirected")
    g <- igraph::simplify(g)
    
    res <- igraph::components(g)
    if(res$no == 1) break()
    counter <- counter+1
  }
  
  n <- igraph::vcount(g)
  igraph::V(g)$color <- colorRampPalette(c("red", "blue"))(n)
  graphics::plot(g, vertex.label=NA, ...)
}

png("../../out/fig/briefnote1/sigmoid.png", height = 900, width = 3200, res = 300, units = "px")
par(mfrow = c(1,4), mar = c(4,4,4,0.5))
image(.rotate(mat_sigmoid), asp = T, ylab = "Different cells", xlab = "Across genome",
      main = "Blueprint")
par(mar = c(0.5,0.5,4,0.5))
for(i in 1:3){
  set.seed(i)
  generate_graph(mat_sigmoid, main = paste0("Simulation ", i))
}
graphics.off()

###########################

# # gaussian density
# 
# len <- 100; seq_vec <- seq(-1, 1, length.out = len)
# mat_bump <- sapply(1:len, function(i){
#   vec <- rep(0, length(seq_vec))
#   val <- seq_vec + seq_vec[i]
#   idx <- intersect(which(val > -1), which(val < 1))
#   if(length(idx) > 0) vec[idx] <- exp(-1/(1-val[idx]^2))
#   vec
# })
# mat_bump <- mat_bump-min(mat_bump)/(diff(range(mat_bump)))
# plot(mat_bump[50,])
# image(.rotate(mat_bump), asp = T)
# 
# par(mfrow = c(1,3), mar = c(0.5,0.5,0.5,0.5))
# for(i in 1:3){
#   set.seed(i)
#   generate_graph(mat_bump)
# }

##########################

len <- 100; seq_vec <- seq(0, 1.5*pi, length.out = len)
mat_sin <- sapply(1:len, function(i){
  sin(seq_vec + seq_vec[i])
})
mat_sin <- (mat_sin+1)/2; range(mat_sin)
plot(mat_sin[50,])
image(.rotate(mat_sin), asp = T)

png("../../out/fig/briefnote1/sin.png", height = 900, width = 3200, res = 300, units = "px")
par(mfrow = c(1,4), mar = c(4,4,4,0.5))
image(.rotate(mat_sin), asp = T, ylab = "Different cells", xlab = "Across genome",
      main = "Blueprint")
for(i in 1:3){
  set.seed(i)
  generate_graph(mat_sin, main = paste0("Simulation ", i))
}
graphics.off()

############################

# low values
len <- 100; seq_vec <- seq(0, 1.5*pi, length.out = len)
mat_sin <- sapply(1:len, function(i){
  sin(seq_vec + seq_vec[i])
})
mat_sin <- (mat_sin-min(mat_sin))/(diff(range(mat_sin)))
mat_sin <- 0.1*mat_sin

png("../../out/fig/briefnote1/sin_low.png", height = 900, width = 3200, res = 300, units = "px")
par(mfrow = c(1,4), mar = c(4,4,4,0.5))
image(.rotate(mat_sin), zlim = c(0,1), asp = T, ylab = "Different cells", xlab = "Across genome",
      main = "Blueprint")
for(i in 1:3){
  set.seed(i)
  generate_graph(mat_sin, main = paste0("Simulation ", i))
}
graphics.off()


# cyclical values
len <- 100; seq_vec <- seq(0, 4*pi, length.out = len)
mat_sin <- sapply(1:len, function(i){
  sin(seq_vec + seq_vec[i])
})
mat_sin <- (mat_sin-min(mat_sin))/(diff(range(mat_sin)))

png("../../out/fig/briefnote1/sin_cylical.png", height = 900, width = 3200, res = 300, units = "px")
par(mfrow = c(1,4), mar = c(4,4,4,0.5))
image(.rotate(mat_sin), zlim = c(0,1), asp = T, ylab = "Different cells", xlab = "Across genome",
      main = "Blueprint")
for(i in 1:3){
  set.seed(i)
  generate_graph(mat_sin, main = paste0("Simulation ", i))
}
graphics.off()

# not enough variables
len <- 100; seq_vec <- seq(0, 4*pi, length.out = len)
mat_sin <- sapply(1:len, function(i){
  sin(seq_vec + seq_vec[i])
})
mat_sin <- (mat_sin-min(mat_sin))/(diff(range(mat_sin)))
mat_sin <- mat_sin[,1:5]

png("../../out/fig/briefnote1/sin_notenough.png", height = 900, width = 3200, res = 300, units = "px")
par(mfrow = c(1,4), mar = c(4,4,4,0.5))
image(.rotate(mat_sin), zlim = c(0,1), asp = 100/5, ylab = "Different cells", xlab = "Across genome",
      main = "Blueprint")
for(i in 1:3){
  set.seed(i)
  generate_graph(mat_sin, main = paste0("Simulation ", i))
}
graphics.off()





