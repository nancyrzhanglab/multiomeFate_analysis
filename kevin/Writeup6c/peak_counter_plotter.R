peak_plotting_func <- function(active_gene_vec, 
                               fitness_vec2, 
                               survive_idx, die_idx, unknown_idx,
                               survive_idx2, die_idx2,
                               filename,
                               condition_name,
                               main4, xlab, ylab4){
  len1 <- length(active_gene_vec)
  len2 <- length(fitness_vec2)
  
  stopifnot(all(survive_idx <= len1), all(die_idx <= len1), all(unknown_idx <= len1),
            all(survive_idx2 <= len2), all(die_idx2 <= len2),
            length(survive_idx) == length(survive_idx2),
            length(die_idx) == length(die_idx2))
  
  vec_1 <- active_gene_vec[survive_idx]
  vec_2 <- active_gene_vec[die_idx]
  vec_3 <- active_gene_vec[c(die_idx, unknown_idx)]
  
  range_vec <- c(0, quantile(c(vec_1, vec_3), probs = 0.99))
  break_vec <- seq(min(range_vec), max(range_vec), length.out = 21)
  
  vec_1_tmp <- vec_1[vec_1 <= range_vec[2]]
  vec_2_tmp <- vec_2[vec_2 <= range_vec[2]]
  vec_3_tmp <- vec_3[vec_3 <= range_vec[2]]
  
  png(filename, height = 1000, width = 3000, units = "px", res = 500)
  par(mfrow = c(1,4), mar = c(4,4,4,0.5))
  hist(vec_1_tmp, breaks = break_vec, col = "gray", 
       main = paste0(condition_name, ", surviving cells\n(w/ lineage),\n",
                     "Mean: ", round(mean(vec_1)), ", n=", length(vec_1)),
       ylab = "Frequency", xlab = xlab,
       cex.main = 0.6, cex.lab = 0.6)
  mean_val <- mean(vec_1)
  median_val <- median(vec_1)
  lines(rep(mean_val, 2), c(-1e6,1e6), col = 2, lwd = 2)
  lines(rep(median_val, 2), c(-1e6,1e6), col = 3, lwd = 2, lty = 2)
  
  ###################
  pval <- stats::wilcox.test(x = vec_1, y = vec_2)
  hist(vec_2_tmp, breaks = break_vec, col = "gray", 
       main = paste0(condition_name, ", dying cells\n(w/ lineage),\n",
                     "Mean: ", round(mean(vec_2)), ", n=", length(vec_2), ",\n-Log10 pval: ", round(-log10(pval$p.value),3)),
       ylab = "Frequency", xlab = xlab,
       cex.main = 0.6, cex.lab = 0.6)
  mean_val <- mean(vec_2)
  median_val <- median(vec_2)
  lines(rep(mean_val, 2), c(-1e6,1e6), col = 2, lwd = 2)
  lines(rep(median_val, 2), c(-1e6,1e6), col = 3, lwd = 2, lty = 2)
  
  pval <- stats::wilcox.test(x = vec_1, y = vec_3)
  hist(vec_3_tmp, breaks = break_vec, col = "gray", 
       main = paste0(condition_name, ", dying cells\n(w/ and w/o lineage),\n",
                     "Mean: ", round(mean(vec_3)), ", n=", length(vec_3), ",\n-Log10 pval: ", round(-log10(pval$p.value),3)),
       ylab = "Frequency", xlab = xlab,
       cex.main = 0.6, cex.lab = 0.6)
  mean_val <- mean(vec_3)
  median_val <- median(vec_3)
  lines(rep(mean_val, 2), c(-1e6,1e6), col = 2, lwd = 2)
  lines(rep(median_val, 2), c(-1e6,1e6), col = 3, lwd = 2, lty = 2)
  
  tmp_x <- fitness_vec2[c(survive_idx2, die_idx2)]
  tmp_y <- active_gene_vec[c(survive_idx, die_idx)]
  df <- data.frame(x = tmp_x, y = tmp_y)
  lm_res <- stats::lm(y ~ x, data = df)
  cor_val <- stats::cor(tmp_x, tmp_y)
  cor_pval <- stats::cor.test(tmp_x, tmp_y)
  
  plot(x = custom_jitter(fitness_vec2[c(survive_idx2, die_idx2)]), 
       y = active_gene_vec[c(survive_idx, die_idx)],
       main = paste0("Regressing ", main4, "\nonto fitness, Corr: ", round(cor_val,2),"\n",
                     "-Log10 pval of corr: ", round(-log10(cor_pval$p.value), 3)),
       ylab = ylab4, xlab = "Fitness of day0 cell in day10",
       cex.main = 0.6, cex.lab = 0.6,
       pch = 16, col = rgb(0.5, 0.5, 0.5, 0.2))
  lines(rep(log1p(20), 2), c(-1e6,1e6), col = 3, lty = 2)
  abline(a = stats::coef(lm_res)[1], b = stats::coef(lm_res)[2],
         col = 2, lwd = 2)
  graphics.off()
  
  invisible()
}

custom_jitter <- function(vec, threshold = 1.2, width = 0.15){
  idx <- which(vec <= threshold)
  if(length(idx) > 0){
    vec[idx] <- vec[idx] + stats::runif(length(idx), min = 0, max = width)
  }
  
  vec
}