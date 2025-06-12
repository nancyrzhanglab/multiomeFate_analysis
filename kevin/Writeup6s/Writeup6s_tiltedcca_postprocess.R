rm(list=ls())
set.seed(10)

treatment_vec <- c("CIS", "COCL2", "DABTRAM")
day_early_vec <- c("day0", "day10")
day_later_vec <- c("day10", "week5")

percentage_list <- vector("list", 0)

for(treatment in treatment_vec){
  for(day_early in day_early_vec){
    print(paste0(treatment, "_", day_early))
    
    load(paste0("../../../../out/kevin/Writeup6s/Writeup6s_", treatment, "_", day_early, "_tiltedCCA.RData"))
    
    y <- scale(cell_imputed_score)
    basis_list <- list(
      common = multiSVD_obj$tcca_obj$common_score,
      distinct_rna = multiSVD_obj$tcca_obj$distinct_score_1,
      distinct_atac = multiSVD_obj$tcca_obj$distinct_score_2,
      common_and_rna = cbind(multiSVD_obj$tcca_obj$common_score, 
                             multiSVD_obj$tcca_obj$distinct_score_1),
      common_and_atac = cbind(multiSVD_obj$tcca_obj$common_score, 
                              multiSVD_obj$tcca_obj$distinct_score_2)
    )
    
    percentage_vec <- sapply(1:length(basis_list), function(ii){
      df <- as.data.frame(cbind(y, basis_list[[ii]]))
      colnames(df)[1] <- "y"
      
      lm_res <- stats::lm(y ~ . - 1, data = df)
      summary(lm_res)$r.squared
    })
    names(percentage_vec) <- names(basis_list)
    
    percentage_list[[paste0(treatment, "_", day_early)]] <- percentage_vec
  }
}

###################

rsquare_list <- vector("list", 0)

for(treatment in treatment_vec){
  for(ii in 1:2){
    day_early <- day_early_vec[ii]
    day_later <- day_later_vec[ii]
    
    print(paste0(treatment, "_", day_early))
    
    load(paste0("../../../../out/kevin/Writeup6r/Writeup6r_", treatment, "_", day_early, "_lineage-imputation_postprocess.RData"))
    
    rna_idx <- grep("^fastTopic*", colnames(cell_features))
    atac_idx <- grep("^peakVI*", colnames(cell_features))
    
    mat_1 <- cell_features[,rna_idx,drop = F]
    mat_2 <- cell_features[,atac_idx,drop = F]
    
    df <- as.data.frame(cbind(scale(cell_imputed_score), mat_1))
    colnames(df)[1] <- "y"
    lm_res <- stats::lm(y ~ . - 1, data = df)
    rna_rsq <- summary(lm_res)$r.squared
    
    df <- as.data.frame(cbind(scale(cell_imputed_score), mat_2))
    colnames(df)[1] <- "y"
    lm_res <- stats::lm(y ~ . - 1, data = df)
    atac_rsq <- summary(lm_res)$r.squared
    
    r2_vec <- c(
      common = as.numeric(percentage_list[[paste0(treatment, "_", day_early)]]["common"]),
      rna = rna_rsq,
      atac = atac_rsq
    )
    
    rsquare_list[[paste0(treatment, "_", day_early)]] <- r2_vec
    
    df <- data.frame(modality = names(r2_vec),
                     r2 = r2_vec)

    p1 <- ggplot2::ggplot(df, ggplot2::aes(x=modality, y=r2)) +
      ggplot2::geom_bar(stat = "identity")
    p1 <- p1 + ggplot2::ggtitle(paste0(
      treatment, "\n", day_later, " growth potential of ", day_early,
      " cells\nImportance of each modality")
    )
    p1 <- p1 + ggplot2::scale_x_discrete(limits = c("common", "rna", "atac")) + ggplot2::ylim(0, 1)
    ggplot2::ggsave(filename = paste0("../../../../out/figures/kevin/Writeup6s/Writeup6s_",
                                      treatment, "-", day_early, "_modality-weights.png"),
                    p1, device = "png", width = 3.5, height = 3.5, units = "in")
  }
}

