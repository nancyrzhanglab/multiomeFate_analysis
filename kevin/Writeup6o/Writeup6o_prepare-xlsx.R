rm(list=ls())
library(openxlsx)
# see https://cran.r-project.org/web/packages/openxlsx/vignettes/Introduction.html

ttest_gene_prefix <- "~/project/Multiome_fate/out/emilia/task_identify_genes_corr_growth_and_lineage_specific/"
ttest_gene_treatment <- c("day10_CIS", "day10_COCL2", "day10_DABTRAM")
ttest_gene_suffix <- "_adaptation_genes_diff_t_tests_on_day0_cells.csv"

# load cis_cor_vec, cocl2_cor_vec, dabtram_cor_vec
# prefix_corr <- "~/project/Multiome_fate/out/kevin/Writeup6n/Writeup6n_correlation-with-growthpotential_"
# timepoint_vec_out <- c("day0.RData", "day10.RData")
corr_gene_file <- "~/project/Multiome_fate/out/kevin/Writeup6n/Writeup6n_lineage-imputation_day0-day10-export.RData"

anova_gene_prefix <- "~/nzhanglab/project/emiliac/lineage_trace/task0_explore_lineage_variability/outs/ANOVA/"
anova_gene_treatment <- c("day0", "day10_CIS", "day10_COCL2", "day10_DABTRAM")
anova_gene_suffix <- "_processed_RNA_pvals.csv"

adaptation_gene_prefix <- "~/project/Multiome_fate/out/emilia/task_identify_genes_corr_growth_and_lineage_specific/"
adaptation_gene_treatment <- c("day10_CIS", "day10_COCL2", "day10_DABTRAM")
adaptation_gene_suffix <- "_adaptation_genes.csv"

file_list <- list(
  day0_CIS = list(
    ttest_file = paste0(ttest_gene_prefix, "day10_CIS", ttest_gene_suffix),
    correlation_file = corr_gene_file,
    correlation_name = "day0_CIS",
    anova_file = paste0(anova_gene_prefix, "day0", anova_gene_suffix),
    adaptation_file = paste0(adaptation_gene_prefix, adaptation_gene_treatment[1], adaptation_gene_suffix)
  ),
  day0_COCL2 = list(
    ttest_file = paste0(ttest_gene_prefix, "day10_COCL2", ttest_gene_suffix),
    correlation_file = corr_gene_file,
    correlation_name = "day0_COCL2",
    anova_file = paste0(anova_gene_prefix, "day0", anova_gene_suffix),
    adaptation_file = paste0(adaptation_gene_prefix, adaptation_gene_treatment[2], adaptation_gene_suffix)
  ),
  day0_DABTRAM = list(
    ttest_file = paste0(ttest_gene_prefix, "day10_DABTRAM", ttest_gene_suffix),
    correlation_file = corr_gene_file,
    correlation_name = "day0_DABTRAM",
    anova_file = paste0(anova_gene_prefix, "day0", anova_gene_suffix),
    adaptation_file = paste0(adaptation_gene_prefix, adaptation_gene_treatment[3], adaptation_gene_suffix)
  ),
  day10_CIS = list(
    ttest_file = NA,
    correlation_file = corr_gene_file,
    correlation_name = "day10_CIS",
    anova_file = paste0(anova_gene_prefix, "day10_CIS", anova_gene_suffix),
    adaptation_file = paste0(adaptation_gene_prefix, adaptation_gene_treatment[1], adaptation_gene_suffix)
  ),
  day10_COCL2 = list(
    ttest_file = NA,
    correlation_file = corr_gene_file,
    correlation_name = "day10_COCL2",
    anova_file = paste0(anova_gene_prefix, "day10_COCL2", anova_gene_suffix),
    adaptation_file = paste0(adaptation_gene_prefix, adaptation_gene_treatment[2], adaptation_gene_suffix)
  ),
  day10_DABTRAM = list(
    ttest_file = NA,
    correlation_file = corr_gene_file,
    correlation_name = "day10_DABTRAM",
    anova_file = paste0(anova_gene_prefix, "day10_DABTRAM", anova_gene_suffix),
    adaptation_file = paste0(adaptation_gene_prefix, adaptation_gene_treatment[3], adaptation_gene_suffix)
  )
)

wb <- openxlsx::createWorkbook()
for(sheet_name in names(file_list)){
  print(sheet_name)
  
  lis <- file_list[[sheet_name]]
  # load cis_cor_vec, cocl2_cor_vec, dabtram_cor_vec
  load(lis$correlation_file)
  cor_mat <- correlation_list[[lis$correlation_name]]
  cor_mat <- as.data.frame(cor_mat)
  
  anova_df <- read.csv(lis$anova_file)
  rownames(anova_df) <- anova_df$feature
  
  # gene names are based on anova and correlation
  gene_vec <- sort(unique(c(rownames(anova_df), rownames(cor_mat))))
  
  # add values to cor_mat as needed
  genes_to_add <- setdiff(gene_vec,  rownames(cor_mat))
  if(length(genes_to_add) > 0){
    mat <- data.frame(correlation = NA,
                      p.value = NA)
    rownames(mat) <- genes_to_add
    cor_mat <- rbind(cor_mat, mat)
  }
  cor_mat <- cor_mat[gene_vec,]
  
  # add values to anova_df as needed
  genes_to_add <- setdiff(gene_vec,  rownames(anova_df))
  if(length(genes_to_add) > 0){
    mat <- data.frame(feature = genes_to_add,
                      F_val = NA,
                      p_val = NA)
    rownames(mat) <- genes_to_add
    anova_df <- rbind(anova_df, mat)
  }
  anova_df <- anova_df[gene_vec,]
  
  if(all(is.na(lis$ttest_file))) {
    ttest_df <- data.frame(gene = anova_df$feature,
                           winner_minus_others = NA,
                           t_statistic = NA,
                           p_val = NA,
                           p_val_adjust = NA)
  } else {
    ttest_df <- read.csv(lis$ttest_file)
    ttest_df <- ttest_df[,c("gene", 
                            "winner_minus_others", 
                            "t_statistic",
                            "p_val",
                            "p_val_adjust")]
  }
  rownames(ttest_df) <- ttest_df$gene
  
  # add values to ttest as needed
  genes_to_add <- setdiff(gene_vec,  rownames(ttest_df))
  if(length(genes_to_add) > 0){
    mat <- data.frame(gene = genes_to_add,
                      winner_minus_others = NA,
                      t_statistic = NA,
                      p_val = NA,
                      p_val_adjust = NA)
    rownames(mat) <- genes_to_add
    ttest_df <- rbind(ttest_df, mat)
  }
  ttest_df <- ttest_df[gene_vec,]
  
  if(!all(is.na(lis$adaptation_file))){
    adapation_df <- read.csv(lis$adaptation_file)
    adaptation_vec <- rep(FALSE, length(gene_vec))
    names(adaptation_vec) <- gene_vec
    adaptation_vec[adapation_df$gene] <- TRUE
  } else {
    adaptation_vec <- rep(NA, length(gene_vec))
    names(adaptation_vec) <- gene_vec
  }
  
  df <- data.frame(gene = gene_vec,
                   cor_statistic = cor_mat[,"correlation"],
                   cor_pvalue = cor_mat[,"p.value"],
                   cor_pvalue_adj = stats::p.adjust(cor_mat[,"p.value"], method = "BH"),
                   ttest_difference = ttest_df$winner_minus_others,
                   ttest_tstat = ttest_df$t_statistic,
                   ttest_pvalue = ttest_df$p_val,
                   ttest_pvalue_adj = ttest_df$p_val_adjust,
                   anova_statistc = anova_df$F_val,
                   anova_pvalue = anova_df$p_val,
                   anova_pvalue_adj = stats::p.adjust(anova_df$p_val, method = "BH"),
                   adaptation_gene = adaptation_vec)
  rownames(df) <- gene_vec
  
  df <- df[order(abs(df$cor_statistic), decreasing = T),]
  if(!all(is.na(df$adaptation_gene))){
    df <- df[c(which(df$adaptation_gene == TRUE), which(df$adaptation_gene == FALSE)),]
  }
  na_vec <- is.na(df$anova_statistc)
  if(any(na_vec == FALSE)){
    df <- df[c(which(na_vec == FALSE), which(na_vec == TRUE)),]
  }
  
  openxlsx::addWorksheet(wb = wb, sheetName = sheet_name)
  writeData(wb = wb, sheet = sheet_name, x = df)
}
openxlsx::saveWorkbook(wb, 
                       file = "../../../../out/kevin/Writeup6o/Writeup6o_gene-expression_statistics.xlsx", 
                       overwrite = TRUE)


##################
# now for Chromvar things
corr_chromvar_prefix <- "~/nzhanglab/project/emiliac/lineage_trace/task_identify_features_predictive_growth_and_lineage_specific/outs/"
corr_chromvar_treatment <- c("day0_chromVar_day0_growth_potential_for_day10_correlation.RData",
              "day10_chromVar_day10_growth_potential_for_week5_correlation.RData")

anova_chromvar_prefix <- "~/nzhanglab/project/emiliac/lineage_trace/task_identify_features_predictive_growth_and_lineage_specific/outs/identify_lineage_specificity_ANOVA/"
anova_chromvar_treatment <- c("day0", "day10_CIS", "day10_COCL2", "day10_DABTRAM")
anova_chromvar_suffix <- "_ANOVA_ChromVar_pvals.csv"

adaptation_chromvar_prefix <- "~/project/Multiome_fate/out/emilia/task_identify_genes_corr_growth_and_lineage_specific/"
adaptation_chromvar_treatment <- c("day10_CIS", "day10_COCL2", "day10_DABTRAM")
adaptation_chromvar_suffix <- "_adaptation_motifs.csv"

file_list <- list(
  day0_CIS = list(
    ttest_file = NA,
    correlation_file = paste0(corr_chromvar_prefix, corr_chromvar_treatment[1]),
    correlation_name = "cis_cor_vec",
    anova_file = paste0(anova_chromvar_prefix, anova_chromvar_treatment[1], anova_chromvar_suffix),
    adaptation_file = paste0(adaptation_chromvar_prefix, adaptation_chromvar_treatment[1], adaptation_chromvar_suffix)
  ),
  day0_COCL2 = list(
    ttest_file = NA,
    correlation_file = paste0(corr_chromvar_prefix, corr_chromvar_treatment[1]),
    correlation_name = "cocl2_cor_vec",
    anova_file = paste0(anova_chromvar_prefix, anova_chromvar_treatment[1], anova_chromvar_suffix),
    adaptation_file = paste0(adaptation_chromvar_prefix, adaptation_chromvar_treatment[2], adaptation_chromvar_suffix)
  ),
  day0_DABTRAM = list(
    ttest_file = NA,
    correlation_file = paste0(corr_chromvar_prefix, corr_chromvar_treatment[1]),
    correlation_name = "dabtram_cor_vec",
    anova_file = paste0(anova_chromvar_prefix, anova_chromvar_treatment[1], anova_chromvar_suffix),
    adaptation_file = paste0(adaptation_chromvar_prefix, adaptation_chromvar_treatment[3], adaptation_chromvar_suffix)
  ),
  day10_CIS = list(
    ttest_file = NA,
    correlation_file = paste0(corr_chromvar_prefix, corr_chromvar_treatment[2]),
    correlation_name = "cis_cor_vec",
    anova_file = paste0(anova_chromvar_prefix, anova_chromvar_treatment[2], anova_chromvar_suffix),
    adaptation_file = paste0(adaptation_chromvar_prefix, adaptation_chromvar_treatment[1], adaptation_chromvar_suffix)
  ),
  day10_COCL2 = list(
    ttest_file = NA,
    correlation_file = paste0(corr_chromvar_prefix, corr_chromvar_treatment[2]),
    correlation_name = "cocl2_cor_vec",
    anova_file = paste0(anova_chromvar_prefix, anova_chromvar_treatment[3], anova_chromvar_suffix),
    adaptation_file = paste0(adaptation_chromvar_prefix, adaptation_chromvar_treatment[2], adaptation_chromvar_suffix)
  ),
  day10_DABTRAM = list(
    ttest_file = NA,
    correlation_file = paste0(corr_chromvar_prefix, corr_chromvar_treatment[2]),
    correlation_name = "dabtram_cor_vec",
    anova_file = paste0(anova_chromvar_prefix, anova_chromvar_treatment[4], anova_chromvar_suffix),
    adaptation_file = paste0(adaptation_chromvar_prefix, adaptation_chromvar_treatment[3], adaptation_chromvar_suffix)
  )
)

wb <- openxlsx::createWorkbook()
for(sheet_name in names(file_list)){
  print(sheet_name)
  
  lis <- file_list[[sheet_name]]
  load(lis$correlation_file)
  if(lis$correlation_name == "cis_cor_vec"){
    cor_mat <- cis_cor_vec
  } else if(lis$correlation_name == "cocl2_cor_vec"){
    cor_mat <- cocl2_cor_vec
  } else {
    cor_mat <- dabtram_cor_vec
  }
  
  anova_df <- read.csv(lis$anova_file)
  rownames(anova_df) <- anova_df$motif_names
  
  # gene names are based on anova and correlation
  motif_vec <- sort(unique(c(rownames(anova_df), rownames(cor_mat))))
  
  # add values to cor_mat as needed
  motifs_to_add <- setdiff(motif_vec,  rownames(cor_mat))
  if(length(motifs_to_add) > 0){
    mat <- data.frame(correlation = NA,
                      p.value = NA)
    rownames(mat) <- genes_to_add
    cor_mat <- rbind(cor_mat, mat)
  }
  cor_mat <- cor_mat[motif_vec,]
  
  # add values to anova_df as needed
  motifs_to_add <- setdiff(motif_vec,  rownames(anova_df))
  if(length(genes_to_add) > 0){
    mat <- data.frame(feature = motifs_to_add,
                      F_val = NA,
                      p_val = NA)
    rownames(mat) <- genes_to_add
    anova_df <- rbind(anova_df, mat)
  }
  anova_df <- anova_df[motif_vec,]
  
  if(all(is.na(lis$ttest_file))) {
    ttest_df <- data.frame(motif = anova_df$feature,
                           winner_minus_others = NA,
                           t_statistic = NA,
                           p_val = NA,
                           p_val_adjust = NA)
  } else {
    ttest_df <- read.csv(lis$ttest_file)
    ttest_df <- ttest_df[,c("motif", 
                            "winner_minus_others", 
                            "t_statistic",
                            "p_val",
                            "p_val_adjust")]
  }
  rownames(ttest_df) <- ttest_df$motif
  
  # add values to ttest as needed
  motifs_to_add <- setdiff(motif_vec,  rownames(ttest_df))
  if(length(motifs_to_add) > 0){
    mat <- data.frame(motif = motifs_to_add,
                      winner_minus_others = NA,
                      t_statistic = NA,
                      p_val = NA,
                      p_val_adjust = NA)
    rownames(mat) <- mat$motif
    ttest_df <- rbind(ttest_df, mat)
  }
  ttest_df <- ttest_df[motif_vec,]
  
  if(!all(is.na(lis$adaptation_file))){
    adapation_df <- read.csv(lis$adaptation_file)
    adaptation_vec <- rep(FALSE, length(motif_vec))
    names(adaptation_vec) <- motif_vec
    adaptation_vec[adapation_df$motif_names] <- TRUE
  } else {
    adaptation_vec <- rep(NA, length(motif_vec))
    names(adaptation_vec) <- motif_vec
  }
  
  df <- data.frame(motif = motif_vec,
                   cor_statistic = cor_mat[,"correlation"],
                   cor_pvalue = cor_mat[,"p.value"],
                   cor_pvalue_adj = stats::p.adjust(cor_mat[,"p.value"], method = "BH"),
                   ttest_difference = ttest_df$winner_minus_others,
                   ttest_tstat = ttest_df$t_statistic,
                   ttest_pvalue = ttest_df$p_val,
                   ttest_pvalue_adj = ttest_df$p_val_adjust,
                   anova_statistc = anova_df$F_val,
                   anova_pvalue = anova_df$p_val,
                   anova_pvalue_adj = stats::p.adjust(anova_df$p_val, method = "BH"),
                   adaptation_motif = adaptation_vec)
  rownames(df) <- motif_vec
  
  df <- df[order(abs(df$cor_statistic), decreasing = T),]
  if(!all(is.na(df$adaptation_motif))){
    df <- df[c(which(df$adaptation_motif == TRUE), which(df$adaptation_motif == FALSE)),]
  }
  na_vec <- is.na(df$anova_statistc)
  if(any(na_vec == FALSE)){
    df <- df[c(which(na_vec == FALSE), which(na_vec == TRUE)),]
  }
  
  openxlsx::addWorksheet(wb = wb, sheetName = sheet_name)
  writeData(wb = wb, sheet = sheet_name, x = df)
}
openxlsx::saveWorkbook(wb, 
                       file = "../../../../out/kevin/Writeup6o/Writeup6o_motif-expression_statistics.xlsx", 
                       overwrite = TRUE)
