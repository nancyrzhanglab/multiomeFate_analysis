rm(list=ls())
set.seed(10)

treatment_vec <- c("CIS", "COCL2", "DABTRAM")
day_early_vec <- c("day0", "day10")

for(treatment in treatment_vec){
  for(ii in 1:2){
    day_early <- day_early_vec[ii]
    
    print(paste0("Working on ", day_early, " for ", treatment))
    load(paste0("../../../../out/kevin/Writeup6r/Writeup6r_", treatment, "_", day_early, "_lineage-imputation_postprocess.RData"))
    
    rna_idx <- grep("^fastTopic*", colnames(cell_features))
    atac_idx <- grep("^peakVI*", colnames(cell_features))
    
    mat_1 <- cell_features[,rna_idx,drop = F]
    mat_2 <- cell_features[,atac_idx,drop = F]
    n <- nrow(mat_1)
    p1 <- ncol(mat_1); p2 <- ncol(mat_2)
    p <- min(p1, p2)
    num_metacells <- max(round(sqrt(n)), 200)
    
    set.seed(10)
    multiSVD_obj <- tiltedCCA::create_multiSVD(mat_1 = mat_1, mat_2 = mat_2,
                                               dims_1 = 1:p1, dims_2 = 1:p2,
                                               center_1 = T, center_2 = T,
                                               normalize_row = T,
                                               normalize_singular_value = T,
                                               recenter_1 = F, recenter_2 = F,
                                               rescale_1 = F, rescale_2 = F,
                                               scale_1 = T, scale_2 = T,
                                               verbose = 0)
    multiSVD_obj <- tiltedCCA::form_metacells(input_obj = multiSVD_obj,
                                              large_clustering_1 = NULL, 
                                              large_clustering_2 = NULL, 
                                              num_metacells = num_metacells,
                                              verbose = 0)
    multiSVD_obj <- tiltedCCA::compute_snns(input_obj = multiSVD_obj,
                                            latent_k = min(10, round(p/2)),
                                            num_neigh = max(num_metacells/10, 15),
                                            bool_cosine = T,
                                            bool_intersect = F,
                                            min_deg = 2,
                                            verbose = 0)
    multiSVD_obj <- tiltedCCA::tiltedCCA(input_obj = multiSVD_obj,
                                         verbose = 0)
    multiSVD_obj <- tiltedCCA::fine_tuning(input_obj = multiSVD_obj,
                                           verbose = 0)
    multiSVD_obj <- tiltedCCA::tiltedCCA_decomposition(input_obj = multiSVD_obj,
                                                       verbose = 0,
                                                       bool_modality_1_full = T,
                                                       bool_modality_2_full = T)
    
    date_of_run <- Sys.time()
    session_info <- devtools::session_info()
    save(multiSVD_obj, 
         date_of_run, session_info,
         file = paste0("../../../../out/kevin/Writeup6s/Writeup6s_", treatment, "_", day_early, "_tiltedCCA.RData"))
    
  }
}