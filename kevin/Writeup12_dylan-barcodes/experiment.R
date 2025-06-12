lin_mat,
cell_lower_limit = 10
cor_threshold = 0.55
warn_merging = TRUE
verbose = 4


stopifnot(inherits(lin_mat, "dgCMatrix"))
stopifnot(length(rownames(lin_mat)) > 0)

lin_mat_t <- Matrix::t(lin_mat)
nlineages <- nrow(lin_mat)
num_cells <- sapply(1:nlineages, function(j){
  length(.nonzero_col(lin_mat_t, col_idx = j, bool_value = F))
})

if(verbose > 0) print("Compute correlation matrix")
# discard any lineages in question (for the purposes of merging) that are too small
lin_mat_t <- lin_mat_t[,which(num_cells >= cell_lower_limit)]
# compute a correlation matrix (# rows/columsn = number of lineages)
cor_mat <- .custom_correlation(lin_mat_t)
cor_mat[lower.tri(cor_mat, diag = T)] <- NA
# determine all the lineages to merge. arr_idx is a 2-column matrix
arr_idx <- which(cor_mat >= cor_threshold, arr.ind = TRUE)
if(verbose > 2){
  print(paste0("There are ", nrow(arr_idx), " number of highly correlated lineages to resolve."))
}

# tabulate uniq_lineage, which is going to keep track of which lineage is 
#  assigned to which "cluster"
lineage_list <- vector("list", length = 0)
uniq_lineage <- sort(unique(as.numeric(arr_idx)))
uniq_lineage <- data.frame(Lineage = uniq_lineage, 
                           Cluster = rep(NA, length(uniq_lineage)))
lineage_name <- rownames(cor_mat)

if(verbose > 0) print("Determining how to merge lineages")
# determine the clusters of lineages to merge
for(i in 1:nrow(arr_idx)){
  if(verbose > 1 && verbose < 4 && nrow(arr_idx) > 10 && i %% floor(nrow(arr_idx)/10) == 0) cat('*')
  if(verbose == 4) print(paste0("Working on lineage ", i))
  
  # find the rows in the uniq_lineage table on which we're currently working on
  # this represents a pair of lineages
  lineage_idx <- which(uniq_lineage[,"Lineage"] %in% arr_idx[i,])
  
  if(length(lineage_list) == 0) {
    # if we have not yet merged any lineages, then this is straight-forward
    # create a new cluster
    lineage_list[[1]] <- sort(arr_idx[1,])
    names(lineage_list)[[1]] <- "c1"
    uniq_lineage[lineage_idx,"Cluster"] <- "c1"
    
  } else {
    # otherwise...
    
    if(all(is.na(uniq_lineage[lineage_idx,2]))) {
      # if this is a completely new cluster (i.e., all unassigned), also pretty straight-forward
      # create a new cluster
      
      new_list <- list(sort(arr_idx[i,])); names(new_list) <- paste0("c", length(lineage_list)+1)
      lineage_list <- c(lineage_list, new_list)
      uniq_lineage[lineage_idx,"Cluster"] <- paste0("c", length(lineage_list))
      if(verbose > 2 && length(lineage_list) %% 100 == 0) {
        print(paste0("There are currently ", length(lineage_list), " cluster of lineages"))
        
        if(verbose > 3) {
          print("The size of the lineages are: ")
          print(quantile(sapply(lineage_list, length)))
        }
      }
      
    } else {
      # the difficulty is this step. The new lineage in question has a high
      #  correlation with an existing cluster of lineages
      
      # first find all the lineages in this existing cluster
      val <- uniq_lineage[lineage_idx,"Cluster"]
      val <- val[!is.na(val)]
      val <- unique(val)
      if(length(val) > 1) {
        if(warn_merging) {warning("Merging happening")}
        
        # just merge the two lineages
        ## to do this, we need to remove the old lineage from lineage_list 
        stopifnot(length(val) == 2)
        val <- sort(val)
        val_pick <- val[1]
        val_dump <- val[2]
        lineage_list[[val_pick]] <- sort(unique(c(lineage_list[[val_pick]], lineage_list[[val_dump]], arr_idx[i,])))
        lineage_list[[val_dump]] <- NA
        
        ## we also need to update all the values in uniq_lineage
        uniq_lineage[lineage_idx,"Cluster"] <- val_pick
        changing_idx <- which(uniq_lineage[,"Cluster"] == val_dump)
        uniq_lineage[changing_idx,"Cluster"] <- val_pick
        
      } else {
        # just add this lineage to the existing cluster
        lineage_list[[val]] <- sort(unique(c(lineage_list[[val]], arr_idx[i,])))
        uniq_lineage[lineage_idx,"Cluster"] <- paste0("c", val)
      }
    }
  }
}

# rename all the indices with their lineage name
lineage_list_name <- lapply(lineage_list, function(vec){
  if(all(is.na(vec))) return(NA)
  lineage_name[vec]
})
uniq_lineage[,"Lineage"] <- lineage_name[uniq_lineage[,"Lineage"]]

######################

# compute the correlation values
min_cor_vec <- sapply(lineage_list_name, function(vec){
  if(all(is.na(vec))) return(NA)
  min(cor_mat[vec,vec], na.rm = TRUE)
})
quantile(min_cor_vec, na.rm = TRUE)
