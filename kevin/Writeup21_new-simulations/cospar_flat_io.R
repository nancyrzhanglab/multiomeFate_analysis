# cospar_flat_io.R
# Kevin Z. Lin, 2026-09-22
#
# Two functions that carry a lineage-barcoded dataset from R to the CoSPAR
# Python runner (`run_cospar.py`) and bring its results back. The hand-off is a
# directory of flat files, never an in-memory bridge. Copy this file next to
# the analysis script and `source()` it; a reader without the skill installed
# still has everything they need.
#
# Flat contract written by `export_for_cospar()` into `export_dir`:
#   counts.mtx   MatrixMarket, rows = cells, columns = genes, non-negative
#   cells.txt    one cell ID per line, the row order of counts.mtx
#   genes.txt    one gene name per line, the column order of counts.mtx
#   obs.csv      cell_id, time_info, clone_id, state_info (one row per cell,
#                same order as cells.txt; clone_id may be NA)
#   X_pca.csv    optional, rows = cells in the same order, columns = dims
#   manifest.json  sizes, labels, and whether X_pca.csv was written
#
# `run_cospar.py` writes its results under `<export_dir>/cospar_out/`, which
# `import_cospar_results()` reads.

library(Matrix)
library(jsonlite)

#' Export a dataset for CoSPAR
#'
#' @param count_mat numeric matrix or `dgCMatrix`, rows are cells, columns are
#'   genes, non-negative and not log-transformed. Row and column names are
#'   required.
#' @param time_info character vector, length `nrow(count_mat)`, the time point
#'   label of every cell. Exactly the two labels passed to `run_cospar.py`.
#' @param clone_id character vector, length `nrow(count_mat)`, the clone of
#'   every cell, `NA` for an unbarcoded cell. A clone must appear at both time
#'   points to constrain the transition map.
#' @param state_info character vector, length `nrow(count_mat)`, the "fate"
#'   vocabulary CoSPAR will use. The two later-time-point fates named to
#'   `run_cospar.py` must be values of this vector.
#' @param export_dir directory to write into; created when absent.
#' @param pca_mat optional numeric matrix, rows are cells in the same order,
#'   columns are embedding dimensions. When supplied, CoSPAR builds its
#'   similarity graph on it and never touches `count_mat` beyond DE testing.
#'   When `NULL`, `run_cospar.py` computes highly variable genes and PCA itself.
#' @param verbose numeric.
#'
#' @returns the manifest as a list, invisibly.
export_for_cospar <- function(count_mat,
                              time_info,
                              clone_id,
                              state_info,
                              export_dir,
                              pca_mat = NULL,
                              verbose = 0){
  stopifnot(length(rownames(count_mat)) == nrow(count_mat),
            length(colnames(count_mat)) == ncol(count_mat),
            anyDuplicated(rownames(count_mat)) == 0,
            anyDuplicated(colnames(count_mat)) == 0,
            length(time_info) == nrow(count_mat),
            length(clone_id) == nrow(count_mat),
            length(state_info) == nrow(count_mat),
            !anyNA(time_info), !anyNA(state_info))
  if(min(count_mat) < 0){
    stop("`count_mat` has negative entries; CoSPAR expects non-negative counts")
  }
  if(length(unique(time_info)) != 2){
    stop("`time_info` must carry exactly two labels, found: ",
         paste(unique(time_info), collapse = ", "))
  }

  # A clone seen at only one time point contributes no transition constraint,
  # and CoSPAR aborts when no clone spans both. Report the count up front.
  clone_tab <- table(clone_id[!is.na(clone_id)], time_info[!is.na(clone_id)])
  num_multitime <- sum(apply(clone_tab > 0, 1, all))
  if(num_multitime == 0){
    stop("no clone is observed at both time points; CoSPAR cannot run")
  }
  if(verbose > 0){
    print(paste0(num_multitime, " of ", nrow(clone_tab),
                 " clones are observed at both time points"))
  }

  dir.create(export_dir, showWarnings = FALSE, recursive = TRUE)

  count_sparse <- methods::as(methods::as(count_mat, "CsparseMatrix"),
                              "generalMatrix")
  Matrix::writeMM(count_sparse, file = file.path(export_dir, "counts.mtx"))
  writeLines(rownames(count_mat), file.path(export_dir, "cells.txt"))
  writeLines(colnames(count_mat), file.path(export_dir, "genes.txt"))

  obs_df <- data.frame(cell_id = rownames(count_mat),
                       time_info = as.character(time_info),
                       clone_id = as.character(clone_id),
                       state_info = as.character(state_info),
                       stringsAsFactors = FALSE)
  utils::write.csv(obs_df, file.path(export_dir, "obs.csv"),
                   row.names = FALSE, quote = TRUE)

  bool_pca <- !is.null(pca_mat)
  if(bool_pca){
    stopifnot(nrow(pca_mat) == nrow(count_mat),
              all(rownames(pca_mat) == rownames(count_mat)))
    utils::write.csv(pca_mat, file.path(export_dir, "X_pca.csv"),
                     row.names = TRUE, quote = FALSE)
  }

  manifest <- list(bool_pca = bool_pca,
                   num_cells = nrow(count_mat),
                   num_clones_multitime = num_multitime,
                   num_genes = ncol(count_mat),
                   r_version = R.version.string,
                   time_labels = sort(unique(as.character(time_info))),
                   written = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  jsonlite::write_json(manifest, file.path(export_dir, "manifest.json"),
                       auto_unbox = TRUE, pretty = TRUE)

  invisible(manifest)
}

#' Read CoSPAR results back
#'
#' @param export_dir the directory passed to `export_for_cospar()`.
#'
#' @returns a list with elements `dge_a_df` and `dge_b_df` (data frames from
#'   `cs.tl.differential_genes`, genes up in the fate-A and fate-B progenitors
#'   respectively, with `gene`, `ratio`, `Qvalue` columns), `fate_bias` (named
#'   numeric over the earlier-time-point cells that CoSPAR scored, in `[0, 1]`,
#'   above 0.5 meaning biased toward fate A), `obs_df` (the full CoSPAR obs
#'   table), `progenitor_a` and `progenitor_b` (named logical over the same
#'   cells), and `run` (the parameters and versions the Python side recorded).
import_cospar_results <- function(export_dir){
  out_dir <- file.path(export_dir, "cospar_out")
  stopifnot(dir.exists(out_dir))

  run <- jsonlite::read_json(file.path(out_dir, "run.json"))
  # `read.csv` rewrites the literal `*` in cospar's `fate_bias_<src>_A*B`
  # column name to `.`; `check.names = FALSE` keeps the Python spelling.
  obs_df <- utils::read.csv(file.path(out_dir, "obs.csv"),
                            row.names = 1, check.names = FALSE,
                            stringsAsFactors = FALSE)

  t1_idx <- which(obs_df$time_info == run$t1)
  fate_bias <- obs_df[t1_idx, run$fate_bias_column]
  names(fate_bias) <- rownames(obs_df)[t1_idx]
  progenitor_a <- as.logical(obs_df[t1_idx, run$progenitor_a_column])
  progenitor_b <- as.logical(obs_df[t1_idx, run$progenitor_b_column])
  names(progenitor_a) <- names(fate_bias)
  names(progenitor_b) <- names(fate_bias)

  dge_a_df <- utils::read.csv(file.path(out_dir, "dge_A.csv"),
                              stringsAsFactors = FALSE)
  dge_b_df <- utils::read.csv(file.path(out_dir, "dge_B.csv"),
                              stringsAsFactors = FALSE)

  list(dge_a_df = dge_a_df,
       dge_b_df = dge_b_df,
       fate_bias = fate_bias,
       obs_df = obs_df,
       progenitor_a = progenitor_a,
       progenitor_b = progenitor_b,
       run = run)
}
