gene = "NRG1"

gene <- gene_vec[i]
print(gene)
print("Computing cutmat")

cutmat_winning <- extract_cutmatrix(
  object = all_data,
  gene = gene,
  cells = winning_cells
)
cutmat_dying <- extract_cutmatrix(
  object = all_data,
  gene = gene,
  cells = dying_cells
)
cutmat_all <- rbind(cutmat_winning, cutmat_dying)

peak_mat <- multiomeFate:::extract_peaks(
  object = all_data,
  gene = gene
)
bin_midpoints <- multiomeFate:::compute_bin_midpoints(peak_mat)
bin_limits <- c(bin_midpoints[1] - abs(bin_midpoints[2]-bin_midpoints[1]),
                bin_midpoints[7] + abs(bin_midpoints[2]-bin_midpoints[1]))
peak_locations <- multiomeFate:::compute_peak_locations(peak_mat)
peak_prior <- compute_peak_prior(mat = cutmat_all,
                                 peak_mat = peak_mat)

res_winning <- multiomeFate:::peak_mixture_modeling(
  bin_limits = bin_limits,
  bin_midpoints = bin_midpoints, 
  cutmat = cutmat_winning, 
  peak_locations = peak_locations,
  peak_prior = peak_prior,
  bool_freeze_prior = T,
  verbose = 3
)
round(res_winning$theta_vec,2)

###############################


