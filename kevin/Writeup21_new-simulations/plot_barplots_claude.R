# plot_barplots_claude.R
# Kevin Z. Lin (drafted by Claude), 2026-09-22
#
# Figure 4 barplots from the CSVs of `make_csvs_claude.R`
# (simulation_design_claude.md, Section 8): per axis, one bar per method at
# each level, height the mean Spearman metric across replicates with the
# replicates overplotted as points, x-axis labelled with the target
# statistic and the realized mean beneath it (for the AED axis also the
# per-clone 5th-95th percentile range). A second figure per axis plots the
# metric against the realized statistic as a check on the calibration, and
# a third shows the score-level diagnostics (correlation of each method's
# fate score with the true fate potential).
#
# Run:  Rscript plot_barplots_claude.R
# Dry:  WRITEUP21_DRY=1 Rscript plot_barplots_claude.R
# PCs:  WRITEUP21_PCA=20 Rscript plot_barplots_claude.R  (reads the `_pca20` CSVs)

library(ggplot2)

rm(list = ls())

# Paths ------------------------------------------------------------------------

repo_dir <- file.path("/Users/kevinlin/Library/CloudStorage/Dropbox",
                      "Collaboration-and-People/archive/Nancy/multiomeFate",
                      "git/multiomeFate_analysis")
csv_dir <- file.path(repo_dir, "csv", "kevin", "Writeup21")
fig_dir <- file.path(repo_dir, "fig", "kevin", "Writeup21")

bool_dry <- Sys.getenv("WRITEUP21_DRY") == "1"
d_pca <- as.numeric(Sys.getenv("WRITEUP21_PCA", unset = "10"))
suffix <- paste0(if(d_pca != 10) paste0("_pca", d_pca) else "",
                 if(bool_dry) "_dry" else "")

# Settings ---------------------------------------------------------------------

method_level_vec <- c("CYFER", "CoSPAR", "Lineage-DE")
method_color_vec <- c(CYFER = "#D55E00", CoSPAR = "#0072B2",
                      "Lineage-DE" = "#009E73")
axis_label_list <- list(
  gini = "Gini of t2 clone sizes (target / realized)",
  aed = "Mean squared AED (target / realized; per-clone 5th-95th pct)")

# Helpers ----------------------------------------------------------------------

# x-axis label per level: the target with the realized mean beneath, and for
# the AED axis the per-clone 5th-95th percentile range on a third line.
.level_label_vec <- function(summary_df, axis){
  level_df <- unique(summary_df[, c("level_idx", "level", "gini_t2_mean",
                                    "aed_mean_mean", "aed_q05_mean",
                                    "aed_q95_mean")])
  level_df <- level_df[order(level_df$level_idx), ]
  if(axis == "gini"){
    label_vec <- paste0(level_df$level, "\n(",
                        sprintf("%.2f", level_df$gini_t2_mean), ")")
  } else {
    label_vec <- paste0(level_df$level, "\n(",
                        sprintf("%.2f", level_df$aed_mean_mean), ")\n[",
                        sprintf("%.2f", level_df$aed_q05_mean), ", ",
                        sprintf("%.2f", level_df$aed_q95_mean), "]")
  }
  names(label_vec) <- level_df$level_idx
  label_vec
}

.plot_barplot <- function(summary_df, detail_df, axis, metric, ylab){
  summary_df$method <- factor(summary_df$method, levels = method_level_vec)
  detail_df$method <- factor(detail_df$method, levels = method_level_vec)
  summary_df$level_factor <- factor(summary_df$level_idx)
  detail_df$level_factor <- factor(detail_df$level_idx)
  mean_col <- paste0(metric, "_mean")
  label_vec <- .level_label_vec(summary_df, axis)

  plot1 <- ggplot2::ggplot(summary_df,
                           ggplot2::aes(x = level_factor, y = .data[[mean_col]],
                                        fill = method))
  plot1 <- plot1 + ggplot2::geom_col(position = ggplot2::position_dodge(0.8),
                                     width = 0.7)
  plot1 <- plot1 + ggplot2::geom_point(
    data = detail_df,
    ggplot2::aes(x = level_factor, y = .data[[metric]], group = method),
    position = ggplot2::position_dodge(0.8), size = 1.2, color = "black",
    inherit.aes = FALSE)
  plot1 <- plot1 + ggplot2::geom_hline(yintercept = 0, linewidth = 0.3)
  plot1 <- plot1 + ggplot2::scale_fill_manual(values = method_color_vec)
  plot1 <- plot1 + ggplot2::scale_x_discrete(labels = label_vec)
  plot1 <- plot1 + ggplot2::labs(x = axis_label_list[[axis]], y = ylab,
                                 fill = NULL)
  plot1 <- plot1 + ggplot2::theme_bw()
  plot1 <- plot1 + ggplot2::theme(legend.position = "top")
  plot1
}

.plot_vs_realized <- function(detail_df, axis){
  detail_df$method <- factor(detail_df$method, levels = method_level_vec)
  realized_col <- if(axis == "gini") "gini_t2" else "aed_mean"
  plot1 <- ggplot2::ggplot(detail_df,
                           ggplot2::aes(x = .data[[realized_col]], y = spearman,
                                        color = method))
  plot1 <- plot1 + ggplot2::geom_line(ggplot2::aes(group = method),
                                      alpha = 0.4)
  plot1 <- plot1 + ggplot2::geom_point(size = 2)
  plot1 <- plot1 + ggplot2::scale_color_manual(values = method_color_vec)
  plot1 <- plot1 + ggplot2::labs(
    x = if(axis == "gini") "Realized Gini of t2 clone sizes" else
      "Realized mean squared AED",
    y = "Spearman(method, truth) over genes", color = NULL)
  plot1 <- plot1 + ggplot2::theme_bw()
  plot1 <- plot1 + ggplot2::theme(legend.position = "top")
  plot1
}

# Figures ----------------------------------------------------------------------

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

for(axis in c("gini", "aed")){
  summary_file <- file.path(csv_dir, paste0("sim_", axis, "_summary", suffix,
                                            ".csv"))
  detail_file <- file.path(csv_dir, paste0("sim_", axis, "_replicate_details",
                                           suffix, ".csv"))
  if(!file.exists(summary_file) || !file.exists(detail_file)){
    print(paste0("Skipping ", axis, ": CSVs not found"))
    next
  }
  summary_df <- utils::read.csv(summary_file, stringsAsFactors = FALSE)
  detail_df <- utils::read.csv(detail_file, stringsAsFactors = FALSE)
  # Sort by level so lines in the realized-statistic plot are monotone.
  detail_df <- detail_df[order(detail_df$method, detail_df[[
    if(axis == "gini") "gini_t2" else "aed_mean"]]), ]

  plot_spearman <- .plot_barplot(summary_df, detail_df, axis, "spearman",
                                 "Spearman(method, truth) over genes")
  ggplot2::ggsave(file.path(fig_dir, paste0("sim_", axis, "_barplot_spearman",
                                            suffix, ".png")),
                  plot_spearman, width = 8, height = 4.5, dpi = 200)

  plot_jaccard <- .plot_barplot(summary_df, detail_df, axis, "jaccard",
                                "Jaccard at BH q < 0.05")
  ggplot2::ggsave(file.path(fig_dir, paste0("sim_", axis, "_barplot_jaccard",
                                            suffix, ".png")),
                  plot_jaccard, width = 8, height = 4.5, dpi = 200)

  # Lineage-DE has no per-cell fate score, so it has no bar here.
  plot_score <- .plot_barplot(summary_df[summary_df$method != "Lineage-DE", ],
                              detail_df[detail_df$method != "Lineage-DE", ],
                              axis, "score_cor",
                              "cor(fate score, true fate potential), t1 cells")
  ggplot2::ggsave(file.path(fig_dir, paste0("sim_", axis, "_barplot_score",
                                            suffix, ".png")),
                  plot_score, width = 8, height = 4.5, dpi = 200)

  plot_realized <- .plot_vs_realized(detail_df, axis)
  ggplot2::ggsave(file.path(fig_dir, paste0("sim_", axis, "_vs_realized",
                                            suffix, ".png")),
                  plot_realized, width = 6, height = 4.5, dpi = 200)

  print(paste0("Wrote ", axis, " figures to ", fig_dir))
}
