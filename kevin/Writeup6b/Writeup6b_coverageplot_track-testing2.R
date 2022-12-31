rm(list=ls())
library(Seurat)
library(Signac)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")

##############

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

Seurat::Idents(all_data) <- "dataset"
Seurat::DefaultAssay(all_data) <- "ATAC"

treatment_vec <- c("CIS", "COCL2", "DABTRAM")
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
gene_vec <- sort(unique(unlist(keygenes)))
gene_vec <- gene_vec[which(gene_vec %in% rownames(all_data[["RNA"]]))]

treatment <- treatment_vec[1]
ident_vec <- all_data$dataset
names(ident_vec) <- colnames(all_data)
lineage_names <- rownames(tab_mat)[which(tab_mat[,paste0("day10_", treatment)] >= 20)]
cell_names1 <- colnames(all_data)[which(all_data$assigned_lineage %in% lineage_names)]
cell_names2 <- colnames(all_data)[which(all_data$dataset == "day0")]
cell_names_winning <- intersect(cell_names1, cell_names2)
cell_names_losing <- setdiff(cell_names2, cell_names1)
ident_vec[cell_names_winning] <- paste0("day0_win_", treatment)
ident_vec[cell_names_losing] <- paste0("day0_lose_", treatment)
all_data$ident <- ident_vec
Seurat::Idents(all_data) <- "ident"

###

gene <- "COL6A2"
object <- all_data
region <- gene
features <- gene
extend.upstream <- 1000
extend.downstream <- 1000

assay = NULL
split.assays = FALSE
assay.scale = "common"
show.bulk = FALSE
expression.assay = NULL
expression.slot = "data"
annotation = TRUE
peaks = TRUE
peaks.group.by = NULL
ranges = NULL
ranges.group.by = NULL
ranges.title = "Ranges"
region.highlight = NULL
links = TRUE
tile = FALSE
tile.size = 100
tile.cells = 100
bigwig = NULL
bigwig.type = "coverage"
bigwig.scale = "common"
group.by = NULL
window = 100
ymax = NULL
scale.factor = NULL
cells = NULL
idents = NULL
sep = c("-", "-")
heights = NULL
max.downsample = 3000
downsample.rate = 0.1

##############

assay <- Signac:::SetIfNull(x = assay, y = Seurat::DefaultAssay(object = object))

region <- Signac:::FindRegion(
  object = object,
  region = region,
  sep = sep,
  assay = assay[[1]],
  extend.upstream = extend.upstream,
  extend.downstream = extend.downstream
)
cells.per.group <- Signac:::CellsPerGroup(
  object = object,
  group.by = group.by
)
obj.groups <- Signac:::GetGroups(
  object = object,
  group.by = group.by,
  idents = idents
)
cm.list <- list()
sf.list <- list()
gsf.list <- list()

##

i <- 1
reads.per.group <- Signac:::AverageCounts(
  object = object,
  assay = assay[[i]],
  group.by = group.by,
  verbose = FALSE
)
cutmat <- Signac:::CutMatrix(
  object = object,
  region = region,
  assay = assay[[i]],
  cells = cells,
  verbose = FALSE
)
colnames(cutmat) <- (BiocGenerics::start(x = region)):(BiocGenerics::end(x = region))
group.scale.factors <- suppressWarnings(reads.per.group * cells.per.group)
scale.factor <- Signac:::SetIfNull(
  x = scale.factor, y = median(x = group.scale.factors)
)
cm.list[[i]] <- cutmat
sf.list[[i]] <- scale.factor
gsf.list[[i]] <- group.scale.factors
names(x = cm.list) <- unlist(x = assay)

# p <- Signac:::CoverageTrack(
#   cutmat = cm.list,
#   region = region,
#   group.scale.factors = gsf.list,
#   scale.factor = sf.list,
#   window = window,
#   ymax = ymax,
#   split.assays = split.assays,
#   assay.scale = assay.scale,
#   obj.groups = obj.groups,
#   region.highlight = region.highlight,
#   downsample.rate = downsample.rate,
#   max.downsample = max.downsample
# )

###########

cutmat = cm.list
group.scale.factors = gsf.list
scale.factor = sf.list

window.size <- IRanges::width(x = region)
levels.use <- levels(x = obj.groups)
chromosome <- as.character(x = seqnames(x = region))
start.pos <- BiocGenerics::start(x = region)
end.pos <- BiocGenerics::end(x = region)
multicov <- length(x = cutmat) > 1

cov.df <- data.frame()

i <- 1
coverages <- Signac:::ApplyMatrixByGroup(
  mat = cutmat[[i]],
  fun = colSums,
  groups = obj.groups,
  group.scale.factors = group.scale.factors[[i]],
  scale.factor = scale.factor[[i]],
  normalize = TRUE
)
if (!is.na(x = window)) {
  coverages <- dplyr::group_by(.data = coverages, group)
  coverages <- dplyr::mutate(.data = coverages, coverage = RcppRoll::roll_sum(
    x = norm.value, n = window, fill = NA, align = "center"
  )) ## KEVIN: THIS SEEMS TO BE THE MAIN FUNCTION... JUST ADDING?
  coverages <- dplyr::ungroup(x = coverages)
} else {
  coverages$coverage <- coverages$norm.value
}

coverages <- coverages[!is.na(x = coverages$coverage), ]
coverages <- dplyr::group_by(.data = coverages, group)
sampling <- min(max.downsample, window.size * downsample.rate)
set.seed(seed = 1234)
coverages <- dplyr::slice_sample(.data = coverages, n = as.integer(x = sampling))
coverages$Assay <- names(x = cutmat)[[i]]

cov.df <- rbind(cov.df, coverages)
coverages <- cov.df
coverages$Assay <- factor(x = coverages$Assay, levels = names(x = cutmat))
coverages$assay_group <- paste(coverages$group, coverages$Assay, sep = "_")

########

# restore factor levels
if (!is.null(x = levels.use)) {
  colors_all <- scales::hue_pal()(length(x = levels.use))
  names(x = colors_all) <- levels.use
  coverages$group <- factor(x = coverages$group, levels = levels.use)
}
covmax <- signif(x = max(coverages$coverage, na.rm = TRUE), digits = 2)
if (is.null(x = ymax)) {
  ymax <- covmax
} else if (is.character(x = ymax)) {
  if (!startsWith(x = ymax, prefix = "q")) {
    stop("Unknown ymax requested. Must be NULL, a numeric value, or 
           a quantile denoted by 'qXX' with XX the desired quantile value,
           e.g. q95 for 95th percentile")
  }
  percentile.use <- as.numeric(
    x = sub(pattern = "q", replacement = "", x = as.character(x = ymax))
  ) / 100
  ymax <- covmax * percentile.use
}
ymin <- 0

# perform clipping
coverages$coverage[coverages$coverage > ymax] <- ymax 
