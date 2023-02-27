rm(list=ls())
library(Seurat)
library(Signac)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")

## from https://github.com/stuart-lab/signac/blob/master/R/visualization.R#L888

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

Seurat::Idents(all_data) <- "dataset"
Seurat::DefaultAssay(all_data) <- "ATAC"

gene <- "ANXA1"
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

##################3

cells <- Signac:::SetIfNull(x = cells, y = colnames(x = object))
assay <- Signac:::SetIfNull(x = assay, y = Seurat::DefaultAssay(object = object))
if (!inherits(x = assay, what = "list")) {
  assay <- list(assay)
}
if (!is.null(x = group.by)) {
  Seurat::Idents(object = object) <- group.by
}
if (!is.null(x = idents)) {
  ident.cells <- Signac:::WhichCells(object = object, idents = idents)
  cells <- intersect(x = cells, y = ident.cells)
}

#####

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

####

i <- 1
reads.per.group <- Signac:::AverageCounts(
  object = object,
  assay = assay[[i]],
  group.by = group.by,
  verbose = FALSE
)

################################################

# zoom in on the CutMatrix function
assay = assay[[i]]
verbose = FALSE
assay <- Signac:::SetIfNull(x = assay, y = DefaultAssay(object = object))
cells <- Signac:::SetIfNull(x = cells, y = colnames(x = object))
fragments <- Signac:::Fragments(object = object[[assay]])

res <- list()
for (i in seq_along(along.with = fragments)) {
  fragment.path <- Signac:::GetFragmentData(object = fragments[[i]], slot = "path")
  cellmap <- Signac:::GetFragmentData(object = fragments[[i]], slot = "cells")
  tabix.file <- Rsamtools:::TabixFile(file = fragment.path)
  open(con = tabix.file)
  # remove regions that aren't in the fragment file
  seqnames.in.both <- intersect(
    x = GenomeInfoDb:::seqnames(x = region),
    y = Rsamtools:::seqnamesTabix(file = tabix.file)
  )
  region <- GenomeInfoDb:::keepSeqlevels(
    x = region,
    value = seqnames.in.both,
    pruning.mode = "coarse"
  )
  if (length(x = region) != 0) {
    cm <- Signac:::SingleFileCutMatrix(
      region = region,
      cellmap = cellmap,
      tabix.file = tabix.file,
      cells = cells,
      verbose = FALSE
    )
    res[[i]] <- cm
  }
  close(con = tabix.file)
}

####
# let's look into Signac:::SingleFileCutMatrix
fragments <- Signac:::GetReadsInRegion(
  region = region,
  cellmap = cellmap,
  cells = cells,
  tabix.file = tabix.file,
  verbose = verbose
)
start.lookup <- BiocGenerics:::start(x = region)
names(start.lookup) <- seq_along(region)
fragstarts <- start.lookup[fragments$ident] + 1
cut.df <- data.frame(
  position = c(fragments$start, fragments$end) - fragstarts,
  cell = c(fragments$cell, fragments$cell),
  stringsAsFactors = FALSE
)

################################################

cutmat <- Signac:::CutMatrix(
  object = object,
  region = region,
  assay = assay[[i]],
  cells = cells,
  verbose = FALSE
)
colnames(cutmat) <- BiocGenerics::start(x = region):BiocGenerics::end(x = region)
group.scale.factors <- suppressWarnings(reads.per.group * cells.per.group)
scale.factor <- Signac:::SetIfNull(
  x = scale.factor, y = median(x = group.scale.factors)
)
cm.list[[i]] <- cutmat
sf.list[[i]] <- scale.factor
gsf.list[[i]] <- group.scale.factors
names(x = cm.list) <- unlist(x = assay)
p <- Signac:::CoverageTrack(
  cutmat = cm.list,
  region = region,
  group.scale.factors = gsf.list,
  scale.factor = sf.list,
  window = window,
  ymax = ymax,
  split.assays = split.assays,
  assay.scale = assay.scale,
  obj.groups = obj.groups,
  region.highlight = region.highlight,
  downsample.rate = downsample.rate,
  max.downsample = max.downsample
)

ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6b/Writeup6b_coverageplot_track-test.png"),
                p, device = "png", width = 10, height = 10, units = "in")

#######################

# digging into the coverage track
# https://github.com/stuart-lab/signac/blob/master/R/visualization.R

cutmat = cm.list
group.scale.factors = gsf.list
scale.factor = sf.list

window.size <- IRanges::width(x = region)
levels.use <- levels(x = obj.groups)
chromosome <- as.character(x = seqnames(x = region))
start.pos <- IRanges::start(x = region)
end.pos <- IRanges::end(x = region)
multicov <- length(x = cutmat) > 1

cov.df <- data.frame()
i <- 1
## this is one the heavy-lifting functions
coverages <- Signac:::ApplyMatrixByGroup(
  mat = cutmat[[i]],
  fun = colSums,
  groups = obj.groups,
  group.scale.factors = group.scale.factors[[i]],
  scale.factor = scale.factor[[i]],
  normalize = TRUE
)
coverages <- dplyr::group_by(.data = coverages, group)
coverages <- dplyr::mutate(.data = coverages, 
                           coverage = RcppRoll::roll_sum(
                             x = norm.value, 
                             n = window, 
                             fill = NA, 
                             align = "center"
                           ))
coverages <- dplyr::ungroup(x = coverages)
dim(coverages)
dim(coverages)
coverages[45:55,]

coverages <- coverages[!is.na(x = coverages$coverage), ]
coverages <- dplyr::group_by(.data = coverages, group)
dim(coverages)
# sampling <- min(max.downsample, window.size * downsample.rate)
# set.seed(seed = 1234)
# coverages <- dplyr::slice_sample(.data = coverages, n = as.integer(x = sampling))
# coverages$Assay <- names(x = cutmat)[[i]]
# dim(coverages)
# cov.df <- rbind(cov.df, coverages)
# 
# coverages <- cov.df
# coverages$Assay <- factor(x = coverages$Assay, levels = names(x = cutmat))
# coverages$assay_group <- paste(coverages$group, coverages$Assay, sep = "_")
# 
# if (!is.null(x = levels.use)) {
#   colors_all <- scales::hue_pal()(length(x = levels.use))
#   names(x = colors_all) <- levels.use
#   coverages$group <- factor(x = coverages$group, levels = levels.use)
# }
# covmax <- signif(x = max(coverages$coverage, na.rm = TRUE), digits = 2)
# ymax <- covmax
# ymin <- 0

tmp <- as.data.frame(coverages)
idx <- which(tmp[,"group"] == "week5_DABTRAM")
x_vec <- tmp[idx,"position"]
y_vec <- tmp[idx,"coverage"]

png(file = "../../../../out/figures/Writeup6b/Writeup6b_coverageplot_track-test2.png",
    height = 1500, width = 4500, units = "px", res = 300)
par(mar = c(4,4,0.5,0.5))
plot(x_vec, y_vec, pch = 16,
     xlab = "Position", ylab = "Normalized coverage")
graphics.off()
