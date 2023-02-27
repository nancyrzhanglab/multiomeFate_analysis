rm(list=ls())
library(Seurat)
library(Signac)

load("../../../../out/kevin/Writeup6b/Writeup6b_genomebinmatrix.RData")
load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")

genome_bin_matrix[5000:5005,1:5]

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

# region2 <- Signac:::FindRegion(
#   object = object,
#   region = region,
#   sep = sep,
#   assay = assay[[1]],
#   extend.upstream = extend.upstream,
#   extend.downstream = extend.downstream
# )
region <- GenomicRanges::GRanges(Rle("chr1"),
                                 IRanges::IRanges(29450001, 29500000))

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

cutmat <- Signac:::CutMatrix(
  object = object,
  region = region,
  assay = assay[[i]],
  cells = cells,
  verbose = FALSE
)
colnames(cutmat) <- start(x = region):end(x = region)
group.scale.factors <- suppressWarnings(reads.per.group * cells.per.group)
scale.factor <- Signac:::SetIfNull(
  x = scale.factor, y = median(x = group.scale.factors)
)
cm.list[[i]] <- cutmat
sf.list[[i]] <- scale.factor
gsf.list[[i]] <- group.scale.factors

########################

# now let's see if the two match up
vec1 <- genome_bin_matrix[5000,]
vec1[1:10]

vec2 <- rowSums(cutmat)
vec2[1:10]

stats::cor(vec1, vec2) # almost very very close

##################################
##################################
##################################
