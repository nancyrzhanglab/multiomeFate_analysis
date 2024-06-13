library(Seurat)
library(ggplot2)

in_dir <- '~/Dropbox/Thesis/Lineage_trace/data/Watermelon/'
out_dir <- '~/Dropbox/Thesis/Lineage_trace/output/Watermelon/'

# ==============================================================================
# Read data
# ==============================================================================
load(paste0(in_dir, 'PC9_time_course_fasttopics.RData'))
metadat <- seurat_object@meta.data

metadat_d14 <- metadat[metadat$time_point == 14, ]
lineage_fate_d14 <- metadat_d14 %>% 
  group_by(lineage_barcode, sample_name) %>% 
  summarise(n_cells_d14 = n())

metadat_d7 <- metadat[metadat$time_point == 7, ]
lineage_fate_d7 <- metadat_d7 %>% 
  group_by(lineage_barcode) %>% 
  summarise(n_cells_d7 = n())

d7_d14 <- merge(lineage_fate_d7, lineage_fate_d14, by = 'lineage_barcode', all.y = T)
d7_d14[is.na(d7_d14)] <- 0
d7_d14 <- d7_d14[d7_d14$lineage_barcode != 0, ]

d7_d14_w <- reshape(d7_d14, idvar = c('lineage_barcode', 'n_cells_d7'), timevar = "sample_name", direction = "wide")
rownames(d7_d14_w) <- d7_d14_w$lineage_barcode
d7_d14_w <- subset(d7_d14_w, select = -c(lineage_barcode))

ggpairs(d7_d14_w)
ggsave('~/Downloads/fate_specific_size_comp.png', dpi = 300, width = 15.5, height = 15.5)

