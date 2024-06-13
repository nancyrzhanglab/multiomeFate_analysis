library(tidyr)
library(ggplot2)

# ==============================================================================
# Read data
# ==============================================================================
in_dir <- '~/Dropbox/Thesis/Lineage_trace/data/Watermelon/'
load(paste0(in_dir,'PC9_time_course_fasttopics.RData'))
watermelon_all <- seurat_object@meta.data
watermelon_all$Cell.ID <- rownames(watermelon_all)
watermelon_all <- watermelon_all %>% drop_na()
w_day0 <- watermelon_all[watermelon_all$time_point == 0, ]
w_day3 <- watermelon_all[watermelon_all$time_point == 3, ]
w_day7 <- watermelon_all[watermelon_all$time_point == 7, ]
w_day14 <- watermelon_all[watermelon_all$time_point == 14, ]

# ==============================================================================
# Helper
# ==============================================================================
compute_dominance_index<-function(lineage_sizes, time = '', doplot=TRUE,
                                  xlab="Cumulative share of lineages from smallest to largest",
                                  ylab="Cumulative share of cells"){
  # compute Lorenz curve
  lineage_sizes=sort(lineage_sizes, decreasing=FALSE)
  sizeprop = cumsum(lineage_sizes/sum(lineage_sizes))
  increment = 1/length(lineage_sizes)
  B=sum(increment*sizeprop)-increment/2
  gini = 1-2*B
  png(paste0('~/Downloads/', time, '.png'),units="in",  width = 5, height = 5, res=200)
  if(doplot){
    plot(seq(increment, 1, increment), sizeprop,
         xlab=xlab, ylab=ylab, xlim=c(0,1), ylim=c(0,1), main=paste("Gini = ", format(gini, digits=3), '( ', time, ' )'),
         xaxt="n", yaxt="n")
    abline(0,1, col="blue")
  }
  dev.off()
}

w_day0_lineage_size <- w_day0 %>% 
  group_by(lineage_barcode) %>% 
  summarise(lineage_size = n_distinct(Cell.ID))
w_day3_lineage_size <- w_day3 %>% 
  group_by(lineage_barcode) %>% 
  summarise(lineage_size = n_distinct(Cell.ID))
w_day7_lineage_size <- w_day7 %>% 
  group_by(lineage_barcode) %>% 
  summarise(lineage_size = n_distinct(Cell.ID))
w_day14_lineage_size <- w_day14 %>% 
  group_by(lineage_barcode) %>% 
  summarise(lineage_size = n_distinct(Cell.ID))


compute_dominance_index(w_day0_lineage_size$lineage_size, time = 'day0')
compute_dominance_index(w_day3_lineage_size$lineage_size, time = 'day3')
compute_dominance_index(w_day7_lineage_size$lineage_size, time = 'day7')
compute_dominance_index(w_day14_lineage_size$lineage_size, time = 'day14')
