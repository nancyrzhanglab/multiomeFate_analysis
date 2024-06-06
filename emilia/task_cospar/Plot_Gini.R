library(ggplot2)

# ==============================================================================
# Read data
# ==============================================================================
in_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/LARRY_hematopoiesis/'
load(paste0(in_dir,'/Writeup8_larry-dataset_step2_lineage-plotting.RData'))
cospar_all <- seurat_object@meta.data
cospar_all$Cell.ID <- rownames(cospar_all)
cospar_day2 <- cospar_all[cospar_all$Time.point == '2', ]
cospar_day4 <- cospar_all[cospar_all$Time.point == '4', ]
cospar_day6 <- cospar_all[cospar_all$Time.point == '6', ]

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

cospar_day2_lineage_size <- cospar_day2 %>% 
  group_by(assigned_lineage) %>% 
  summarise(lineage_size = n_distinct(Cell.ID))
cospar_day4_lineage_size <- cospar_day4 %>% 
  group_by(assigned_lineage) %>% 
  summarise(lineage_size = n_distinct(Cell.ID))
cospar_day6_lineage_size <- cospar_day6 %>% 
  group_by(assigned_lineage) %>% 
  summarise(lineage_size = n_distinct(Cell.ID))

compute_dominance_index(cospar_day2_lineage_size$lineage_size, time = 'day2')
compute_dominance_index(cospar_day4_lineage_size$lineage_size, time = 'day4')
compute_dominance_index(cospar_day6_lineage_size$lineage_size, time = 'day6')
