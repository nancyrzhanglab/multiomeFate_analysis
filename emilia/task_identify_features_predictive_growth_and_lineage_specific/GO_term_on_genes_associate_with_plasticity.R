library(tidyverse)
library(msigdbr)
library(clusterProfiler)
library(fgsea)
library(org.Hs.eg.db)
library(AnnotationDbi)

# ==============================================================================
# Read data
# ==============================================================================
input_dir <- '~/Dropbox/Thesis/Lineage_trace/outputs/task0_explore_lineage_variability/features_associated_with_plasticity/'
genes_to_test <- read.csv(paste0(input_dir, 'day10_DABTRAM_gene_exp_mean_day10_plasticity_correlation_TOP1000.csv'))

input_dir <- '~/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/'
# genes_to_test <- read.csv(paste0(input_dir, 'common_genes_in_corr_with_growth_v2.csv'))
genes_to_test <- read.csv(paste0(input_dir, 'common_genes_in_corr_with_growth_v3_cis_cocl2_clust3.csv'))


#Hallmark
H <- as.data.frame(msigdbr(species = "Homo sapiens",
                           category = "H"))


genes_to_test.entrez <-mapIds(org.Hs.eg.db, keys = genes_to_test$gene,
       column = "ENTREZID", keytype = "SYMBOL")
H.entrez <- H[, c('gs_name', 'entrez_gene')]

enrich.H <- enricher(gene = genes_to_test.entrez, TERM2GENE = H.entrez)

enrich.H.df <- enrich.H@result %>%
  #separate ratios into 2 columns of data
  separate(BgRatio, into=c("size.term","size.category"), sep="/") %>%
  separate(GeneRatio, into=c("size.overlap.term", "size.overlap.category"),
           sep="/") %>%
  #convert to numeric
  mutate_at(vars("size.term","size.category",
                 "size.overlap.term","size.overlap.category"),
            as.numeric) %>%
  #Calculate k/K
  mutate("k.K"=size.overlap.term/size.term)

enrich.H.df %>%
  filter(p.adjust <= 0.05) %>%
  #Beautify descriptions by removing _ and HALLMARK
  mutate(Description = gsub("HALLMARK_","", Description),
         Description = gsub("_"," ", Description)) %>%
  ggplot(aes(x=reorder(Description, k.K), #Reorder gene sets by k/K values
             y=k.K)) +
  geom_col() +
  theme_classic() +
#Some more customization to pretty it up
#Flip x and y so long labels can be read
  coord_flip() +
  #fix labels
  labs(y="Significant genes in set / Total genes in set",
       x="Gene set")
