library(Seurat)
setwd("../multiome_fate/")
source("./scripts/precursor_functions.R")

brain=readRDS("/Volumes/Cloud/Projects/multiome_fate/brain/data/obj_seurat_links.rds")
linkpeak_result=brain@assays$ATAC@links
linkpeak_result=linkpeak_result[which(linkpeak_result$zscore>0)]
brain=readRDS("/Volumes/Cloud/Projects/multiome_fate/brain/var_selection/mbrain.rds")
table(Idents(brain))

# prep matrices for RNA and ATAC
matrices=PrepMatrices(obj = brain, linkpeak_result = linkpeak_result)
gene_count=matrices$gene_count
gene_score_atac=matrices$gene_score_atac



# find var genes
S1 =c("Oligodendrocyte")
S0=NULL
#c( "Cortical or hippocampal glutamatergic","Mixed region GABAergic","Forebrain GABAergic","Forebrain glutamatergic","Midbrain glutamatergic","Hindbrain glycinergic","Oligodendrocyte")
pval.thresh.DE=1e-100
S1.genes=FindVarGenes(obj=brain,S1=S1, S0=S0, pval.thresh.DE = 1e-100)


## compute precursor scores
scores=ComputeScores(obj=brain, gene_count = gene_count, gene_score_atac = gene_score_atac, S1=S1, S1.genes = S1.genes)
rna_score=scores$rna_score
precs_score=scores$atac_score


## visualization
library(cowplot)
library(ggplot2)

ind=which(rownames(gene_count) %in% S1.genes)
df_plot=data.frame(UMAP1=brain@reductions$umap@cell.embeddings[,1],UMAP2=brain@reductions$umap@cell.embeddings[,2] 
                   ,gene_count=t(gene_count[ind,]), gene_score_atac=t(gene_score_atac[ind,]))

# plot each gene rna and atac
pdf("/Volumes/Cloud/Projects/multiome_fate/brain/precursors/precursor_Oligodendrocyte_glio_only.pdf")

for(ii in 1:length(S1.genes)){
  gene=S1.genes[ii]
  
  tryCatch({
    
    p1=ggplot(data=df_plot,aes(x=UMAP1, y=UMAP2))+
      geom_point(aes(color= get(paste0("gene_count.",gene))), alpha=0.4)+
      scale_color_gradient(name = paste0(gene," RNA"), low = "grey", high = "blue")
    #scale_color_grey() + theme_classic()
    
    p2=ggplot(data=df_plot,aes(x=UMAP1, y=UMAP2))+
      geom_point(aes(color= get(paste0("gene_score_atac.",gene))), alpha=0.4)+
      scale_color_gradient(name = paste0(gene," ATAC"), low = "grey", high = "blue")
    #scale_color_grey() + theme_classic()
    
    brain[[paste0("gene_score_atac.",gene)]]=df_plot[paste0("gene_score_atac.",gene)]
    brain[[paste0("gene_count.",gene)]]=df_plot[paste0("gene_count.",gene)]
    
    
    p3=DotPlot(brain, features=c(paste0("gene_count.",gene),paste0("gene_score_atac.",gene)))
    
    
    pp_combine=plot_grid(p1|p2,p3,ncol = 1, nrow = 2)
    
    print(pp_combine)
    print(ii)
    
  }, error=function(e){})
  
  
}

dev.off()



## plot precursor scores

df_plot$rna_score=rna_score
p1=ggplot(data=df_plot,aes(x=UMAP1, y=UMAP2))+
  geom_point(aes(color= rna_score), alpha=0.4)+
  scale_color_gradient(name = paste0("Score"), low = "grey", high = "blue")+
  ggtitle(paste0(S1 ," RNA across ",length(ind)," genes"))


df_plot$precs_score=precs_score

p2=ggplot(data=df_plot,aes(x=UMAP1, y=UMAP2))+
  geom_point(aes(color= precs_score), alpha=0.4)+
  scale_color_gradient(name = paste0("Score"), low = "grey", high = "blue")+
  ggtitle(paste0(S1 ," ATAC across ",length(ind)," genes"))



p1|p2
dim_plot=DimPlot(brain, reduction = "umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
print(dim_plot)


