#precursor_functions
################
FindVarGenes=function(obj=NULL, S1=NULL, S0=NULL, pval.thresh.DE=1e-100){
par(mfrow=c(2,1))
DefaultAssay(obj)="SCT"
markers = FindMarkers(obj, ident.1=S1, ident.2=S0)
plot(markers$pct.1, markers$pct.2, xlab=paste("Proportion detected in", S1), ylab=paste("Proportion detected in", S0))

# Do more stringent selection.
markers=markers[which(markers$avg_log2FC>0),]
plot(markers$pct.1, markers$pct.2, xlab=paste("Proportion detected in", S1), ylab=paste("Proportion detected in", S0), main='After filtering')
markergenes=row.names(markers)
S1.genes=markergenes[which(markers$p_val<pval.thresh.DE)]
cat(length(S1.genes),"remaining after DE selection.\n")

return(S1.genes)}


##################
PrepMatrices=function(obj=NULL, linkpeak_result=NULL){
library(Matrix)

gene_uni=linkpeak_result$gene[!duplicated(linkpeak_result$gene)]
rna=as.matrix(obj@assays$RNA@counts)
atac=as.matrix(obj@assays$ATAC@counts)

gene_count=matrix(nrow=length(gene_uni),ncol=ncol(rna))
gene_score_atac=matrix(nrow=length(gene_uni),ncol=ncol(rna))
for(ii in 1:length(gene_uni)){
  gene_count[ii,]=rna[which(rownames(rna)==gene_uni[ii]),]
  gene_score_atac[ii,]=colSums(atac[which(rownames(atac) %in% linkpeak_result$peak[which(linkpeak_result$gene==gene_uni[ii])]),, drop=F])
  print(ii)
}

rownames(gene_count)=gene_uni
rownames(gene_score_atac)=gene_uni
colnames(gene_count)=colnames(rna)
colnames(gene_score_atac)=colnames(rna)

mats=list(gene_count=gene_count, gene_score_atac=gene_score_atac)
return(mats)
}


##############
#Compute Precursor scores
ComputeScores=function(obj=NULL, gene_count=NULL, gene_score_atac=NULL, S1=NULL, S1.genes=NULL){
  ind=which(rownames(gene_count) %in% S1.genes)
  rna_score=colSums(gene_count[ind,])#/colSums(gene_score_atac)
  rna_score=pmin(rna_score, quantile(rna_score[which(Idents(obj) %in% c(S1))], 0.9))
  rna_score=rna_score/max(rna_score)
  
  atac_score=colSums(gene_score_atac[ind,])#/colSums(gene_score
  atac_score=pmin(atac_score, quantile(atac_score[which(Idents(brain) %in% c(S1))], 0.9))
  atac_score=atac_score/max(atac_score)
  
  scores=list(rna_score=rna_score, atac_score=atac_score)
  return(scores)
}


# ComputeScores=function(gene_count=NULL, gene_score_atac=NULL, S1=NULL,ident=NULL, S1.genes=NULL){
#   ind=which(rownames(gene_count) %in% S1.genes)
#   rna_score=colSums(gene_count[ind,])#/colSums(gene_score_atac)
#   rna_score=pmin(rna_score, quantile(rna_score[which(ident %in% c(S1))], 0.9))
#   rna_score=rna_score/max(rna_score)
#   
#   atac_score=colSums(gene_score_atac[ind,])#/colSums(gene_score
#   atac_score=pmin(atac_score, quantile(atac_score[which(ident %in% c(S1))], 0.9))
#   atac_score=atac_score/max(atac_score)
#   
#   scores=list(rna_score=rna_score, atac_score=atac_score)
#   return(scores)
# }







