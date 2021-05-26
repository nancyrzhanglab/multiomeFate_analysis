#### embrain data using share-seq method
#1 generate gene_count and gene_score_atac matrices from seurat linpeaks
brain=readRDS("/Volumes/Cloud/Projects/multiome_fate/brain/data/obj_seurat_links.rds")

rna=as.matrix(brain@assays$RNA@counts)
atac=as.matrix(brain@assays$ATAC@counts)
df=data.frame(gene=brain@assays$ATAC@links$gene, peak=brain@assays$ATAC@links$peak, zscore=brain@assays$ATAC@links$zscore, pvalue=brain@assays$ATAC@links$pvalue, stringsAsFactors = F)
gene_uni=brain@assays$ATAC@links$gene[!duplicated(brain@assays$ATAC@links$gene)]
#rm(brain)
gc()
gene_count=matrix(nrow=length(gene_uni),ncol=ncol(rna))
gene_score_atac=matrix(nrow=length(gene_uni),ncol=ncol(rna))
for(ii in 1:length(gene_uni)){
  gene_count[ii,]=rna[which(rownames(rna)==gene_uni[ii]),]
  gene_score_atac[ii,]=colSums(atac[which(rownames(atac) %in% df$peak[which(df$gene==gene_uni[ii])]),, drop=F])
  print(ii)
}

rownames(gene_count)=gene_uni
rownames(gene_score_atac)=gene_uni
colnames(gene_count)=colnames(rna)
colnames(gene_score_atac)=colnames(rna)

rm(atac) 
gc()
saveRDS(gene_count,"/Volumes/Cloud/Projects/multiome_fate/brain/data/gene_count.rds")
saveRDS(gene_score_atac,"/Volumes/Cloud/Projects/multiome_fate/brain/data/gene_score_atac.rds")

#2 generate gene_count and gene_score_atac matrices from Nancy's peak selection method
gene_score_atac=t(myobj$Yhat.cell$regress.significant)
gene_count=t(myobj$Y.cell)

rownames(gene_count)=myobj$gene.names
rownames(gene_score_atac)=myobj$gene.names
#colnames(gene_count)=colnames(rna)
colnames(gene_score_atac)=colnames(gene_count)


genes.hockey=c("Tnc", "Slit3", "Npas3", "Npy", "Foxp2",
               "Tenm2", "Top2a", "Sox6","Opcml","Apoe", "Sparc",
               "Fabp7", "Ntng1", "Mki67")
gene_score_atac=gene_score_atac[which(rownames(gene_score_atac) %in% genes.hockey),]
gene_count=gene_count[which(rownames(gene_count) %in% genes.hockey),]


sel_gene=which(apply(gene_score_atac, 1, function(x) !any(is.na(x))))
gene_score_atac=gene_score_atac[sel_gene,]
gene_count=gene_count[sel_gene,]

saveRDS(gene_count,"/Volumes/Cloud/Projects/multiome_fate/brain/data/gene_count_sum.sig.hockey.rds")
saveRDS(gene_score_atac,"/Volumes/Cloud/Projects/multiome_fate/brain/data/gene_score_atac_sum.sig.regress.hockey.rds")



# nearest neighbors based on ATAC UMAP embedding from scVelo
library(FNN)
k=15
umap_coor=read.table("~/cellrank/umap_coor.txt", stringsAsFactors = F, sep='\t', header=T)
umap.old <-  data.frame(umap1=umap_coor$X0, umap2=umap_coor$X1)
knn.smooth = get.knn(umap.old, k = k)
knn.smooth <- knn.smooth$nn.index


# smooth the two matrix using wnn nearest neighbors
gene_counts=readRDS("~/gene_count.rds")
gene_score_atac=readRDS("~/gene_score_atac..rds")

gene_count_smooth=matrix(nrow=nrow(gene_count), ncol=ncol(gene_count))
gene_score_atac_smooth=matrix(nrow=nrow(gene_count), ncol=ncol(gene_count))
for(ii in 1:ncol(gene_count_smooth)){
  gene_count_smooth[,ii]=rowSums(gene_count[,c(ii, knn.smooth[ii,]), drop=F])
  gene_score_atac_smooth[,ii]=rowSums(gene_score_atac[,c(ii, knn.smooth[ii,]), drop=F])
  print(ii)
}

colnames(gene_count_smooth)=colnames(gene_count)
colnames(gene_score_atac_smooth)=colnames(gene_score_atac)
rownames(gene_count_smooth)=rownames(gene_count)
rownames(gene_score_atac_smooth)=rownames(gene_score_atac)

# find high var genes
x <- log10(rowMeans(gene_score_atac_smooth))
std <- apply(gene_score_atac_smooth, 1, sd)
y = std^2/rowMeans(gene_score_atac_smooth) 
plot(x, y, ylim=c(0,60))
genes <- names(y[y > 10]) 
length(genes) # 429 genes are selected

#celltype=as.character(atac.se$assign[which(cellindex)])
celltype=readRDS("~/data_tenx_labels_Jane.rds")
celltype0=celltype$savercatLable
celltype=celltype0[match(colnames(gene_count_smooth),names(celltype0))]

## knn.index.dist
library(KernelKnn)
# this step is very slow
ATAC.RNA.KNN3 <- knn.index.dist(t(gene_score_atac_smooth[genes,]),t(gene_count_smooth[genes,]), method="pearson_correlation", k=10, threads = 3) # pearson_correlation
dim(ATAC.RNA.KNN3$test_knn_idx)
saveRDS(ATAC.RNA.KNN3, "~/ATAC.RNA.KNN3.umap.rds")
# ATAC.RNA.KNN3 <- readRDS('ATAC.RNA.KNN3.rds')
# smooth in ATAC umap with 10 KNNs
umap_coor=read.table("~/cellrank/umap_coor.txt", stringsAsFactors = F, sep='\t', header=T)
umap.old <-  data.frame(umap1=umap_coor$X0, umap2=umap_coor$X1)
umap.new <- umap.old
umap.new$mean.dist <- 0
for (i in 1:nrow(umap.new)){
  umap.new[i, 1:2] <- colMeans(umap.old[ATAC.RNA.KNN3$test_knn_idx[i, ], ])
  umap.new$mean.dist[i] <- max(dist(umap.old[ATAC.RNA.KNN3$test_knn_idx[i, ], ],method = "euclidean"))
}

#pdf("./skin/plots/umap_old_new_col.pdf")
plot(umap.old[, 1:2], col=as.factor(celltype),pch=20)
plot(umap.new[, 1:2], col=as.factor(celltype), pch=20)
#dev.off()

# remove ORS cells for visulization
temp <- data.frame(umap.old, umap.new)

# filter out top 5% long arrow
temp$arrow.legnth <- ((temp$umap1.1- temp$umap1)^2 + (temp$umap2.1- temp$umap2)^2)^0.5
hist(temp$arrow.legnth, breaks = 100)
temp$umap1.1[temp$arrow.legnth > quantile(temp$arrow.legnth, probs=0.95)] <- temp$umap1[temp$arrow.legnth > quantile(temp$arrow.legnth, probs=0.95)]
temp$umap2.1[temp$arrow.legnth > quantile(temp$arrow.legnth, probs=0.95)] <- temp$umap2[temp$arrow.legnth > quantile(temp$arrow.legnth, probs=0.95)]

## smooth in umap space
umap.raw <- temp; head(umap.raw)
pcs <- data.frame(umap1=temp$umap1, umap2=temp$umap2)
k = 15 
library(FNN)
knn.norm = get.knn(pcs, k = k)
knn.norm <- knn.norm$nn.index

umap.smooth <- umap.raw
for (i in 1:nrow(umap.smooth)){
  umap.smooth[i,3] <- mean(umap.raw[knn.norm[i, ],3])
  umap.smooth[i,4] <- mean(umap.raw[knn.norm[i, ],4])
}


umap.smooth$arrow.legnth <- ((umap.smooth$umap1.1- umap.smooth$umap1)^2 + (umap.smooth$umap2.1- umap.smooth$umap2)^2)^0.5

chromatin.potential <- umap.smooth 
saveRDS(chromatin.potential, "~/chromatin.potential.rds")


library(RColorBrewer)
library(metR)
library(RColorBrewer)
library(BuenColors)

scale.factor = 0.1
pp=ggplot(umap.smooth, aes(x=umap1, y=umap2)) +
  geom_point(aes(color=arrow.legnth), size=1, stroke=0, alpha = 1) +
  geom_arrow(aes(dx = (umap1.1-umap1),dy = (umap2.1-umap2)), skip = 1,
             color = "black",lineend="round", size=0.08, arrow.length = 0.5) +
  scale_mag(max_size = 0.1, guide = "none") +
  scale_color_gradientn(colors = rev(colorRampPalette(brewer.pal(11,'Spectral')[-6])(100))) +
  pretty_plot(fontsize = 5) + L_border() + 
  guides(colour = FALSE)

pp

#ggsave("./skin/plots/my.chromatin.ptime_rnrmarrow.png", width=3.8, height=3, units = "cm", dpi=1200)



