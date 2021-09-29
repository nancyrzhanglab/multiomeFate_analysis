rm(list=ls())
nn <- readRDS("../../../../data/SHAREseq_mouseskin/ATAC.RNA.KNN3.rds")
names(nn)

atac <- readRDS("../../../../data/SHAREseq_mouseskin/atac.se.rds")
dim(atac@assays$data$counts)
dim(atac@colData)
colnames(atac@colData)
table(atac@colData$predicted.id)
table(atac$L1)
table(atac$assign)
atac@colData[1:10,1:10]
head(atac@rowRanges)
dim(as.data.frame(atac@rowRanges))

rna <- readRDS("../../../../data/SHAREseq_mouseskin/RNA_smooth.rds")
dim(rna)
rna[1:5,1:5]

# according to chromatin.potential.script.R, this is the DORC file
sesum <- readRDS("../../../../data/SHAREseq_mouseskin/SEsum_smooth.rds")
dim(sesum)
sesum[1:5,1:5]

topic_csv <- read.csv("../../../../data/SHAREseq_mouseskin/topic.matrix.csv",
                      sep = "\t")
dim(topic_csv)
topic_csv[1:5,1:5]
saveRDS(topic_csv, file = "../../../../data/SHAREseq_mouseskin/topic.matrix.rds")

