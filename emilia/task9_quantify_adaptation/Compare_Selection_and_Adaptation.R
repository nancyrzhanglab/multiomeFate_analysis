library(eulerr)

theme_Publication<- function(base_size=12, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            plot.subtitle = element_text(face = "bold", hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            # panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold"),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            # axis.line.x = element_line(colour="black"),
            # axis.line.y = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#F0F0F0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.box.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "mm"),
            legend.direction = "vertical",
            legend.key.size= unit(0.5, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"),
            strip.background=element_rect(colour="#F0F0F0",fill="#F0F0F0"),
            strip.text = element_text(face="bold")
    ))
}

de_d0 <- read.csv('~/Downloads/de.csv', row.names = 1)
de_d10_w5 <- read.csv('~/Downloads/de_d10_w5_adapting_cells.csv', row.names = 1)
# de_d10_w5$logfc <- (-1) * de_d10_w5$logfc

colnames(de_d0) <- c('logfc.d0_high_vs_low_FB', 'p.value.d0_high_vs_low_FB', 'padj.d0_high_vs_low_FB', 'neglog10_pval.d0_high_vs_low_FB')
colnames(de_d10_w5) <- c('logfc.d10_vs_w5', 'p.value.d10_vs_w5', 'padj.d10_vs_w5', 'neglog10_pval.d10_vs_w5')

comp <- merge(de_d0, de_d10_w5, by = 'row.names', all = TRUE)
comp <- comp[comp$padj.d0_high_vs_low_FB < 0.05, ]
comp <- comp[comp$padj.d10_vs_w5 < 0.05, ]
colnames(comp)[1] <- 'gene'
comp$category <- ifelse(comp$gene %in% jackpot, 'Jackpot', 'Other')
comp$category <- ifelse(comp$gene %in% isg.rs, 'ISG.RS', comp$category)
comp <- comp[order(comp$category, decreasing = T), ]

ggplot(comp, aes(x = logfc.d0_high_vs_low_FB, y = logfc.d10_vs_w5)) +
  geom_point(aes(color = category)) +
  scale_color_manual(values = c("blue", "red", "gray")) +
  ggrepel::geom_text_repel(data = subset(comp, category == 'Jackpot'), color = 'red',
                           aes(label = gene), size = 3, max.overlaps = 10) +
  ggrepel::geom_text_repel(data = subset(comp, category == 'ISG.RS'), color = 'blue',
                           aes(label = gene), size = 3, max.overlaps = 10) +
  stat_cor() +
  theme_bw()


de_d0_up <- subset(de_d0, padj.d0_high_vs_low_FB < 0.05 & logfc.d0_high_vs_low_FB > 0.5)
de_d10_w5_up <- subset(de_d10_w5, padj.d10_vs_w5 < 0.05 & logfc.d10_vs_w5 > 0.5)
de_up_intersect <- intersect(rownames(de_d0_up), rownames(de_d10_w5_up))



draw.pairwise.venn(area1=nrow(de_d0_up), area2=nrow(de_d10_w5_up),cross.area=length(de_up_intersect),
                   category=c("Select. d0 (Up)", "Adpatat. d10-to-w5 (Up)"),fill=c("#B6B09F","#FE5D26"))


venn.list <- list(Set1 = rownames(de_d0_up),
                  Set2 = rownames(de_d10_w5_up))
set.seed(12345)
pdf('~/Downloads/venn_selection_adaptation.pdf', width = 3, height = 3)
plot(euler(venn.list), 
     fills = list(fill = c("Set1" = "#FF4F0F",
                           "Set2" =  "#FFE3BB"), 
                  alpha = 0.9),
     labels = list(col = "black", font = 2),
     quantities = list(col = "black", font = 2))
dev.off()

de_d10_w5_up_unique <- setdiff(rownames(de_d10_w5_up), de_up_intersect)

# write de_d10_w5_up_unique to a file
write.csv(de_d10_w5_up_unique, file = '~/Downloads/de_d10_w5_up_unique.csv', row.names = TRUE)

write.csv(de_up_intersect, file = '~/Downloads/de_d10_w5_up_intersect.csv', row.names = TRUE)


de_d0_up_relax <- subset(de_d0, padj.d0_high_vs_low_FB < 0.05 & logfc.d0_high_vs_low_FB > 0.1)
de_up_intersect_relax <- intersect(rownames(de_d0_up_relax), rownames(de_d10_w5_up))


# ==============================================================================
# Read signatures
# ==============================================================================
ref_dir <- '/Users/emiliac/Dropbox/Thesis/resources/GSEA_pathways/'
hallmark <- read.gmt(paste0(ref_dir, "h.all.v2024.1.Hs.symbols.gmt"))
reactome <- read.gmt(paste0(ref_dir, "c2.cp.reactome.v2024.1.Hs.symbols.gmt"))
threeca <- read.gmt(paste0(ref_dir, "c4.3ca.v2024.1.Hs.symbols.gmt"))


de_d10_w5_unique <- de_d10_w5[!rownames(de_d10_w5) %in% de_up_intersect, ]
de_d10_w5 <- de_d10_w5[order(de_d10_w5$logfc.d10_vs_w5, decreasing = TRUE), ]
de_d10_w5$gene <- rownames(de_d10_w5)

# GSEA DABTRAM d0 vs w5
gsea_input.DABTRAM.adaptation <- de_d10_w5$logfc
names(gsea_input.DABTRAM.adaptation) <- de_d10_w5$gene

set.seed(123)
GSEA_res.DABTRAM.adaptation <- GSEA(geneList = gsea_input.DABTRAM.adaptation, 
                                   TERM2GENE = threeca, 
                                   pvalueCutoff = 0.2,
                                   seed = T,
                                   verbose = F)
GSEA_res.DABTRAM.adaption.df <- as_tibble(GSEA_res.DABTRAM.adaptation@result)
GSEA_res.DABTRAM.adaption.df <- GSEA_res.DABTRAM.adaption.df[, c('ID', 'NES', 'p.adjust', 'qvalue', 'setSize')]

GSEA_res.DABTRAM.adaption.df <- GSEA_res.DABTRAM.adaption.df[GSEA_res.DABTRAM.adaption.df$qvalue < 0.05, ]
GSEA_res.DABTRAM.adaption.df$Category <- 'Adaptation'
GSEA_res.DABTRAM.adaption.df <- GSEA_res.DABTRAM.adaption.df[order(GSEA_res.DABTRAM.adaption.df$NES, decreasing = T), ]

# GSEA DABTRAM d0 high vs low adaptation
de_d0 <- de_d0[order(de_d0$logfc.d0_high_vs_low_FB, decreasing = TRUE), ]
de_d0$gene <- rownames(de_d0)
gsea_input.DABTRAM.selection <- de_d0$logfc
names(gsea_input.DABTRAM.selection) <- de_d0$gene

set.seed(123)
GSEA_res.DABTRAM.selection <- GSEA(geneList = gsea_input.DABTRAM.selection, 
                                   TERM2GENE = threeca, 
                                   pvalueCutoff = 1,
                                   seed = T,
                                   verbose = F)
GSEA_res.DABTRAM.selection.df <- as_tibble(GSEA_res.DABTRAM.selection@result)
GSEA_res.DABTRAM.selection.df <- GSEA_res.DABTRAM.selection.df[, c('ID', 'NES', 'p.adjust', 'qvalue', 'setSize')]

GSEA_res.DABTRAM.selection.df <- GSEA_res.DABTRAM.selection.df[GSEA_res.DABTRAM.selection.df$qvalue < 0.05, ]
GSEA_res.DABTRAM.selection.df$Category <- 'Selection'


comp.df <- rbind(GSEA_res.DABTRAM.adaption.df, 
                 GSEA_res.DABTRAM.selection.df)


comp.df <- comp.df[!comp.df$ID %in% c('GAVISH_3CA_MALIGNANT_METAPROGRAM_26_NPC_GLIOMA', 'GAVISH_3CA_MALIGNANT_METAPROGRAM_34_PLATELET_ACTIVATION', 
                                      'GAVISH_3CA_MALIGNANT_METAPROGRAM_16_MES_GLIOMA', 'GAVISH_3CA_MALIGNANT_METAPROGRAM_28_OLIGO_NORMAL',
                                      'GAVISH_3CA_MALIGNANT_METAPROGRAM_40_PDAC_RELATED', 'GAVISH_3CA_MALIGNANT_METAPROGRAM_25_ASTROCYTES', 
                                      'GAVISH_3CA_MALIGNANT_METAPROGRAM_29_NPC_OPC'), ]
comp.df$ID <- gsub('GAVISH_3CA_MALIGNANT_METAPROGRAM_', '', comp.df$ID)
comp.df$ID <- paste0('MP', comp.df$ID)

comp.df <- comp.df[order(comp.df$Category ,comp.df$NES, decreasing = F), ]

comp.df$ID <- gsub('_', ' ', comp.df$ID)
comp.df$ID <- factor(comp.df$ID, levels = unique(comp.df$ID))


ggplot(comp.df, aes(x = Category, y = ID)) +
  geom_point(aes(size = setSize, color = NES)) +
  scale_color_gradient2(low = "#3674B5", mid = 'white', high = "#B22222", midpoint = 0) +
  scale_size_continuous(range = c(1, 7), breaks = c(10, 20, 30)) +
  labs( x = "",
       y = "",
       size = "Set Size",
       color = "Normalized Enrichment Score (NES)") +
  theme_Publication()
ggsave('~/Downloads/Compare_Selection_vs_Adaptation.pdf', width = 8, height = 5)

