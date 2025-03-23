rm(list = ls())

library(Biobase)
library(GSVA)
library(survival)
library(randomForestSRC)
library(ggplot2)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data_other/Hatzis_et_al_breast_cancer_clinical_microarray_GSE25065/'
result_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task8_validate_week5_predictors_with_clinical_data/'
figure_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/figures/MultiomeFate_V3/Fig5/'


irds <- c("FLNB", "GLUD1", "IFI44", "IFIT1", "IFIT3", "IGF1R", "ISG15", 
          "MX1", "OAS1", "STAT1", "UBE2D3", "ZNF273", "CTDSPL", "ABCD3")

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
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold"),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line.x = element_line(colour="black"),
            axis.line.y = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="white"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.box.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "mm"),
            legend.direction = "horizontal",
            legend.key.size= unit(0.5, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"),
            strip.background=element_rect(colour="#F0F0F0",fill="#F0F0F0"),
            strip.text = element_text(face="bold")
    ))
}

use.pal <- c('#4285F4', '#C00000')
names(use.pal) <- c('pCR', 'RD')
# ==============================================================================
# Load data
# ==============================================================================
load(paste0(data_dir, 'MDACC JAMA expression set.RData'))
load(paste0(data_dir, 'MDACC JAMA rma expression set.RData'))
load(paste0(data_dir, "Prognostic signature probes and centroids.RData"))

# ==============================================================================
# Load helpers
# ==============================================================================
source(paste0(data_dir, "Probes to genes functions.R"))
source(paste0(data_dir, "Functions to estimate IRDS and prognostic signatures.R"))
source(paste0(data_dir, "RSF-RF functions.R"))

# ==============================================================================
# Process data and assemble metadata
# ==============================================================================
eset <- MDACCrmaSet

## Collapse probes to genes
unqGSmc <- calc.probesTOgenes(exprs(eset), gs=as.character(fData(eset)$Symbol))

## Calc RS status
RS <- calcRS(eset, RS.probes, probe.id="affy.probeID")$RS$RSclass

## Calc IRDS score
irds.edat <- medCenter(unqGSmc[irds, ])
tsp.irds <- tspIRDS(irds.edat)

## Prepare and clean-up data frame
cen <- pData(eset)$drfs_1_event_0_censored
cen <- ifelse(is.na(cen), 0, cen)
rfs <- pData(eset)$drfs_even_time_years
rfs <- ifelse(rfs <= 0, 0.001, rfs)

use.index <- c(4:7, 9:10, 12:13, 21) # clinical variables to extract
rsf.dat <- data.frame(cen=cen, rfs=rfs, pData(eset)[, use.index], RS=RS)
rsf.dat[, "grade"] <- as.factor(ifelse(rsf.dat$grade=="4=Indeterminate", "Indeterminate", rsf.dat$grade))

# ==============================================================================
# Gene signatures
# ==============================================================================
genes.all <- rownames(unqGSmc)

ref_dir <- '/Users/emiliac/Dropbox/Epi Evolution Paper/Data/Gavish Clinical Analysis/Resources/'
ref_dir2 <- '/Users/emiliac/Dropbox/Epi Evolution Paper/Data/Persistent IFN ISG Groups/'

# ISG.RS
isg.rs <- read.table(paste0(ref_dir, 'ISG.RS.txt'), header=FALSE)
colnames(isg.rs) <- c("Gene")
isg.rs <- unique(na.omit(isg.rs$Gene[isg.rs$Gene %in% genes.all]))

# ISG.mem
isg.mem <- read.csv(paste0(ref_dir2, 'Memory ISGs Human.csv'), header = FALSE)
colnames(isg.mem) <- c("Gene")
isg.mem <- unique(na.omit(isg.mem$Gene[isg.mem$Gene %in% genes.all]))

# Gavish MP
gavish.mp <- read.csv(paste0(ref_dir, 'Gavish_Malignant_Meta_Programs.csv'))
gavish.mp.list <- list()
for (mp in colnames(gavish.mp)) {
  gavish.mp.list[[mp]] <- c(gavish.mp[[mp]] %>% as.character())
}

gavish.mp.list <- lapply(1:length(gavish.mp.list), function(x) {
  gs <- gavish.mp.list[[x]]
  gs <- unique(na.omit(gs[gs %in% genes.all]))
  return(gs)
})
names(gavish.mp.list) <- colnames(gavish.mp)

# Custom signatures
custom.cis <- c('NUCKS1','TUBB','HMGB2', 'ICMT', 'CBX5', 'TUBA1B', 'ANP32B','TYMS',
                'GMNN', 'USP1', 'NASP', 'TMPO', 'NCAPH', 'TK1', 'TUBG1', 'PRC1',
                'PBK', 'SMC3', 'RRM2', 'RAD51AP1')
custom.dabtram <- c('ACTB', 'TMEM43', 'TPM4', 'CALM2', 'FN1', 'PALLD', 'LMO7',
                    'ACTN1', 'HSPG2', 'MYOF', 'TNFRSF12A', 'TUBB', 'RCN1', 'CRIM1',
                    'COL5A2', 'SAMD5', 'TPM1', 'OXSR1', 'CBX5')
custom.cocl2 <- c('GXYLT2', 'ANTXR1', 'CADM1', 'ITGB3', 'BICC1', 'SLC1A4',
                  'CADPS', 'HMGA2', 'TIMP3', 'PTPRG', 'SERPINE2', 'IMMP2L', 'LRMDA',
                  'MFSD12', 'SOX5', 'EPHA3', 'PRKG2', 'IL1RAP', 'SLC44A1', 'KCNQ5')


custom.cis.acute <- c('FRMD4B', 'GHR', 'MYO1D', 'TDRD3', 'CDYL', 'ACBD6','SLCO5A1',
                      'ADCY2', 'HNRNPA2B1', 'PLCB4', 'SEMA5A', 'PDE3B', 'PMEL', 'SLC7A2',
                      'B4GALT5', 'BNC2', 'GPM6B', 'PRKG2', 'IL16', 'NEDD4L')

# custom.acute <- c('ADCY2', 'MYO1D', 'SLCO5A1', 'HIBCH', 'GHR', 'PDE3B', 'PLCB4', 'STK32A',
#                   'TMEM163', 'BNC2', 'NEDD4L', 'IGSF11', 'RAB38', 'SLC39A11', 'TDRD3',
#                   'PRKG2', 'HMCN1', 'GRK3', 'SEMA6A', 'UBL3', 'PIK3R3', 'IL16', 'NRG3',
#                   'TRIM2')

reported.pcr <- c('AGA', 'ANGEL1', 'ANKRD11', 'ARID3A', 'ASTN2', 'BID', 'BPI', 'CCHCR1',
                  'CGREF1', 'COMT', 'DCAF8', 'DLGAP4', 'FTSJ3', 'GAD2', 'GCNT1', 'GOSR1',
                  'HMGXB4', 'KIAA0406', 'KIF24', 'KLF11', 'KPNA1', 'LMO4', 'LPCAT1', 'MEF2A', 'MYCL1',
                  'NAGA', 'PAK1IP1', 'PI4KA', 'PPARA')

custom.list <- list(custom.dabtram, custom.cocl2, custom.cis, custom.cis.acute, reported.pcr)
custom.list <- lapply(1:length(custom.list), function(x) {
  gs <- custom.list[[x]]
  gs <- unique(na.omit(gs[gs %in% genes.all]))
  return(gs)
})
names(custom.list) <- c('custom.dabtram', 'custom.cocl2', 'custom.cis', 'custom.cis.acute', 'reported.pcr')


gene_list <- c(list(isg.rs), list(isg.mem), custom.list, gavish.mp.list)
names(gene_list) <- c('ISG.RS', 'ISGmem', names(custom.list), names(gavish.mp.list))

# ==============================================================================
# Run GSVA
# ==============================================================================
gsvaPar <- gsvaParam(unqGSmc, gene_list)
gsva.es <- gsva(gsvaPar, verbose=FALSE)

rownames(gsva.es) <- names(gene_list)

# ==============================================================================
# Calculate relapse risk and plot boxplots
# ==============================================================================

ntree <- 1000
nodesize <- 1
nsplit <- 2

ntree.null    <- c(10, 100, 1000, 10000)[3]
risk.perc <- c(0.1, 0.8)  # low and hi risk mortality cut-offs for RF-RSF ensemble survival plots


## Divide data into training and test sets based on original publication
##----------------------------------------------------------------------------
traindat <- which(pData(eset)$train.test=="discovery")

## Multivariable analysis using RSF for relapse and RF for pCR
##------------------------------------------------------------------
rsf.dat$RS <- as.factor(rsf.dat$RS)
rfsrc.f <- as.formula(Surv(rfs, cen) ~ .)
rfit <- rfsrc(rfsrc.f, ntree=ntree, nsplit=nsplit, nodesize=nodesize, forest=TRUE, 
              split.depth="by.tree", data=rsf.dat[traindat, ], na.action = "na.impute")


cp.dat <- data.frame(oob.surv=rfit$predicted.oob/rfit$n, rfit$xvar)


# ==============================================================================
# DE between pCR and RD
# ==============================================================================
gex.mat <- unqGSmc

pt.pCR <- cp.dat$Sample[which(cp.dat$pathologic_response_pcr_rd == "pCR")]
pt.RD <- cp.dat$Sample[which(cp.dat$pathologic_response_pcr_rd == "RD")]

pt.pCR <- intersect(pt.pCR, colnames(gex.mat))
pt.RD <- intersect(pt.RD, colnames(gex.mat))

pvalue_list <- lapply(rownames(gex.mat), function(gene){
  x_vec <- gex.mat[gene, pt.pCR]
  y_vec <- gex.mat[gene, pt.RD]
  
  if(diff(range(x_vec)) <= 1e-6 || diff(range(y_vec)) <= 1e-6 ) return(NA)
  
  test_res <- stats::wilcox.test(x = x_vec, y = y_vec)
  list(stat = mean(x_vec) - mean(y_vec), pvalue = test_res$p.value)
})

names(pvalue_list) <- rownames(gex.mat)
pvalue_list <- pvalue_list[sapply(pvalue_list, function(x){!all(is.na(x))})]
pvalue_vec <- sapply(pvalue_list, function(x){x$pvalue})
names(pvalue_vec) <- names(pvalue_list)
pvalue_vec <- stats::p.adjust(pvalue_vec, method = "BH")

diff_vec <- sapply(pvalue_list, function(x){x$stat})
logpval_vec <- sapply(pvalue_list, function(x){-log10(x$pvalue)})

df <- data.frame(difference = diff_vec,
                 log10pval = logpval_vec,
                 name = names(logpval_vec))

result_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/'

cor.dabtram <- read.csv(paste0(result_dir, 'dabtram_saver_cor_d0_d10.csv'))
# cor.cocl2 <- read.csv(paste0(result_dir, 'cocl2_saver_cor_d0_d10.csv'))
cor.cis <- read.csv(paste0(result_dir, 'cis_saver_cor_d0_d10.csv')) 

df <- merge(df, cor.cis[, c('gene', 'correlation.CIS_d0')], by.x='name', by.y = 'gene')
df <- merge(df, cor.dabtram[, c('gene', 'correlation.DABTRAM_d0')], by.x='name', by.y = 'gene')
df.2 <- df[df$name %in% set$gene, ]

ggplot(df, aes(x = difference, y = correlation.CIS_d0)) +
  geom_point( alpha = 0.5) +
  geom_point(data = df.2, aes(x = difference, y = correlation.CIS_d0), color = 'red') +
  theme_minimal() +
  stat_cor() +
  xlab('Difference in expression between pCR and RD') +
  ylab('Correlation with CIS fate potential') +
  ggtitle('CIS fate potential correlation with DE between pCR and RD') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed')

# ==============================================================================
# Validate.1 
# ==============================================================================
rsf.dat.use <- rsf.dat[, c('rfs', 'cen', 'pathologic_response_pcr_rd', 'RS')]
rsf.dat.use$Sample <- rownames(rsf.dat.use)
rsf.dat.use$cen <- as.factor(rsf.dat.use$cen)

gsva.es.df <- as.data.frame(t(gsva.es))
gsva.es.df$Sample <- rownames(gsva.es.df)

gsva.es.df <- merge(gsva.es.df, rsf.dat.use, by='Sample')

cp.dat$Sample <- rownames(cp.dat)
cp.dat <- merge(cp.dat, gsva.es.df, by=c("Sample", 'pathologic_response_pcr_rd', 'RS'), all.x=TRUE)
cp.dat <- cp.dat %>% 
  filter(pathologic_response_pcr_rd %in% c('pCR', 'RD'))



# ==============================================================================
# Validate.2
# ==============================================================================

cp.dat$relapse.risk.hi <- ifelse(cp.dat$oob.surv > median(cp.dat$oob.surv, na.rm=TRUE), "High", "Low")
cp.dat$relapse.risk.hi <- factor(cp.dat$relapse.risk.hi, levels=c("Low", "High"))

cp.dat$isg.rs.hi <- ifelse(cp.dat$ISG.RS > median(cp.dat$ISG.RS, na.rm=TRUE), "High", "Low")
cp.dat$isg.rs.hi <- factor(cp.dat$isg.rs.hi, levels=c("Low", "High"))

cp.dat$dabtram.hi <- ifelse(cp.dat$custom.dabtram > median(cp.dat$custom.dabtram, na.rm=TRUE), "High", "Low")
cp.dat$dabtram.hi <- factor(cp.dat$dabtram.hi, levels=c("Low", "High"))

cp.dat$cocl2.hi <- ifelse(cp.dat$custom.cocl2 > median(cp.dat$custom.cocl2, na.rm=TRUE), "High", "Low")
cp.dat$cocl2.hi <- factor(cp.dat$cocl2.hi, levels=c("Low", "High"))

cp.dat$cis.hi <- ifelse(cp.dat$custom.cis > median(cp.dat$custom.cis, na.rm=TRUE), "High", "Low")
cp.dat$cis.hi <- factor(cp.dat$cis.hi, levels=c("Low", "High"))

# Cisplatin long-term
p1.cis <- ggplot(cp.dat, aes(x=pathologic_response_pcr_rd, y=custom.cis)) + 
  geom_violin(aes(fill = pathologic_response_pcr_rd), scale = 'width', alpha = 0.2) +
  geom_jitter(aes(color = pathologic_response_pcr_rd), width = 0.1, size = 0.5) +
  geom_boxplot(outlier.shape = NA, width = 0.2, alpha = 0.8, lwd = 0.8) + 
  stat_compare_means(comparisons = list(c('pCR', 'RD')),
                     paired = FALSE, 
                     method = "wilcox", aes(label =  after_stat(p.signif))) +
  scale_color_manual(values = use.pal) +
  scale_fill_manual(values = use.pal) +
  labs(x = "Pathologic Response", 
       y = "CIS w5 top20 predictors",
       title = "CIS and Pathologic Response") +
  theme_Publication() + 
  theme(legend.position="bottom")

p2.cis <- ggplot(cp.dat, aes(x=custom.cis, y=oob.surv)) + 
  geom_point(color = 'gray') +
  stat_cor() +
  facet_wrap(~pathologic_response_pcr_rd) +
  xlab("CIS w5 top20 predictors") + 
  ylab('Relapse Risk') +
  ggtitle('CIS w5 top20 predictors (GSVA) vs Relapse Risk') +
  theme_bw()


p3.cis <- ggplot(cp.dat, aes(x=relapse.risk.hi, y=custom.cis)) + 
  geom_violin(aes(fill = pathologic_response_pcr_rd), scale = 'width', alpha = 0.2) +
  geom_jitter(aes(color = pathologic_response_pcr_rd), width = 0.1, size = 0.5) +
  geom_boxplot(outlier.shape = NA, width = 0.2, alpha = 0.8, lwd = 0.8) + 
  stat_compare_means(comparisons = list(c('Low', 'High')),
                     paired = FALSE, 
                     method = "wilcox", aes(label =  after_stat(p.signif))) +
  facet_wrap(. ~ pathologic_response_pcr_rd) +
  scale_color_manual(values = use.pal) +
  scale_fill_manual(values = use.pal) +
  labs(x = "Relapse Risk Status", 
       y = "CIS w5 top20 predictors",
       title = "CIS and Relapse Risk") +
  theme_Publication() + 
  theme(legend.position="bottom")

p.cis <- ggarrange(
  p1.cis, p3.cis,  widths = c(0.6, 1),
  common.legend = TRUE, legend = "bottom")

ggsave(p.cis, filename = paste0(figure_dir, "Fig5H.cisplatin_d10ToWeek5Top20Genes.pdf"), width = 7, height = 4)

# COCL2 long term
p1.cocl2 <- ggplot(cp.dat, aes(x=pathologic_response_pcr_rd, y=custom.cocl2)) + 
  geom_violin(aes(fill = pathologic_response_pcr_rd), scale = 'width', alpha = 0.2) +
  geom_jitter(aes(color = pathologic_response_pcr_rd), width = 0.1, size = 0.5) +
  geom_boxplot(outlier.shape = NA, width = 0.2, alpha = 0.8, lwd = 0.8) + 
  stat_compare_means(comparisons = list(c('pCR', 'RD')),
                     paired = FALSE, 
                     method = "wilcox", aes(label =  after_stat(p.signif))) +
  scale_color_manual(values = use.pal) +
  scale_fill_manual(values = use.pal) +
  labs(x = "Pathologic Response", 
       y = "COCL2 w5 top20 predictors",
       title = "COCL2 and Pathologic Response") +
  theme_Publication() + 
  theme(legend.position="bottom")

p3.cocl2 <- ggplot(cp.dat, aes(x=relapse.risk.hi, y=custom.cocl2)) + 
  geom_violin(aes(fill = pathologic_response_pcr_rd), scale = 'width', alpha = 0.2) +
  geom_jitter(aes(color = pathologic_response_pcr_rd), width = 0.1, size = 0.5) +
  geom_boxplot(outlier.shape = NA, width = 0.2, alpha = 0.8, lwd = 0.8) + 
  stat_compare_means(comparisons = list(c('Low', 'High')),
                     paired = FALSE, 
                     method = "wilcox", aes(label =  after_stat(p.signif))) +
  facet_wrap(. ~ pathologic_response_pcr_rd) +
  scale_color_manual(values = use.pal) +
  scale_fill_manual(values = use.pal) +
  labs(x = "Relapse Risk Status", 
       y = "COCL2 w5 top20 predictors",
       title = "COCL2 and Relapse Risk") +
  theme_Publication() + 
  theme(legend.position="bottom")

p.cocl2 <- ggarrange(
  p1.cocl2, p3.cocl2,  widths = c(0.6, 1),
  common.legend = TRUE, legend = "bottom")
ggsave(p.cocl2, filename = paste0(figure_dir, "Fig5H.cocl2_d10ToWeek5Top20Genes.pdf"), width = 7, height = 4)


ggplot(cp.dat, aes(x=pathologic_response_pcr_rd, y=reported.pcr)) + 
  geom_violin(aes(fill = pathologic_response_pcr_rd), scale = 'width', alpha = 0.2) +
  geom_jitter(aes(color = pathologic_response_pcr_rd), width = 0.1, size = 0.5) +
  geom_boxplot(outlier.shape = NA, width = 0.2, alpha = 0.8, lwd = 0.8) + 
  stat_compare_means(comparisons = list(c('pCR', 'RD')),
                     paired = FALSE, 
                     method = "wilcox", aes(label =  after_stat(p.signif))) +
  scale_color_manual(values = use.pal) +
  scale_fill_manual(values = use.pal) +
  labs(x = "Pathologic Response", 
       y = "COCL2 w5 top20 predictors",
       title = "COCL2 and Pathologic Response") +
  theme_Publication() + 
  theme(legend.position="bottom")
