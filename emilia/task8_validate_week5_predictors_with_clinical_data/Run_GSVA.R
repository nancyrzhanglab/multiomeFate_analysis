rm(list = ls())

library(Biobase)
library(GSVA)
library(survival)
library(randomForestSRC)
library(ggplot2)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data_other/Hatzis_et_al_breast_cancer_clinical_microarray_GSE25065/'

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


# custom.cocl2 <- c('MT-ND4L', 'MT-ND3', 'MT-ATP8', 'MT-CYB', 'ZNF704', 'MT-ND5',
#                   'CADPS', 'HMGA2', 'TIMP3', 'PTPRG', 'SERPINE2', 'MT-ND2', 'MDGA2',
#                   'FMN1', 'MT-ATP6', 'ROBO2', 'MT-ND1', 'IL1RAP', 'SLC44A1', 'MTRNR2L10')
custom.cis.acute <- c('FRMD4B', 'GHR', 'MYO1D', 'TDRD3', 'CDYL', 'ACBD6','SLCO5A1',
                      'ADCY2', 'HNRNPA2B1', 'PLCB4', 'SEMA5A', 'PDE3B', 'PMEL', 'SLC7A2',
                      'B4GALT5', 'BNC2', 'GPM6B', 'PRKG2', 'IL16', 'NEDD4L')

custom.list <- list(custom.dabtram, custom.cocl2, custom.cis, custom.cis.acute)
custom.list <- lapply(1:length(custom.list), function(x) {
  gs <- custom.list[[x]]
  gs <- unique(na.omit(gs[gs %in% genes.all]))
  return(gs)
})
names(custom.list) <- c('custom.dabtram', 'custom.cocl2', 'custom.cis', 'custom.cis.acute')


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

ggplot(cp.dat, aes(x = rfs, y = oob.surv)) + 
  geom_point(aes(color = cen)) + 
  geom_smooth(method='lm') + 
  stat_cor() +
  xlab('Relapse free survival (years)') +
  ylab('Relapse Risk') +
  theme_Publication()

ggplot(cp.dat, aes(x = rfs, y = oob.surv)) + 
  geom_point(aes(color = cen)) + 
  geom_smooth(method='lm') + 
  stat_cor() +
  facet_wrap(. ~ pathologic_response_pcr_rd) +
  xlab('Relapse free survival (years)') +
  ylab('Relapse Risk') +
  theme_Publication()



# ggplot(gsva.es.df, aes(x=rfs, y=Proteasomal.degradation)) + 
#   geom_point() + 
#   geom_smooth(method='lm') + 
#   theme_bw()
# 
# ggplot(gsva.es.df, aes(x=pathologic_response_pcr_rd, y=Proteasomal.degradation)) + 
#   geom_boxplot() + 
#   facet_wrap(~RS) +
#   theme_bw()
# 
# ggplot(gsva.es.df, aes(x=RS, y=Proteasomal.degradation)) + 
#   geom_boxplot() + 
#   facet_wrap(~pathologic_response_pcr_rd) +
#   theme_bw()

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

cp.dat$cis.acute.hi <- ifelse(cp.dat$custom.cis.acute > median(cp.dat$custom.cis.acute, na.rm=TRUE), "High", "Low")
cp.dat$cis.acute.hi <- factor(cp.dat$cis.acute.hi, levels=c("Low", "High"))

p1 <- ggplot(cp.dat, aes(x=pathologic_response_pcr_rd, y=ISG.RS)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, color = 'gray') +
  theme_bw() + 
  xlab("pathologic_response_pcr_rd") + 
  ylab("ISG.RS") + 
  ggtitle("ISG.RS and Pathologic Response") + theme(legend.position="bottom")

p2 <- ggplot(cp.dat, aes(x=ISG.RS, y=oob.surv)) + 
  geom_point(color = 'gray') +
  stat_cor() +
  facet_wrap(~pathologic_response_pcr_rd) +
  ylab('Relapse Risk') +
  ggtitle('ISG.RS (GSVA) vs Relapse Risk') +
  theme_bw()

p3 <- ggplot(cp.dat, aes(x=isg.rs.hi, y=oob.surv)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, color = 'gray') +
  facet_wrap(~pathologic_response_pcr_rd) +
  theme_bw() + 
  xlab("ISG.RS status") + ylab("Relapse Risk") + 
  ggtitle("ISG.RS and Relapse Risk") + 
  theme(legend.position="bottom")

p3 <- ggplot(cp.dat, aes(x=relapse.risk.hi, y=ISG.RS)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, color = 'gray') +
  facet_wrap(~pathologic_response_pcr_rd) +
  theme_bw() + 
  xlab("Relapse Risk status") + ylab("ISG.RS") + 
  ggtitle("ISG.RS and Relapse Risk") + 
  theme(legend.position="bottom")

ggarrange(
  ggarrange(
  p1, p3, labels = c("A", "B"), widths = c(0.6, 1),
  common.legend = TRUE, legend = "bottom"), 
  p2, labels = c("", "C"), ncol = 1)

ggarrange(
  p1, p3, labels = c("A", "B"), widths = c(1, 0.8),
  common.legend = TRUE, legend = "bottom")

wilcox.test(cp.dat$ISG.RS ~ cp.dat$pathologic_response_pcr_rd)
wilcox.test(cp.dat$ISG.RS ~ cp.dat$relapse.risk.hi)


p1.dabtram <- ggplot(cp.dat, aes(x=pathologic_response_pcr_rd, y=custom.dabtram)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, color = 'gray') +
  theme_bw() + 
  xlab("pathologic_response_pcr_rd") + 
  ylab("DABTRAM w5 top20 predictors") + 
  ggtitle("DABTRAM and Pathologic Response") + theme(legend.position="bottom")

p2.dabtram <- ggplot(cp.dat, aes(x=custom.dabtram, y=oob.surv)) + 
  geom_point(color = 'gray') +
  stat_cor() +
  facet_wrap(~pathologic_response_pcr_rd) +
  xlab("DABTRAM w5 top20 predictors") + 
  ylab('Relapse Risk') +
  ggtitle('DABTRAM w5 top20 predictors (GSVA) vs Relapse Risk') +
  theme_bw()

p3.dabtram <- ggplot(cp.dat, aes(x=dabtram.hi, y=oob.surv)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, color = 'gray') +
  facet_wrap(~pathologic_response_pcr_rd) +
  theme_bw() + 
  xlab("DABTRAM top20 predictive genes D10 to W5") + ylab("Relapse Risk") + 
  ggtitle("custom.DABTRAM and Relapse Risk") + theme(legend.position="bottom")

p3.dabtram <- ggplot(cp.dat, aes(x=relapse.risk.hi, y=custom.dabtram)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, color = 'gray') +
  facet_wrap(~pathologic_response_pcr_rd) +
  theme_bw() + 
  xlab("Relapse Risk Status") + ylab("DABTRAM w5 top20 predictors") + 
  ggtitle("custom.DABTRAM and Relapse Risk") + theme(legend.position="bottom")

ggarrange(
  ggarrange(
    p1.dabtram, p3.dabtram, labels = c("A", "B"), widths = c(0.6, 1),
    common.legend = TRUE, legend = "bottom"), 
  p2.dabtram, labels = c("", "C"), ncol = 1)

cp.dat.pCR <- cp.dat[cp.dat$pathologic_response_pcr_rd == 'RD', ]
cp.dat.pCR.dabtram.hi <- cp.dat.pCR[cp.dat.pCR$dabtram.hi == 'High', ]
cp.dat.pCR.dabtram.lo <- cp.dat.pCR[cp.dat.pCR$dabtram.hi == 'Low', ]
wilcox.test(cp.dat.pCR.dabtram.hi$oob.surv, cp.dat.pCR.dabtram.lo$oob.surv)




p1.cocl2 <- ggplot(cp.dat, aes(x=pathologic_response_pcr_rd, y=custom.cocl2)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, color = 'gray') +
  theme_bw() + 
  xlab("pathologic_response_pcr_rd") + 
  ylab("COCL2 w5 top20 predictors") + 
  ggtitle("COCL2 and Pathologic Response") + theme(legend.position="bottom")

p2.cocl2 <- ggplot(cp.dat, aes(x=custom.cocl2, y=oob.surv)) + 
  geom_point(color = 'gray') +
  stat_cor() +
  facet_wrap(~pathologic_response_pcr_rd) +
  xlab("COCL2 w5 top20 predictors") + 
  ylab('Relapse Risk') +
  ggtitle('COCL2 w5 top20 predictors (GSVA) vs Relapse Risk') +
  theme_bw()

p3.cocl2 <- ggplot(cp.dat, aes(x=cocl2.hi, y=oob.surv)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, color = 'gray') +
  facet_wrap(~pathologic_response_pcr_rd) +
  theme_bw() + 
  xlab("COCL2 top20 predictive genes D10 to W5") + ylab("Relapse Risk") + 
  ggtitle("custom.COCL2 and Relapse Risk") + theme(legend.position="bottom")

p3.cocl2 <- ggplot(cp.dat, aes(x=relapse.risk.hi, y=custom.cocl2)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, color = 'gray') +
  facet_wrap(~pathologic_response_pcr_rd) +
  theme_bw() + 
  xlab("Relapse Risk Status") + ylab("COCL2 w5 top20 predictors") + 
  ggtitle("custom.COCL2 and Relapse Risk") + theme(legend.position="bottom")

ggarrange(
  ggarrange(
    p1.cocl2, p3.cocl2, labels = c("A", "B"), widths = c(0.6, 1),
    common.legend = TRUE, legend = "bottom"), 
  p2.cocl2, labels = c("", "C"), ncol = 1)

ggplot(cp.dat, aes(x=cocl2.hi, y=oob.surv)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, color = 'gray') +
  facet_wrap(~pathologic_response_pcr_rd) +
  theme_bw() +  
  xlab("COCL2 top20 predictive genes D10 to W5") + ylab("Relapse Risk") + 
  ggtitle("custom COCL2 and Relapse Risk") + theme(legend.position="bottom")
cp.dat.pCR <- cp.dat[cp.dat$pathologic_response_pcr_rd == 'pCR', ]
cp.dat.pCR.cocl2.hi <- cp.dat.pCR[cp.dat.pCR$cocl2.hi == 'High', ]
cp.dat.pCR.cocl2.lo <- cp.dat.pCR[cp.dat.pCR$cocl2.hi == 'Low', ]
wilcox.test(cp.dat.pCR.cocl2.hi$oob.surv, cp.dat.pCR.cocl2.lo$oob.surv)


# Cisplatin long-term
p1.cis <- ggplot(cp.dat, aes(x=pathologic_response_pcr_rd, y=custom.cis)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, color = 'gray') +
  theme_bw() + 
  xlab("pathologic_response_pcr_rd") + 
  ylab("CIS w5 top20 predictors") + 
  ggtitle("CIS and Pathologic Response") + theme(legend.position="bottom")

p2.cis <- ggplot(cp.dat, aes(x=custom.cis, y=oob.surv)) + 
  geom_point(color = 'gray') +
  stat_cor() +
  facet_wrap(~pathologic_response_pcr_rd) +
  xlab("CIS w5 top20 predictors") + 
  ylab('Relapse Risk') +
  ggtitle('CIS w5 top20 predictors (GSVA) vs Relapse Risk') +
  theme_bw()


p3.cis <- ggplot(cp.dat, aes(x=cis.hi, y=oob.surv)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, color = 'gray') +
  facet_wrap(~pathologic_response_pcr_rd) +
  theme_bw() + 
  xlab("CIS top20 predictive genes D10 to W5") + ylab("Relapse Risk") + 
  ggtitle("custom.CIS and Relapse Risk") + theme(legend.position="bottom")

p3.cis <- ggplot(cp.dat, aes(x=relapse.risk.hi, y=custom.cis)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, color = 'gray') +
  facet_wrap(~pathologic_response_pcr_rd) +
  theme_bw() + 
  xlab("Replase Risk Status") + ylab("CIS w5 top20 predictors") + 
  ggtitle("custom.CIS and Relapse Risk") + theme(legend.position="bottom")

ggarrange(
  ggarrange(
    p1.cis, p3.cis, labels = c("A", "B"), widths = c(0.6, 1),
    common.legend = TRUE, legend = "bottom"), 
  p2.cis, labels = c("", "C"), ncol = 1)

wilcox.test(cp.dat$custom.cis ~ cp.dat$pathologic_response_pcr_rd)

cp.dat.RD <- cp.dat[cp.dat$pathologic_response_pcr_rd == 'RD', ]
wilcox.test(cp.dat.RD$custom.cis ~ cp.dat.RD$relapse.risk.hi)

ggplot(cp.dat, aes(x=cis.hi, y=oob.surv)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, color = 'gray') +
  facet_wrap(~pathologic_response_pcr_rd) +
  theme_bw() +
  xlab("CIS top20 predictive genes D10 to W5") + ylab("Relapse Risk") + 
  ggtitle("custom CIS and Relapse Risk") + theme(legend.position="bottom")
cp.dat.pCR <- cp.dat[cp.dat$pathologic_response_pcr_rd == 'RD', ]
cp.dat.pCR.cis.hi <- cp.dat.pCR[cp.dat.pCR$cis.hi == 'High', ]
cp.dat.pCR.cis.lo <- cp.dat.pCR[cp.dat.pCR$cis.hi == 'Low', ]
wilcox.test(cp.dat.pCR.cis.hi$oob.surv, cp.dat.pCR.cis.lo$oob.surv)

# cis acute
p1.cis.acute <- ggplot(cp.dat, aes(x=pathologic_response_pcr_rd, y=custom.cis.acute)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, color = 'gray') +
  theme_bw() + 
  xlab("pathologic_response_pcr_rd") + 
  ylab("CIS d10 top20 predictors") + 
  ggtitle("CIS and Pathologic Response") + theme(legend.position="bottom")

p2.cis.acute <- ggplot(cp.dat, aes(x=custom.cis.acute, y=oob.surv)) + 
  geom_point(color = 'gray') +
  stat_cor() +
  facet_wrap(~pathologic_response_pcr_rd) +
  xlab("CIS d10 top20 predictors") + 
  ylab('Relapse Risk') +
  ggtitle('CIS d10 top20 predictors (GSVA) vs Relapse Risk') +
  theme_bw()

p3.cis.acute <- ggplot(cp.dat, aes(x=cis.acute.hi, y=oob.surv)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, color = 'gray') +
  facet_wrap(~pathologic_response_pcr_rd) +
  theme_bw() + 
  xlab("CIS top20 predictive genes D0 to D10") + ylab("Relapse Risk") + 
  ggtitle("custom.CIS.acute and Relapse Risk") + theme(legend.position="bottom")

p3.cis.acute <- ggplot(cp.dat, aes(x=relapse.risk.hi, y=custom.cis.acute)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, color = 'gray') +
  facet_wrap(~pathologic_response_pcr_rd) +
  theme_bw() + 
  xlab("Replase Risk Status") + ylab("CIS D10 top20 predictors") + 
  ggtitle("custom.CIS and Relapse Risk") + theme(legend.position="bottom")

ggarrange(
  ggarrange(
    p1.cis.acute, p3.cis.acute, labels = c("A", "B"), widths = c(0.5, 1),
    common.legend = TRUE, legend = "bottom"), 
  p2.cis.acute, labels = c("", "C"), ncol = 1)




p1.ISGmem <- ggplot(cp.dat, aes(x=pathologic_response_pcr_rd, y=ISGmem)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, color = 'gray') +
  theme_bw() + 
  xlab("pathologic_response_pcr_rd") + 
  ylab("ISG.mem") + 
  ggtitle("ISG.mem and Pathologic Response") + theme(legend.position="bottom")

p2.ISGmem <- ggplot(cp.dat, aes(x=ISGmem, y=oob.surv)) + 
  geom_point(color = 'gray') +
  stat_cor() +
  facet_wrap(~pathologic_response_pcr_rd) +
  ylab('Relapse Risk') +
  ggtitle('ISG.mem (GSVA) vs Relapse Risk') +
  theme_bw()

p3.ISGmem <- ggplot(cp.dat, aes(x=relapse.risk.hi, y=ISGmem)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, color = 'gray') +
  facet_wrap(~pathologic_response_pcr_rd) +
  theme_bw() + 
  xlab("Relapse Risk status") + ylab("ISG.mem") + 
  ggtitle("ISG.mem and Relapse Risk") + 
  theme(legend.position="bottom")

ggarrange(
  ggarrange(
    p1.ISGmem, p3.ISGmem, labels = c("A", "B"), widths = c(0.5, 1),
    common.legend = TRUE, legend = "bottom"), 
  p2.ISGmem, labels = c("", "C"), ncol = 1)




ggplot(cp.dat, aes(x = relapse.risk.hi, y = `Proteasomal.degradation`)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, color = 'gray') +
  facet_wrap(~pathologic_response_pcr_rd) +
  theme_bw()

cp.dat.melt <- cp.dat[, c('Sample', 'oob.surv', 'relapse.risk.hi', 'pathologic_response_pcr_rd', names(gene_list))]
cp.dat.melt <- melt(cp.dat.melt, id.vars = c('Sample','oob.surv','relapse.risk.hi', 'pathologic_response_pcr_rd'))

cp.dat.melt.pCR <- cp.dat.melt %>% 
  group_by(variable, pathologic_response_pcr_rd) %>%
  summarise(mean.PR = mean(value)) %>% 
  drop_na()
cp.dat.melt.pCR <- spread(cp.dat.melt.pCR, key = pathologic_response_pcr_rd, value = mean.PR)
cp.dat.melt.pCR$diff.PR <- cp.dat.melt.pCR$pCR - cp.dat.melt.pCR$RD

cp.dat.melt.relapseRisk <- cp.dat.melt %>%
  group_by(variable, relapse.risk.hi) %>%
  summarise(mean.RR = mean(value))
cp.dat.melt.relapseRisk <- spread(cp.dat.melt.relapseRisk, key = relapse.risk.hi, value = mean.RR)
cp.dat.melt.relapseRisk$diff.RR <- cp.dat.melt.relapseRisk$Low - cp.dat.melt.relapseRisk$High

cp.dat.comp <- merge(cp.dat.melt.relapseRisk, cp.dat.melt.pCR, by = 'variable')
var <- c('Cell.Cycle...G1.S', 'Cell.Cycle...G2.M', 'Cell.Cycle.HMG.rich', 'Chromatin', 'custom.cis', 'custom.cis.acute',
         'custom.cocl2', 'custom.dabtram', 'EMT.I', 'EMT.II', 'EMT.III', 'EMT.IV', 'Epithelial.Senescence',
         'Hypoxia', 'Interferon.MHC.II..I.', 'Interferon.MHC.II..II.', 'ISG.RS', 'MYC', 'Proteasomal.degradation',
         'Protein.maturation', 'Secreted.I', 'Secreted.II', 'Skin.pigmentation', 'Stress', 'Stress..in.vitro.',
         'Translation.initiation', 'Unfolded.protein.response')
cp.dat.comp <- cp.dat.comp[cp.dat.comp$variable %in% var,]
ggplot(cp.dat.comp, aes(x = diff.PR, y = diff.RR)) + 
  geom_point() + 
  # geom_smooth(method = 'lm') +
  ggrepel::geom_text_repel(aes(label = variable), box.padding = 0.5) +
  stat_cor() +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  labs(x = 'Higher in pCR', y = 'Higher in low relapse Risk') +
  theme_bw()
