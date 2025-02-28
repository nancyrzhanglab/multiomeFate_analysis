rm(list = ls())

library(Biobase)
library(survival)
library(randomForestSRC)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data_other/Hatzis_et_al_breast_cancer_clinical_microarray_GSE25065/'

irds <- c("FLNB", "GLUD1", "IFI44", "IFIT1", "IFIT3", "IGF1R", "ISG15", 
          "MX1", "OAS1", "STAT1", "UBE2D3", "ZNF273", "CTDSPL", "ABCD3")
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
# Process data
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
rsf.dat <- data.frame(cen=cen, rfs=rfs, tsp.irds=tsp.irds, pData(eset)[, use.index], RS=RS)
rsf.dat[, "grade"] <- as.factor(ifelse(rsf.dat$grade=="4=Indeterminate", "Indeterminate", rsf.dat$grade))

# ==============================================================================
# KM
# ==============================================================================
temp.dat <- data.frame(rsf.dat, dlda30=ifelse(pData(MDACCrmaSet)$dlda30_prediction=="pCR", "Pred pCR", "Pred RD"))

# pCR vs RFS
sfit.pcr <- survfit(Surv(rfs, cen) ~ pathologic_response_pcr_rd, data=temp.dat)
sdiff.pcr <- survdiff(Surv(rfs, cen) ~ pathologic_response_pcr_rd, data=temp.dat, na.action="na.omit")

plot(sfit.pcr, main="pCR and RFS", col=c("blue", "brown"), xlab="Years", ylab="Proportion Relapse-Free")
legend("bottomleft", legend=c("pCR", "RD"), lwd=2, col=c("blue", "brown"))
text(3, 0.30, paste("p=", round(1 - pchisq(sdiff.pcr$chisq, 1), 4), sep=""), cex=1.5)

# ISG.RS vs RFS
sfit.irds <- survfit(Surv(rfs, cen) ~ tsp.irds < 2, data=temp.dat)
sdiff.irds <- survdiff(Surv(rfs, cen) ~ tsp.irds < 2, data=temp.dat, na.action="na.omit")

plot(sfit.irds, main="IRDS and RFS", col=c("brown", "blue"), xlab="Years", ylab="Proportion Relapse-Free")
legend("bottomleft", legend=c("IRDS(+)", "IRDS(-)"), lwd=2, col=c("brown", "blue"))
text(3, 0.30, paste("p=", round(1 - pchisq(sdiff.irds$chisq, 1), 4), sep=""), cex=1.5)

# ==============================================================================
# Calculate relapse risk and plot boxplots
# ==============================================================================

ntree <- 1000
nodesize <- 1
nsplit <- 2

ntree.null    <- c(10, 100, 1000, 10000)[3]

n.centroids <- c(4, 5)[1]  # using 4 will exclude "normal" subtype similar to Parker et al
conservative <- FALSE  # used in max.subtree and var.select functions
binIRDS <- TRUE  # used in coplots to show IRDS as binary or continuous variable
coplot.min <- c(1, 5, 10)[3]  # used in coplots to avoid plotting when only a few data points
mcrep <- 20  # for Monte Carlo reps used in var.select

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


## Show coplots
if (binIRDS) temp.irds <-  as.factor(rfit[["xvar"]][, 'tsp.irds'] >= 2)
if (!binIRDS) temp.irds <- rfit$predictors$tsp.irds
cp.dat <- data.frame(oob.surv=rfit$predicted.oob/rfit$n, irds=temp.irds, rfit$xvar)
cp.dat <- data.frame(oob.surv=rfit$predicted.oob/rfit$n, rfit$xvar)

form <- vector("list", 1)
form[[1]] <- as.formula("oob.surv ~ irds | pathologic_response_pcr_rd")

coplot(form[[1]], number=5, overlap=0, axlabels = function(x) NULL, data=cp.dat,
       subscripts = TRUE, xlab="IRDS Status", ylab="Relapse Risk",
       panel = function(x, y, subscripts, ...) {
         err <- get("vimp", "package:randomForestSRC")(rfit, subset=which(subscripts))$err.rate
         plotit <- ifelse(length(y) < coplot.min, FALSE, TRUE)
         boxplot(y ~ x, add=TRUE, boxwex=0.25, at=c(1.25, 1.75), yaxt="n", outline=FALSE, 
                 col=c("royal blue", "red"), names=c("IRDS-", "IRDS+"), plot=plotit)
         if (plotit) legend(1, 0.5, legend=paste("Err:", round(err, 3)), cex=0.8, 
                            bg=rgb(0.5, 0.5, 0, alpha=0.3), xjust=0)
       })



cp.dat$Sample <- rownames(cp.dat)
cp.dat <- merge(cp.dat, gsva.es.df, by="Sample", all.x=TRUE)
cp.dat$ISG.RS.hi <- ifelse(cp.dat$ISG.RS > median(cp.dat$ISG.RS, na.rm=TRUE), "High", "Low")
cp.dat$Stress.hi <- ifelse(cp.dat$Stress > median(cp.dat$Stress, na.rm=TRUE), "High", "Low")
cp.dat$cis.hi <- ifelse(cp.dat$custom.cis > median(cp.dat$custom.cis, na.rm=TRUE), "High", "Low")

ggplot(cp.dat, aes(x=ISG.RS, y=oob.surv)) + 
  geom_point() + 
  geom_smooth(method="lm", se=FALSE, color="red") +
  stat_cor(method = 'spearman') +
  facet_wrap(~pathologic_response_pcr_rd.x) +
  theme_bw() + xlab("IRDS Status") + ylab("Relapse Risk") + 
  ggtitle("IRDS and Relapse Risk") + theme(legend.position="bottom")
ggplot(cp.dat, aes(x=custom.cis, y=oob.surv)) + 
  geom_point() + 
  geom_smooth(method="lm", se=FALSE, color="red") +
  stat_cor(method = 'spearman') +
  facet_wrap(~pathologic_response_pcr_rd) +
  theme_bw() + xlab("IRDS Status") + ylab("Relapse Risk") + 
  ggtitle("IRDS and Relapse Risk") + theme(legend.position="bottom")

ggplot(cp.dat, aes(x=ISG.RS.hi, y=oob.surv)) + 
  geom_boxplot() + 
  facet_wrap(~pathologic_response_pcr_rd.x) +
  theme_bw() + xlab("IRDS Status") + ylab("Relapse Risk") + 
  ggtitle("IRDS and Relapse Risk") + theme(legend.position="bottom")

ggplot(cp.dat, aes(x=Stress.hi, y=oob.surv)) + 
  geom_boxplot() + 
  facet_wrap(~pathologic_response_pcr_rd.x) +
  theme_bw() + xlab("IRDS Status") + ylab("Relapse Risk") + 
  ggtitle("IRDS and Relapse Risk") + theme(legend.position="bottom")

ggplot(cp.dat, aes(x=cis.hi, y=oob.surv)) + 
  geom_boxplot() + 
  facet_wrap(~pathologic_response_pcr_rd) +
  theme_bw() + xlab("IRDS Status") + ylab("Relapse Risk") + 
  ggtitle("IRDS and Relapse Risk") + theme(legend.position="bottom")

cp.dat.pCR <- cp.dat[cp.dat$pathologic_response_pcr_rd == 'pCR', ]
cp.dat.pCR.ISG.hi <- cp.dat.pCR[cp.dat.pCR$ISG.RS.hi == 'High', ]
cp.dat.pCR.ISG.lo <- cp.dat.pCR[cp.dat.pCR$ISG.RS.hi == 'Low', ]

cp.dat.pCR.Stress.hi <- cp.dat.pCR[cp.dat.pCR$Stress.hi == 'High', ]
cp.dat.pCR.Stress.lo <- cp.dat.pCR[cp.dat.pCR$Stress.hi == 'Low', ]

cp.dat.pCR.cis.hi <- cp.dat.pCR[cp.dat.pCR$cis.hi == 'High', ]
cp.dat.pCR.cis.lo <- cp.dat.pCR[cp.dat.pCR$cis.hi == 'Low', ]

wilcox.test(cp.dat.pCR.ISG.hi$oob.surv, cp.dat.pCR.ISG.lo$oob.surv)
wilcox.test(cp.dat.pCR.Stress.hi$oob.surv, cp.dat.pCR.Stress.lo$oob.surv)
wilcox.test(cp.dat.pCR.cis.hi$oob.surv, cp.dat.pCR.cis.lo$oob.surv)
