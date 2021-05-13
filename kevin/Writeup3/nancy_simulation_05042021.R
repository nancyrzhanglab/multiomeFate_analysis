rm(list=ls())
set.seed(10)

library(fields) # for its image.plot function.
library(umap) 



### Parameters governing randomness in the ATAC data
generate_nancy_data <- function(){
  
  # the "active" time interval for each phase.  
  makeActiveIntervals<-function(timesteps, nPhases, phaseGapFraction=0){
    activeInterval=matrix(nrow=nPhases, ncol=2, data=0)
    phaseLength=timesteps/nPhases
    gap=round(phaseLength*phaseGapFraction)
    for(i in 1:nPhases){
      if(i==1){
        activeInterval[i,1]=1
      } else {
        activeInterval[i,1]=phaseLength*(i-1)+gap
      }
      if(i==nPhases){
        activeInterval[i,2]=timesteps
      } else {
        activeInterval[i,2]=phaseLength*i-gap
      }
    }
    activeInterval
  }
  
  # in which phase is a gene active?
  # returns vector of length nGenes, with values belonging to 1,...,nPhases
  makeGenePhases<-function(nGenes,nPhases){
    genePhase=rep(0,0)
    genesPerPhase=nGenes/nPhases
    for(i in 1:nPhases){
      genePhase=c(genePhase, rep(i,genesPerPhase))
    }
    genePhase
  }
  
  # in which branch is a gene active?
  # returns vector of length nGenes, with values "0"= common gene,
  # "1" branch 1 gene, "2" branch 2 gene, ...
  makeGeneBranches<-function(nGenes,nBranches){
    geneBranch=rep(0,0)
    for(i in 1:nPhases){
      geneBranch = c(geneBranch, rep(0,nCommonPhaseGenes[i]))
      for(k in 1:nBranches){
        geneBranch=c(geneBranch, rep(k,nBranchInformativeGenes[k,i]) )
      }
    }
    geneBranch
  }
  
  # For each gene and each peak in enhancer set for gene, how long is the lead time
  # for the peak?
  makeChromatinLead<-function(nGenes, nRegionsPerGene, timesteps,
                              sampleMethod="step",chromatinLeadStep=0.01, maxChromatinLead=0.1){
    if(sampleMethod=="step"){
      chromatinLead=matrix(nrow=nGenes, ncol=nRegionsPerGene, data=0)
      for(k in 1:nRegionsPerGene){
        chromatinLead[,k]=k*round(chromatinLeadStep*timesteps)
      }
    }
    if(sampleMethod=="random"){
      chromatinLead=matrix(nrow=nGenes, ncol=nRegionsPerGene, 
                           data=sample(round(maxChromatinLead*timesteps),nGenes*nRegionsPerGene))
    }
    chromatinLead
  }
  
  # For each gene and each peak in enhancer set for gene, how long is the lag time
  # for the peak?
  makeChromatinLag<-function(nGenes, nRegionsPerGene, timesteps, sampleMethod="step",
                             chromatinLagStep=0.005, maxChromatinLag=0.01){
    if(sampleMethod=="step"){
      chromatinLag=matrix(nrow=nGenes, ncol=nRegionsPerGene, data=0)
      for(k in 1:nRegionsPerGene){
        chromatinLag[,k]=k*round(chromatinLagStep*timesteps)
      }
    }
    if(sampleMethod=="random"){
      chromatinLag=matrix(nrow=nGenes, ncol=nRegionsPerGene, 
                          data=sample(round(maxChromatinLag*timesteps),nGenes*nRegionsPerGene))
    }
    chromatinLag
  }
  
  ## make chromatinOpen matrices.  These can be thought of as the "true" open/close
  ## states from which the ATAC data "x" is sampled.
  ## chromatinOpen is a list of length nBranches.  
  ## Cells in branch b follow the chromatin profile of chromatinOpen[[b]].
  ## chromatinOpen[[b]][g,i,t]=1 if region i of gene g is open at time t for a cell belonging to branch b. 
  makeChromatinOpen<-function(nGenes, nRegionsPerGene, timesteps,geneBranch, genePhase, 
                              chromatinLead, chromatinLag, noise.method=1,
                              noisy.gene.frac=0.05, noisy.gene.duration=0.01,
                              noisy.lead=1, noisy.lag=1){
    chromatinOpen=vector("list",nBranches)
    for(b in 1:nBranches){
      chromatinOpen[[b]]=array(dim=c(nGenes,nRegionsPerGene,timesteps), data=0)  
      for(i in 1:nGenes){
        if(genePhase[i]==0){
          if(noise.method==1){
            # the noisy genes are not smooth in time.
            noisy.num.periods = ceiling(noisy.gene.frac*timesteps)
            noisy.duration=1
            geneStart = sample(timesteps, noisy.num.periods)
            geneEnd=geneStart+rpois(noisy.num.periods, noisy.duration)
            
          } else {
            # the noisy genes are smooth in time and have chromatin lead and lag.
            noisy.num.periods = ceiling(noisy.gene.frac*timesteps)
            noisy.duration = ceiling(noisy.gene.duration*timesteps)
            geneStart = sample(timesteps, noisy.num.periods)
            geneEnd=geneStart+rpois(noisy.num.periods, noisy.duration)
          }
        }
        if(geneBranch[i]==b || geneBranch[i]==0){ # if gene is active in this branch, or if gene is active in all branches
          for(k in 1:nRegionsPerGene){
            if(genePhase[i]!=0){
              # Pseudotime and branch informative gene
              startOpen=max(1,activeInterval[genePhase[i],1]-chromatinLead[i,k])
              endOpen=min(timesteps,activeInterval[genePhase[i],2]+chromatinLag[i,k])
              chromatinOpen[[b]][i,k,startOpen:endOpen]=1
            } else {
              # Pure noise gene
              if(noise.method==1){
                # the noisy genes are not smooth in time, also decrease chromatinLag
                chromatinLead[i,k] = rpois(noisy.lead,1)
                chromatinLag[i,k] = rpois(noisy.lag,1)
              } 
              startOpen=pmax(1, geneStart-chromatinLead[i,k])
              endOpen = pmin(timesteps,geneEnd+chromatinLag[i,k])
              for(j in 1:noisy.num.periods){
                chromatinOpen[[b]][i,k,startOpen[j]:endOpen[j]]=1
              }
            }
          }  
        }
      }
    }
    chromatinOpen  
  }
  
  # Sample the ATAC for *one* cell given its gene-by-region chromatinOpen matrix.
  samplePeaks<-function(chromatinOpen, closed.lambda=0.1, open.lambda=2){
    x=array(dim=dim(chromatinOpen))
    inds= which(chromatinOpen==0, arr.ind=TRUE)
    x[inds] = rpois(nrow(inds), lambda=closed.lambda)
    inds=which(chromatinOpen>0, arr.ind=TRUE)
    x[inds]=rpois(nrow(inds), lambda=open.lambda)
    x
  }
  
  # Compute the chromatin potential given a gene-by-region matrix of peak counts for one cell.
  computePotential<-function(peaks){
    potential=apply(peaks,1,sum)  
    potential
  }
  
  
  # Sample the RNA expression for all genes in one cell (vector of length nGenes)
  # given the chromatin potential.
  sampleRNA<-function(potential, potential.thresh, betaG=1, sigmaG=0.5){
    potential2 = pmax(0,potential-potential.thresh)
    y=rnorm(length(potential), mean=betaG*potential2, sd=sigmaG)
    y  
  } 
  
  # flatten the first two dimensions of a 3-d array
  flatten<-function(arr){ 
    newarr = matrix(nrow = prod(dim(arr)[1:2]), ncol=dim(X)[3])
    for(i in 1:dim(arr)[3]){
      newarr[,i]=as.vector(t(arr[,,i]))
    }
    newarr
  }
  
  # gets any number of atomic variables and put them into a string.  
  # useful for debugging.
  # try makeVariableValueString(sigmaG, betaG)
  makeVariableValueString<-function(...){
    varnames=lapply(substitute(list(...))[-1], deparse)
    li=list(...)
    names(li)=varnames
    str=""
    for(i in 1:length(li)){
      str=paste(str,varnames[i], "=",paste(li[[i]]), sep="")
      if(i<length(li)) str=paste(str,", ", sep="")
    }
    str
  }
  
  # ---------------------------------------------
  # Set up model
  # ---------------------------------------------
  
  ### How many branches?  How many phases?  Active genes per branch per phase
  nPhases = 20   # Number of phases of gene activation.
  nBranches = 2
  nGenesPerBranchPerPhase = 5
  nGenesPerPhase = nGenesPerBranchPerPhase*nBranches  # Number of genes activated per phase, need to be divisible by nBranches
  # number of branch informative genes per phase, per branch
  nBranchInformativeGenes=c(rep(0,floor((nPhases-nGenesPerBranchPerPhase)/2)), 
                            c(1:nGenesPerBranchPerPhase), 
                            rep(nGenesPerBranchPerPhase,ceiling((nPhases-nGenesPerBranchPerPhase)/2)))  
  nBranchInformativeGenes=matrix(nrow=nBranches, ncol=nPhases, data=rep(nBranchInformativeGenes,nBranches), byrow=TRUE)
  rownames(nBranchInformativeGenes)=paste("branch", c(1:nBranches), sep="")
  # number of branch uninformative genes per phase
  nCommonPhaseGenes = rep(nGenesPerPhase, nPhases)-colSums(nBranchInformativeGenes)
  nGenes=nGenesPerPhase*nPhases
  ### Pure noise genes
  nNoisyGenes = nGenes
  nTotalGenes = nGenes+nNoisyGenes
  
  ### How many cells?  What are their true branches, times?
  timesteps = 500  # time t=1,2,...
  ncells=timesteps  # need to be divisible by number of branches
  trueTime=c(1:timesteps) # true time of cells
  trueBranch = rep(c(1,2), ncells/2) # true branch that each cell belongs to.
  
  ### Parameters governing ATAC
  nRegionsPerGene=5 # number of ATAC peak regions per gene
  phaseGapFraction=0.05
  chromatinLeadStep=0.1
  chromatinLagStep=0
  
  ### Set the active interval for each phase.
  activeInterval = makeActiveIntervals(timesteps=timesteps, nPhases=nPhases, phaseGapFraction=phaseGapFraction)
  ### Set the active phase and active branch for each gene.
  # noisy genes have phase "0".
  genePhase = c(makeGenePhases(nGenes, nPhases), rep(0, nNoisyGenes))
  
  # In which branch is a gene active?  "0"=both, "1"=1, "2"=2
  # noisy genes don't distinguish between branch.
  geneBranch = c(makeGeneBranches(nGenes, nBranches), rep(0, nNoisyGenes))
  
  ## Make the chromatin lead, a.k.a. how long before a gene's activation phase
  ## does a peak open
  chromatinLead = makeChromatinLead(nTotalGenes,nRegionsPerGene, timesteps, 
                                    sampleMethod="step", chromatinLeadStep=chromatinLeadStep)
  
  ## Make the chromatin lag, a.k.a. how long after a gene's activation phase
  ## does a peak close?
  chromatinLag = makeChromatinLag(nTotalGenes,nRegionsPerGene, timesteps, 
                                  sampleMethod="step", chromatinLagStep=chromatinLagStep)
  
  ## make chromatinOpen matrices.  These can be thought of as the "true" open/close
  ## states from which the ATAC data "x" is sampled.
  ## chromatinOpen is a list of length nBranches.  
  ## Cells in branch b follow the chromatin profile of chromatinOpen[[b]].
  ## chromatinOpen[[b]][g,i,t]=1 if region i of gene g is open at time t for a cell belonging to branch b. 
  chromatinOpen=makeChromatinOpen(nGenes=nTotalGenes, nRegionsPerGene=nRegionsPerGene, timesteps=timesteps, 
                                  geneBranch=geneBranch, genePhase=genePhase, chromatinLead=chromatinLead, chromatinLag=chromatinLag,
                                  noisy.gene.frac = 0.01, noisy.gene.duration=0.002)
  
  # -----------------------------------
  # Sample X=ATAC, Y=RNA for each cell.
  # -----------------------------------
  
  open.lambda = 2
  closed.lambda= 0.5
  ### Parameters governing the ATAC-RNA relationship
  sigmaG = 1
  betaG = 0.5
  potentialThreshold=0.6
  
  X=array(dim=c(nTotalGenes, nRegionsPerGene, ncells))
  Y=array(dim=c(nTotalGenes, ncells))
  chromPotential=array(dim=c(nTotalGenes,ncells))
  
  potential.thresh=open.lambda*nRegionsPerGene*potentialThreshold
  
  set.seed(10)
  for(i in 1:ncells){
    # get corresponding peak profile according to true time and branch
    X[,,i] = samplePeaks(chromatinOpen[[trueBranch[i]]][,,trueTime[i]],
                         closed.lambda=closed.lambda, open.lambda=open.lambda)  
    chromPotential[,i] = computePotential(X[,,i])
    
    Y[,i] = sampleRNA(chromPotential[,i],potential.thresh, betaG=betaG, sigmaG=sigmaG)
  }
  Xflat=flatten(X)
  
  list(X = Xflat, Y = Y, 
       trueBranch = trueBranch,
       trueTime = trueTime,
       chromatinOpen = chromatinOpen, 
       chromPotential = chromPotential,
       potential.thresh = potential.thresh, nRegionsPerGene = nRegionsPerGene)
}

# -------------------------------
# Visualize X, Y
# -------------------------------

# set.seed(10)
# X.umap=umap(t(Xflat), n_neighbors=30)
# Y.umap=umap(t(Y), n_neighbors=30)
# pal = colorRampPalette(c("blue", "red"))
# 
# atacTitleString=paste("ATAC:",makeVariableValueString(open.lambda,closed.lambda,chromatinLeadStep))
# rnaTitleString=paste("RNA:",makeVariableValueString(potentialThreshold, betaG, sigmaG))
# 
# png("../../out/fig/writeup3/simdata.png", height=800, width=1200)
# par(mfrow=c(2,3),oma = c(0, 0, 2, 0))
# image.plot(1:timesteps, 1:nrow(Xflat), t(Xflat), xlab="True time", ylab="peak", main="ATAC", cex.main=1.5)
# plot(X.umap$layout, main="ATAC, colored by true branch", xlab="Umap 1", ylab="Umap 2",col=trueBranch, cex.main=1.5)
# plot(X.umap$layout, main="ATAC, colored by true time",xlab="Umap 1", ylab="Umap 2", col=pal(timesteps)[trueTime], cex.main=1.5)
# image.plot(1:timesteps, 1:nTotalGenes, t(Y), xlab="True time", ylab="Gene", main="RNA", cex.main=1.5)
# plot(Y.umap$layout, main="RNA, colored by true branch",xlab="Umap 1", ylab="Umap 2", col=trueBranch,cex.main=1.5)
# plot(Y.umap$layout, main="RNA, colored by true time", xlab="Umap 1", ylab="Umap 2",col=pal(timesteps)[trueTime], cex.main=1.5)
# mtext(paste(atacTitleString, rnaTitleString, sep="   "), outer = TRUE, cex = 1)
# 
# dev.off()
# 
# 
