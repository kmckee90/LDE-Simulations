# ---------------------------------------------------------------------
#   LDE Simulations
#   Basic Univ/Multivariate LDEs, no events
#   Author: Kevin McKee
#   Last updated: 9/14/2017
# ---------------------------------------------------------------------

library(deSolve)
library(OpenMx)
library(psych)
library(MASS)

source("lib/GLLAfunctions.R")
source("lib/fitODE.R")
source("lib/EventAnalysis.R")
source("lib/plotSim.R")
source("lib/klmFuncs.R")
source("lib/MS/MSD_simPar.R")
source("lib/MS/MSD_genData.R")

mxOption(NULL, 'Default optimizer', 'NPSOL')
options(width=80)
options(scipen=999)
genPlotting=T

# Settings ----------------------------------------------------------------
genPlotting=T
par(mfrow=c(1,1), mai=c(0.4,1,0.4,1))


defaults<-list(    
  "N"=4*30,
  "nVars"= 1,
  "nTrials"=5,
  "SNR"=8,
  "embedD1"=4,
  "embedD2"=4,
  "theTau" = 1,
  "deltaT" = 1,
  "eventScale"=1,
  "tryHard.tries"=10
  
)

# defaultEnv<-list2env(defaults)
# attach(defaultEnv)


# Generate parameters and initial conditions ------------------------------
source("lib/MS/MSD_matrices.R")
source("lib/MS/MS_model.R")
source("lib/MS/MS_constraints.R")

trueVals<-initCond<-NULL
for(i in 1:defaults$nTrials){
  randIntercept<-runif(defaults$nVars, -1, 1)
  randSlope<-runif(defaults$nVars, -0.05, 0.05)
  randMean<-cbind(randIntercept, randSlope)
  
  randEta<-runif(defaults$nVars, -.7, -0.1)
  randZeta<-runif(defaults$nVars, -0.5, 0)
  randIndDyn<-cbind(randEta, randZeta)
  
  randEtaM<-runif(defaults$nVars, -.09, -0.01)
  randZetaM<-runif(defaults$nVars, -0.05, -0.005)
  randIndDynM<-cbind(randEtaM, randZetaM)

  randGamma.L<-randGamma.S<-diag(1, defaults$nVars)
  randGammas<-which(randGamma.L==0)
  randGamma.L[randGammas]<-runif(length(randGammas), -0.6, 0.2)
  randGamma.S[randGammas]<-runif(length(randGammas), -0.6, 0.2)
  
  randPars<-c(as.vector(randIndDyn), as.vector(randIndDynM), as.vector(randGamma.L), as.vector(randGamma.S), as.vector(randMean))
  names(randPars)<-c("etaX1", "zetaX1", "mEtaX1", "mZetaX1", "gammaX11L", "gammaX11S", as.vector(meanLabs))
  trueVals<-rbind(trueVals, randPars)
  initCond<-rbind(initCond, as.vector( rbind( runif(defaults$nVars*2, 0.5, 2) * (rbinom(defaults$nVars*2, 1, 0.5)*2-1), 0)) )
}


dat<-MSD_genData(N=100, SNR=8, trueVals[1,], initCond[1,], defaults,
            eventScale=1, plotting=T)


tEmbedded<- gllaEmbed(dat, embed=defaults$embedD1*defaults$embedD2, tau=defaults$theTau, idColumn=FALSE)


tEmbedded<-NULL
disConts<-seq(1, defaults$N, defaults$embedD1)
for(i in 1:defaults$nVars){
  dat.emb <- gllaEmbed(dat[,i], embed=defaults$embedD, tau=defaults$theTau, idColumn=FALSE)
  dat.emb<-dat.emb[disConts[disConts<nrow(dat.emb)],]
  tEmbedded <- cbind(tEmbedded, dat.emb)
}
colnames(tEmbedded)<-manifestVars
tEmbedded<-cbind(tEmbedded, "Occasion"=1:nrow(tEmbedded))

model$data<-mxData(tEmbedded, type="raw")
model<-mxTryHard(model, extraTries = defaults$tryHard.tries)

simMSD<-MSD_simPar("pEvent", 1, trueVals=trueVals, useEventAnalysis=F, defaults=defaults)

#To do:
#N/day, embedD, second level discontinuity


