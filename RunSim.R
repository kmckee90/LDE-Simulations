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
source("lib/NL/N_simPar.R")
source("lib/NL/N_genData.R")

mxOption(NULL, 'Default optimizer', 'NPSOL')
options(width=80)
options(scipen=999)
genPlotting=T

# Settings ----------------------------------------------------------------
genPlotting=T
par(mfrow=c(1,1), mai=c(0.4,1,0.4,1))

defaults<-list(    
  "N"=45,
  "nVars"= 1,
  "nTrials"=3,
  "SNR"=4,
  "pEvent"=0.0,
  "eventNormal"=F,
  "embedD"=4,
  "theTau" = 1,
  "deltaT" = 1
)



defaultEnv<-list2env(defaults)
attach(defaultEnv)


# Generate parameters and initial conditions ------------------------------
  source("lib/NL/N_matrices.R")
  source("lib/NL/N_model.R")
  source("lib/NL/N_constraints.R")
    
  trueVals<-initCond<-NULL
  for(i in 1:defaults$nTrials){
    randIntercept<-runif(defaults$nVars, -1, 1)
    randSlope<-runif(defaults$nVars, -0.05, 0.05)
    randMean<-cbind(randIntercept, randSlope)
    
    randEta<-runif(defaults$nVars, -.7, -0.1)
    randZeta<-runif(defaults$nVars, -0.5, 0)
    randIndDyn<-cbind(randEta, randZeta)
    
    randGamma.L<-randGamma.S<-diag(1, defaults$nVars)
    randGammas<-which(randGamma.L==0)
    randGamma.L[randGammas]<-runif(length(randGammas), -0.6, 0.2)
    randGamma.S[randGammas]<-runif(length(randGammas), -0.6, 0.2)
    
    randPars<-c(as.vector(randIndDyn), as.vector(randGamma.L), as.vector(randGamma.S), as.vector(randMean))
    names(randPars)<-c(as.vector(p.l), as.vector(gL.l), as.vector(gS.l), as.vector(meanLabs))
    trueVals<-rbind(trueVals, randPars)
    initCond<-rbind(initCond, as.vector( rbind( runif(defaults$nVars, 0.5, 2) * (rbinom(defaults$nVars, 1, 0.5)*2-1), 0)) )
  }

# rm(defaultEnv)
# Run simulation ----------------------------------------------------------
  par(mfrow=c(1,1), mai=c(1,1,1,1))
  
  
sim.pEvent<-N_simPar("pEvent", 0.1, trueVals=trueVals, useEventAnalysis=T, defaults=defaults)



sim.pEvent.lm<-lm(sim.pEvent$Estimates[,"etaX1"]~sim.pEvent$trueVals[,"etaX1"])
  
  
  
  
  

# Generate data -----------------------------------------------------------
dat<-N_genData(45, SNR=10, numVars=nVars, trueVals[2,], initCond[2,], 
               events=0.1, eventNormal=F, eventType="level", eventScale = 1,
               plotting=T)


tEmbedded<-NULL
for(i in 1:nVars){
  tEmbedded <- cbind(tEmbedded, gllaEmbed(dat[,i], embed=embedD, tau=theTau, idColumn=FALSE))
}
colnames(tEmbedded)<-manifestVars
tEmbedded<-cbind(tEmbedded, "Occasion"=1:nrow(tEmbedded))

ea<-eventAnalysis(tEmbedded, model, min.data = 10, tries=3)
ea<-eventAnalysis(tEmbedded, model, min.data = 10, tries=12, use=ea$info)



source("lib/MS/MS_defaults.R")
source("lib/MS/MS_genData.old.R")
model$data<-mxData(tEmbedded, type="raw")
model<-mxRun(model)
model<-mxTryHard(model, extraTries = 100, greenOK=T)
model<-mxTryHardWideSearch(model, extraTries = 20, greenOK=T)
summary(model)














defaults<-list(    
  "N"=3*50,
  "nVars"= 1,
  "nTrials"=5,
  "SNR"=4,
  "pEvent"=1,
  "eventNormal"=F,
  "embedD1"=3,
  "embedD2"=4,
  "theTau" = 1,
  "deltaT" = 1,
  "eventScale"=1,
  "eventScaleRatio"=1/3,
  "eventTypeX1"="level",
  "eventTypeX2"="level"
  
)

defaultEnv<-list2env(defaults)
attach(defaultEnv)


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
  names(randPars)<-c(c("etaX1", "zetaX1", "etaXM1", "zetaXM1"), as.vector(gL.l), as.vector(gS.l), as.vector(meanLabs))
  trueVals<-rbind(trueVals, randPars)
  initCond<-rbind(initCond, as.vector( rbind( runif(defaults$nVars*2, 0.5, 2) * (rbinom(defaults$nVars*2, 1, 0.5)*2-1), 0)) )
}

dat<-MSD_genData(150, SNR=10, trueVals[3,], initCond[3,], defaults,
                events=0.1, eventNormal=F, eventType="level", eventScale = 1,
                plotting=T)


simMSD<-MSD_simPar("pEvent", 0.5, trueVals=trueVals, useEventAnalysis=F, defaults=defaults)




