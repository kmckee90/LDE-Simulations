# ---------------------------------------------------------------------
#   LDE Simulations
#   Basic Univ/Multivariate LDEs
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
par(mfrow=c(1,1), mai=c(1,1,1,1))

# Settings ----------------------------------------------------------------
embedD.vals<- 3:20
N.vals<- c(15, 30, 45, 60, 90, 120, 180, 240, 300)  #[1:4]
pEvent.vals<- round(seq(0, 0.3, 0.02),2)               #[1:6]
SNR.vals<- c(1/2, 1, 2, 4, 8, 16)        # [c(3,4,5,7,10)]
eventScale.vals<- seq(0.25, 4, 0.25)        # [c(3,4,5,7,10)]

# Defaults, parameters, and initial conditions ------------------------------
defaults<-list(    
  "nTrials"=50,
  
  "nVars"= 4,
  
  "N"=45,
  "SNR"=8,
  "pEvent"=0,
  "embedD"=4,
  "eventScale"=1,
  "eventNormal"=F,
  
  "theTau" = 1,
  "deltaT" = 1,
  "tryHard.tries" = 20,
  
  "plotting"=F
)

  source("lib/NL/N_matrices.R")
  source("lib/NL/N_model.R")
  source("lib/NL/N_constraints.R")
    
  trueVals<-initCond<-NULL
  for(i in 1:defaults$nTrials){
    randIntercept<-runif(defaults$nVars, -1, 1)
    randSlope<-runif(defaults$nVars, -0.05, 0.05)
    randMean<-cbind(randIntercept, randSlope)
    
    randEta<-runif(defaults$nVars, -1, -0.01)
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
  
  
  
# Univariate Random Effects ----------------------------------------------------------
defaults$nVars<-1
plotParams<-names(omxGetParameters(model))[grepl("eta",names(omxGetParameters(model)))][1:(defaults$nVars*2)]

sim.1.embedD<-N_simPar("embedD", embedD.vals, trueVals=trueVals, useEventAnalysis=F, defaults=defaults)
sim.1.N<-N_simPar("N", N.vals, trueVals=trueVals, useEventAnalysis=F, defaults=defaults)
sim.1.SNR<-N_simPar("SNR", SNR.vals, trueVals=trueVals, useEventAnalysis=F, defaults=defaults)
sim.1.pEvent<-N_simPar("pEvent", pEvent.vals, trueVals=trueVals, useEventAnalysis=F, defaults=defaults)

defaults$pEvent<-0.1
sim.1.eventScale<-N_simPar("eventScale", eventScale.vals, trueVals=trueVals, useEventAnalysis=F, defaults=defaults)
defaults$pEvent<-0.0


# Univariate Fixed Effects ------------------------------------------------

trueVals.fixed.X1<-trueVals
trueVals.fixed.X1[,"etaX1"]<- -0.5
trueVals.fixed.X1[,"zetaX1"]<- -0.3

sim.1.embedD<-N_simPar("embedD", embedD.vals, trueVals=trueVals, useEventAnalysis=F, defaults=defaults)
sim.1.N.fixed<-N_simPar("N", N.vals, trueVals=trueVals.fixed.X1, useEventAnalysis=F, defaults=defaults)
sim.1.SNR.fixed<-N_simPar("SNR", SNR.vals, trueVals=trueVals.fixed.X1, useEventAnalysis=F, defaults=defaults)
sim.1.pEvent.fixed<-N_simPar("pEvent", pEvent.vals, trueVals=trueVals.fixed.X1, useEventAnalysis=F, defaults=defaults)

defaults$pEvent<-0.05
sim.1.eventScale.fixed<-N_simPar("eventScale", N.vals, trueVals=trueVals.fixed.X1, useEventAnalysis=F, defaults=defaults)
defaults$pEvent<-0.0


# 2 Vars, fixed individual parameters ------------------------------------------------------------------
defaults$nVars<-2
plotParams<-names(omxGetParameters(model))[grepl("eta",names(omxGetParameters(model)))][1:(defaults$nVars*2)]

trueVals.fixed.X1.X2<-trueVals.fixed.X1
trueVals.fixed.X1.X2[,"etaX2"]<- -0.25
trueVals.fixed.X1.X2[,"zetaX2"]<- -0.15

sim.2.embedD<-N_simPar("embedD", embedD.vals, trueVals=trueVals, useEventAnalysis=F, defaults=defaults)
sim.2.N<-N_simPar("N", N.vals, trueVals=trueVals, useEventAnalysis=F, defaults=defaults)
sim.2.pEvent<-N_simPar("pEvent", pEvent.vals, trueVals=trueVals, useEventAnalysis=T, defaults=defaults)
sim.2.SNR<-N_simPar("SNR", SNR.vals, trueVals=trueVals, useEventAnalysis=F, defaults=defaults)

#2 Vars, fixed gammas
trueVals.fixed.X1.X2.gammas<-trueVals.fixed.X1.X2
trueVals.fixed.X1.X2.gammas[,"gamma12L"]<- -0.3
trueVals.fixed.X1.X2.gammas[,"gamma21L"]<- -0.15
trueVals.fixed.X1.X2.gammas[,"gamma12S"]<- -0.15
trueVals.fixed.X1.X2.gammas[,"gamma21S"]<- -0.3


sim.2.embedD.fixed<-N_simPar("embedD", embedD.vals, trueVals=trueVals.fixed.X1.X2.gammas, useEventAnalysis=F, defaults=defaults)
sim.2.N.fixed<-N_simPar("N", N.vals, trueVals=trueVals.fixed.X1.X2.gammas, useEventAnalysis=F, defaults=defaults)
sim.2.pEvent.fixed<-N_simPar("pEvent", pEvent.vals, trueVals=trueVals.fixed.X1.X2.gammas, useEventAnalysis=T, defaults=defaults)
sim.2.SNR.fixed<-N_simPar("SNR", SNR.vals, trueVals=trueVals.fixed.X1.X2.gammas, useEventAnalysis=F, defaults=defaults)



# Looking at series collinearity ------------------------------------------
#Random zeta, fixed etaX2 = etaX1
trueVals.fixed.X1.gammas.etaX2<-trueVals.fixed.X1.gammas
trueVals.fixed.X1.gammas.etaX2[,"etaX2"]<- -0.25

#Fixed zetam random eta
trueVals.fixed.X1.gammas.zetaX2<-trueVals.fixed.X1.gammas
trueVals.fixed.X1.gammas.zetaX2[,"zetaX2"]<- -0.15





# 4 Vars ------------------------------------------------------------------
defaults$nVars<-4
plotParams<-names(omxGetParameters(model))[grepl("eta",names(omxGetParameters(model)))][1:(defaults$nVars*2)]

sim.4.N<-N_simPar("N", N.vals, trueVals=trueVals, useEventAnalysis=F, defaults=defaults)
sim.4.embedD<-N_simPar("embedD", embedD.vals, trueVals=trueVals, useEventAnalysis=F, defaults=defaults)
sim.4.pEvent<-N_simPar("pEvent", pEvent.vals, trueVals=trueVals, useEventAnalysis=F, defaults=defaults)
sim.4.SNR<-N_simPar("SNR", SNR.vals[c(1,3)], trueVals=trueVals, useEventAnalysis=F, defaults=defaults)



