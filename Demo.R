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
embedD.vals<-c(3,4,5,6,7,8)
N.vals<- c(45, 90, 135, 180, 225, 270, 315, 360)  [1:4]
pEvent.vals<- round(seq(0, 0.2, 0.02),2)               #[1:6]
SNR.vals<- c(1/6, 1/4, 1/2, 1, 2, 3, 4, 5, 6, 8)         [c(3,4,5,7,10)]

# Defaults, parameters, and initial conditions ------------------------------
defaults<-list(    
  "nTrials"=100,
  
  "nVars"= 4,
  
  "N"=60,
  "SNR"=6,
  "pEvent"=0.05,
  "embedD"=4,
  "eventScale"=1,
  "eventNormal"=F,
  
  "theTau" = 1,
  "deltaT" = 1,
  "tryHard.tries" = 20,
  
  "plotting"=T
)

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




# Run demo of data generation and algorithm -------------------------------


defaults$nVars<-1
source("lib/NL/N_matrices.R")
source("lib/NL/N_model.R")
source("lib/NL/N_constraints.R")

dat<-N_genData(defaults$N, SNR=defaults$SNR, numVars=defaults$nVars, trueVals[1,], initCond[1,], 
               events=0.05, eventNormal=defaults$eventNormal, eventType="level", eventScale = 1,
               plotting=defaults$plotting)

tEmbedded<-NULL
for(i in 1:defaults$nVars){
  tEmbedded <- cbind(tEmbedded, gllaEmbed(dat[,i], embed=defaults$embedD, tau=defaults$theTau, idColumn=FALSE))
}
# print(head(tEmbedded))
colnames(tEmbedded)<-manifestVars
tEmbedded<-cbind(tEmbedded, "Occasion"=1:nrow(tEmbedded))

EA<-eventAnalysis(tEmbedded, model, plotting=defaults$plotting, prediction=T, tries=defaults$tryHard.tries, min.data = 20*defaults$nVars)
EA<-eventAnalysis(tEmbedded, model, plotting=defaults$plotting, prediction=T, tries=defaults$tryHard.tries, min.data = 50*defaults$nVars, use=EA$info)



