#SIMULATION SCRIPT
#Test various parameters for LDE fitting
#By Kevin McKee
#setwd("\\\\Kevin-pc/vipbg/Spring 2017/simulations/ParameterRecovery/EventAnalysisV4")

library(deSolve)
library(OpenMx)
library(psych)
library(MASS)

source("../GLLAfunctions.R")
source("fitODE.R")
source("EventAnalysis_trainObj.R")
source("simParEA.R")
source("plotSim.R")
source("klmFuncs.R")
source("N_defaults_basic_events.R")

mxOption(NULL, 'Default optimizer', 'NPSOL')
options(width=80)
options(scipen=999)

# Settings ----------------------------------------------------------------
genPlotting=F

#  ------------------------------------------------------------------------
# Through random simulations, get data to determine ideal linear mixed
# objective weights 
#  ------------------------------------------------------------------------
normTest.dat<-NULL


for(i in 1:100){
  
tSNR<-sample(2:8, 1)
N<-sample(30:60, 1)
nVars<-sample(2,1)
eventScale<-runif(1,min=1,max=5)
pEventTotal<-runif(1, min=0.02, max=0.15)

indDyn<- matrix(c(   -.6,   -.3,
                      -.3,   -.15,
                      -.2,   -.2,
                      -.1,   -.1),
                 4, 2, byrow=T)
gammas<-matrix(c(   1,   0.05,   -0.2,     -0.1,
                     -0.3,   1,     0,        0.2,
                     -.05,  -0,     1,       -0.2,
                     -.05,  -0.2,  -0.2,        1  )
                
                ,4,4, byrow=T)

indDyn<-matrix(indDyn[1:nVars,], nVars, 2, byrow=T)
gammas<-matrix(gammas[1:nVars,1:nVars], nVars, nVars, byrow=T)
initCond<-matrix(0,nVars, 2)
initCond[,1]<-rnorm(nVars)
initCond[,2]<-rnorm(nVars)
source("N_matrices.R")
source("N_model.R")
source("N_genData.R")

NL<-eventAnalysis_trainObj(tEmbedded, model, threshold.sv = -.03,
                  optim=T, plotting=F, tries=4, min.data = 15)

write.table(normTest.dat, file="output/normTest.dat.2.csv", sep=",", append=T, row.names=F, col.names= F)
}



