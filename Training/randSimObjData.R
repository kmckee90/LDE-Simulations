#SIMULATION SCRIPT
#Test various parameters for LDE fitting
#By Kevin McKee

library(deSolve)
library(OpenMx)
library(psych)
library(MASS)

source("../GLLAfunctions.R")
source("fitODE.r")
source("EventAnalysis_trainObj.r")
source("simParEA.r")
source("plotSim.R")
source("klmFuncs.r")
source("expectedValsML.r")
source("N_defaults_basic_events.R")

mxOption(NULL, 'Default optimizer', 'NPSOL')
options(width=80)
options(scipen=999)

# Settings ----------------------------------------------------------------
genPlotting=F
par(mfrow=c(1,1), mai=c(1,1,1,1))

nTrials<- 100
N.vals<- c(45, 90, 135, 180, 225, 270, 315, 360)  [1:3]
pEvent.vals<- round(seq(0, 0.3, 0.02),2)               #[1:2]
tSNR.vals<- c(1/6, 1/4, 1/2, 1, 2, 3, 4, 5, 6, 8)         [c(3,4,5,7,9,10)]



#  ------------------------------------------------------------------------
# Through random simulations, get data to determine ideal linear mixed
# objective weights 
#  ------------------------------------------------------------------------
normTest.dat<-NULL

for(i in 1:3){
tSNR<-sample(2:8, 1)
N<-sample(25:75, 1)
nVars<-sample(c(1,2),1)
eventScale<-runif(1,min=1,max=5)
pEventTotal<-runif(1, min=0.02, max=0.2)

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
                  optim=T, plotting=F, tries=20, min.data = 8)
}

save(normTest.dat, file="output/normTest.dat.rand.1.RData")


