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
source("N_defaults_basic_events.R")

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
nVars<-sample(c(1),1)
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
                  optim=T, plotting=T, tries=4, min.data = 15)

write.csv(normTest.dat, file="output/normTest.dat.csv", append=T, row.names=F, col.names= F)
}



