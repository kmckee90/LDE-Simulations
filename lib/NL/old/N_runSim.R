#SIMULATION SCRIPT
#Test various parameters for LDE fitting
#By Kevin McKee

library(deSolve)
library(OpenMx)
library(psych)
library(MASS)

source("../GLLAfunctions.R")
source("fitODE.r")
source("EventAnalysis.r")
source("simParEA.r")
source("plotSim.R")
source("klmFuncs.r")
source("expectedValsML.r")

mxOption(NULL, 'Default optimizer', 'NPSOL')
options(width=80)
options(scipen=999)

# Settings ----------------------------------------------------------------
genPlotting=T
par(mfrow=c(2,1), mai=c(0.4,1,0.4,1))

nTrials<- 100
N.vals<- c(45, 90, 135, 180, 225, 270, 315, 360)  [1:3]
pEvent.vals<- round(seq(0, 0.3, 0.02),2)               #[1:2]
tSNR.vals<- c(1/6, 1/4, 1/2, 1, 2, 3, 4, 5, 6, 8)         [c(3,4,5,7,9,10)]



# Multiscale model test ---------------------------------------------------
results.1<-NULL
for(i in 1:20){
  source("MS_defaults_basic.R")
  source("MS_matrices.R")
  source("N_model.R")
  source("MS_genData.old.R")
  model<-omxSetParameters(model, labels="gamma", newlabels = "eta")
  model<-omxSetParameters(model, labels="intX1", values=0, free=F)
  model<-omxSetParameters(model, labels="slopeX1", values=0, free=F)
  
  model$data<-mxData(tEmbedded, type="raw")
  model<-mxRun(model)
  model<-mxTryHard(model, extraTries = 40, greenOK=T)
  # results.1<-rbind(results.1, as.vector(omxGetParameters(model)[1:4]))
  summary(model)
}


#One method: calculate similar to a coupled system (current version)
#Two-step: Run broad LDE, run narrow LDE on the residuals
#Other method using narrow LDE means?


# source("MS_defaults_basic.R")
# source("MS_genData.old.R")
# source("N_defaults_basic_events.R")



NL<-eventAnalysis(tEmbedded, model, threshold.sv = -.03,
                  optim=F, plotting=T, prediction=F, tries=2, min.data = 15*nVars)

NL2<-eventAnalysis(tEmbedded, model, threshold.sv = 0.006,
                   optim=F, plotting=T, prediction=F, tries=5, min.data = 20*nVars, use=NL$info)
summary(NL2$fit)


lines(shocksX^2, type="h",lwd=1,col="blue")




NL.s<-expectedValsML(NL)


tData<-tOscData
tEmbedded<- gllaEmbed(tData, embed=embedD, tau=theTau, idColumn=FALSE)
colnames(tEmbedded)<-manifestVars
Occasion<-1:nrow(tEmbedded) + round(embedD/2)
tEmbedded<-cbind(tEmbedded, Occasion)


# #For simple testing disable slope and intercept estimation:
# model<-omxSetParameters(model, labels="intX1", values=0, free=F)
# model<-omxSetParameters(model, labels="slopeX1", values=0, free=F)
# 
# #Disable the broad LDE
# model$L$values[(3*embedD2*embedD+1):length(model$L$values)]<-0
# model<-omxSetParameters(model, labels="etaM", values=0, free=F)
# model<-omxSetParameters(model, labels="zetaM", values=0, free=F)
# model<-omxSetParameters(model, labels="gamma", values=1, free=F)
# model<-omxSetParameters(model, labels="VXM", values=0, free=F)
# model<-omxSetParameters(model, labels="VdXM", values=0, free=F)
# model<-omxSetParameters(model, labels="Vd2XM", values=0, free=F)
# model<-omxSetParameters(model, labels="rVXM_VdXM", values=0, free=F)


# model$data<-mxData(tEmbedded, type="raw")
# model<-mxRun(model)
# model<-mxTryHard(model, extraTries = 30, greenOK=T)
# summary(model)








# Run Simulations ---------------------------------------------------------
# source("N_defaults_basic_events.R")
# plotVars<-names(omxGetParameters(model))[1:2]
# NL.pEvent<-simParEA("pEvent", seq(0, .2, 0.025))


source("N_defaults_basic_events.R")
source("N_genData.R")
NL<-eventAnalysis(tEmbedded, model, threshold.sv = -.03,
                  optim=F, plotting=T, tries=6, min.data = 12*nVars)
NL2<-eventAnalysis(tEmbedded, model, threshold.sv = 0.0275,
                  optim=T, plotting=T, tries=4, min.data = 8*nVars, use=NL$info)
summary(NL$fit)






