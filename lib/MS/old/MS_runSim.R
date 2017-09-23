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
source("MS_simPar.r")
source("plotSim.R")
source("klmFuncs.r")
source("expectedValsML.r")

mxOption(NULL, 'Default optimizer', 'NPSOL')
options(width=80)
options(scipen=999)

# Settings ----------------------------------------------------------------
genPlotting=T
par(mfrow=c(1,1), mai=c(0.75,1,0.75,1))

nTrials<- 500
N.vals<- c(45, 90, 135, 180, 225, 270, 315, 360)  [1:3]
pEvent.vals<- round(seq(0, 0.3, 0.02),2)               #[1:2]
tSNR.vals<- c(1/6, 1/4, 1/2, 1, 2, 3, 4, 5, 6, 8)         [c(3,4,5,7,9,10)]

source("MS_defaults_basic.R")
source("MS_genData.old.R")

# model$L$values[,1:6]<-0
model$L$values[,4:9]<-0
# model$L$values[,c(1:3, 6:9)]<-0


model<-omxSetParameters(model, labels="gamma", newlabels = "eta")
model<-omxSetParameters(model, labels="intX1", values=0, free=F)
model<-omxSetParameters(model, labels="slopeX1", values=0, free=F)

# model$data<-mxData(tEmbedded, type="raw")
# model<-mxRun(model)
# model<-mxTryHard(model, extraTries = 39, greenOK=T)
# summary(model)

MS2<-eventAnalysis(tEmbedded, model, threshold.sv = .014,
                  optim=F, plotting=F, tries=30, min.data = 12*nVars)
MS_o<-eventAnalysis(tEmbedded, model, threshold.sv = .014,
                   optim=F, plotting=F, tries=51, min.data = 12*nVars, use=MS2$info)
summary(MS_o$fit)

plot(MS1$info, type="l")
lines(shocksX^2, col="red", type="h")
plot(MS2$info, type="l")
lines(shocksX^2, col="red", type="h")
plot(MS3$info, type="l")
lines(shocksX^2, col="red", type="h")
plot(MS4$info, type="l")
lines(shocksX^2, col="red", type="h")
plot(MS5$info, type="l")
lines(shocksX^2, col="red", type="h")








out<-MS_simPar("pEvent", 0)
save(out, file="MS_simPar_results.RData")

load("MS_simPar_results.RData")
out.500<-out

load("MS_simPar_results1.RData")
out.250<-out

etaM.trim<-ests$etaM>-0.1 & ests$etaM < 0.1
zetaM.trim<-ests$zetaM>-0.1 & ests$zetaM < 0.1
zeta.trim<-ests$zeta>-2
ests<-rbind(out.500$Estimates, out.250$Estimates)
truVals<-as.data.frame(rbind(out.500$trueVals, out.250$trueVals))


pdf("Multiscale.pdf", height=6, width=9)


par(mfrow=c(1,1), mai=c(1,1,1,1))
plot(c(min(theTimes), max(theTimes[1:N])), c(-4, 4),
     xlab="Time",
     ylab="Score",
     type='n',
     main=paste0("Single Multiscale Oscillator"))
lines(out1$time, tOscData1, type='l', lwd=1, col='lightblue')

lines(out2$time, out2$x/tSD1, type='l', lwd=1, col='blue')
lines(out2$time, out2$z/tSD1, type='l', lwd=1, col='red')
lines(c(min(theTimes), max(theTimes[1:N])), c(-0, 0), type='l', lty=2, col=1)

plot(tOscData1, type='l', xlab="Time", ylab="Score", ylim=c(-4,4),
     main=paste0("Single Multiscale Oscillator"))
# lines(ou

par(mfrow=c(1,1), mai=c(0.75,1,0.75,1))

plot(ests$eta, truVals$eta, xlab="Estimate", ylab="True Value", 
     main=paste0("Eta (Narrow), r=",
                 round(cor(ests$eta, truVals$eta, use="pairwise.complete"),2),
                 " bias=",round(median(truVals$eta-ests$eta),3)))
lines(truVals$eta, truVals$eta)
title(sub="750 runs of uniform random true values")


plot(ests$zeta[zeta.trim], truVals$zeta[zeta.trim], xlab="Estimate", ylab="True Value", 
     main=paste0("Zeta (Narrow), r=",
                 round(cor(ests$zeta[zeta.trim], truVals$zeta[zeta.trim], use="pairwise.complete"),2),
                 " bias=",round(median(truVals$zeta-ests$zeta),3)))
lines(truVals$zeta, truVals$zeta)
title(sub="750 runs of uniform random true values")

plot(ests$etaM[etaM.trim], truVals$etaM[etaM.trim], xlab="Estimate", ylab="True Value", 
     main=paste0("Eta (Broad), r=",
                 round(cor(ests$etaM[etaM.trim], truVals$etaM[etaM.trim], use="pairwise.complete"),2),
                 " bias=",round(median(truVals$etaM-ests$etaM),4)))
lines(truVals$etaM, truVals$etaM)
title(sub="750 runs of uniform random true values")

plot(ests$zetaM[zetaM.trim], truVals$zetaM[zetaM.trim], xlab="Estimate", ylab="True Value", 
     main=paste0("Zeta (Broad), r=",
                 round(cor(ests$zetaM[zetaM.trim], truVals$zetaM[zetaM.trim], use="pairwise.complete"),2),
                 " bias=",round(median(truVals$zetaM-ests$zetaM),4)))
lines(truVals$zetaM, truVals$zetaM)
title(sub="750 runs of uniform random true values")



plot(c(min(theTimes), max(theTimes[1:N])), c(-4, 4),
     xlab="Time",
     ylab="Score",
     type='n',
     main=paste("Single Multiscale Oscillator", sep=""))
lines(out1$time, tOscData, type='l', lwd=1, col='lightblue')
lines(c(0,shocksX/4), type="h",lwd=3,col="blue")
lines(c(0,shocksY/4), type="h",lwd=3,col="red")

lines(out1$time, out1$x/tSD, type='l', lwd=1, col='blue')
# lines(out1$time, out1$y/tSD, type='l', lwd=2, col='blue')
lines(out1$time, out1$z/tSD, type='l', lwd=1, col='red')
# lines(out1$time, out1$w/tSD, type='l', lwd=2, col='red')

lines(c(min(theTimes), max(theTimes[1:N])), c(-0, 0), type='l', lty=2, col=1)

plot(tOscData, type='l', xlab="Time", ylab="Score", ylim=c(-4,4),
     main=paste0("Single Multiscale Oscillator"))

dev.off()

# 
# # Multiscale model test ---------------------------------------------------
# results.1<-NULL
# for(i in 1:20){
  source("MS_defaults_basic.R")
  source("MS_matrices.R")
  source("N_model.R")
  source("MS_genData.old.R")
#   model<-omxSetParameters(model, labels="gamma", newlabels = "eta")
#   model<-omxSetParameters(model, labels="intX1", values=0, free=F)
#   model<-omxSetParameters(model, labels="slopeX1", values=0, free=F)
#   
#   model$data<-mxData(tEmbedded, type="raw")
#   model<-mxRun(model)
#   model<-mxTryHard(model, extraTries = 40, greenOK=T)
#   # results.1<-rbind(results.1, as.vector(omxGetParameters(model)[1:4]))
#   summary(model)
# }
# 
# 
# #One method: calculate similar to a coupled system (current version)
# #Two-step: Run broad LDE, run narrow LDE on the residuals
# #Other method using narrow LDE means?
# 
# 
# # source("MS_defaults_basic.R")
# # source("MS_genData.old.R")
# # source("N_defaults_basic_events.R")
# 
# 
# 
# NL<-eventAnalysis(tEmbedded, model, threshold.sv = -.03,
#                   optim=F, plotting=T, prediction=F, tries=2, min.data = 15*nVars)
# 
# NL2<-eventAnalysis(tEmbedded, model, threshold.sv = 0.006,
#                    optim=F, plotting=T, prediction=F, tries=5, min.data = 20*nVars, use=NL$info)
# summary(NL2$fit)
# 
# 
# lines(shocksX^2, type="h",lwd=1,col="blue")
# 
# 
# 
# 
# NL.s<-expectedValsML(NL)
# 
# 
# tData<-tOscData
# tEmbedded<- gllaEmbed(tData, embed=embedD, tau=theTau, idColumn=FALSE)
# colnames(tEmbedded)<-manifestVars
# Occasion<-1:nrow(tEmbedded) + round(embedD/2)
# tEmbedded<-cbind(tEmbedded, Occasion)
# 
# 
# # #For simple testing disable slope and intercept estimation:
# # model<-omxSetParameters(model, labels="intX1", values=0, free=F)
# # model<-omxSetParameters(model, labels="slopeX1", values=0, free=F)
# # 
# # #Disable the broad LDE
# # model$L$values[(3*embedD2*embedD+1):length(model$L$values)]<-0
# # model<-omxSetParameters(model, labels="etaM", values=0, free=F)
# # model<-omxSetParameters(model, labels="zetaM", values=0, free=F)
# # model<-omxSetParameters(model, labels="gamma", values=1, free=F)
# # model<-omxSetParameters(model, labels="VXM", values=0, free=F)
# # model<-omxSetParameters(model, labels="VdXM", values=0, free=F)
# # model<-omxSetParameters(model, labels="Vd2XM", values=0, free=F)
# # model<-omxSetParameters(model, labels="rVXM_VdXM", values=0, free=F)
# 
# 
# # model$data<-mxData(tEmbedded, type="raw")
# # model<-mxRun(model)
# # model<-mxTryHard(model, extraTries = 30, greenOK=T)
# # summary(model)
# 
# 
# 
# 
# 
# 
# 
# 
# # Run Simulations ---------------------------------------------------------
# # source("N_defaults_basic_events.R")
# # plotVars<-names(omxGetParameters(model))[1:2]
# # NL.pEvent<-simParEA("pEvent", seq(0, .2, 0.025))
# 
# 
# source("N_defaults_basic_events.R")
# source("N_genData.R")
# NL<-eventAnalysis(tEmbedded, model, threshold.sv = -.03,
#                   optim=F, plotting=T, tries=6, min.data = 12*nVars)
# NL2<-eventAnalysis(tEmbedded, model, threshold.sv = 0.0275,
#                   optim=T, plotting=T, tries=4, min.data = 8*nVars, use=NL$info)
# summary(NL$fit)
# 
# 
# 
# 
# 
# 
