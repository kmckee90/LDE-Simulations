#SIMULATION SCRIPT
#Test various parameters for LDE fitting
#By Kevin McKee

library(deSolve)
library(OpenMx)
library(psych)
library(MASS)

source("lib/GLLAfunctions.R")
source("lib/fitODE.R")
source("lib/EventAnalysis.R")
source("lib/plotSim.R")
source("lib/klmFuncs.R")
source("lib/expectedValsML.R")
source("lib/MS/MSD_simPar.R")

mxOption(NULL, 'Default optimizer', 'NPSOL')
options(width=80)
options(scipen=999)

# Settings ----------------------------------------------------------------
genPlotting=T
par(mfrow=c(2,1), mai=c(0.4,1,0.4,1))

nTrials<- 750
#Generate true values
randEta1<-matrix(runif(nTrials,  -0.9, -0.2))
randZeta1<-matrix(runif(nTrials, -0.5, -0.1 ))
randEta2<-matrix(runif(nTrials,  -0.04, -0.01))
randZeta2<-matrix(runif(nTrials, -0.03, -0.001))
trueVals<-cbind(randEta1, randZeta1, randEta2, randZeta2)
colnames(trueVals)<-c("eta", "zeta", "etaM", "zetaM")

#Plot some example series

pdf("results/Multiscale Discontinuous_Discontinuous.pdf", height=6, width=9)
for(i in 1:3){  
  embedD1<<-5
  N<<-7*2*5
  source("lib/MS/MSD_genData.R")  
}
for(i in 1:3){  
  embedD1<<-3
  N<<-7*3*3
  source("lib/MS/MSD_genData.R")  
}
dev.off()



source("lib/MS/MSD_defaults.R")
source("lib/MS/MSD_genData.R")
model$data<-mxData(tEmbedded, type="raw")
model<-mxRun(model)
model<-mxTryHard(model, extraTries=20, greenOK=T)
summary(model)

# 
# embEsts4<-NULL
# for(q in 3:6){
#   embedD2<<-q
#   source("lib/MS/MSD_matrices.R")
#   source("lib/MS/MS_model.R")
#   source("lib/MS/MS_constraints.R")
#   tEmbedded<- gllaEmbed(tData, embed=embedD, tau=theTau, idColumn=FALSE)
#   tEmbedded<-tEmbedded[disConts[disConts<nrow(tEmbedded)],]
#   colnames(tEmbedded)<-manifestVars
#   Occasion<-1:nrow(tEmbedded) + round(embedD/2)
#   tEmbedded<-cbind(tEmbedded, Occasion)
#   
#   model$data<-mxData(tEmbedded, type="raw")
#   model<-mxRun(model)
#   model<-mxTryHard(model, extraTries=50, greenOK=T)
#   embEsts.c<-as.matrix(model$output$estimate)
#   colnames(embEsts.c)<-q
#   embEsts4<-cbind(embEsts4, embEsts.c)
# }
# 
# embEsts3
# embEsts4
# 
# 
# embEsts3/embEsts4
# 
# embEsts3[,1:3]/embEsts3[,2:4]
# embEsts4[,1:3]/embEsts4[,2:4]
# 
# plot(embEsts3['etaM',], type="l", ylim=c(-1,0))
# plot(embEsts3['zetaM',], type="l", ylim=c(-0.1,0))
# 

#Run a sim
MSD_sim<-MSD_simPar("pEvent", 0.5, trueVals)

# save(MSD_sim, file="MSD_sim_750_3x4_forLance.RData")
# MSD_sim<-MSD_sim3

load(file="MSD_sim_1000_3x4_adj.RData")
# load(file="MSD_sim_1000_3x7.RData")



pdf("results/Multiscale Discontinuous Sim 3x21 Occasions.pdf", height=6, width=9)

MSD_sim.qc<-MSD_sim
# MSD_sim.qc$Estimates<-MSD_sim$Estimates[MSD_sim$Estimates$eta<0,]
# MSD_sim.qc$trueVals<-MSD_sim$trueVals[MSD_sim$Estimates$eta<0,]

par(mfrow=c(1,1), mai=c(1,1,1,1))


#Mult-scale basic sim
plot(MSD_sim.qc$trueVals[,'eta'], MSD_sim.qc$Estimates[,'eta'],
     # ylim=c(-0.05, 0.0),      xlim=c(-0.05, 0.0),
     # ylim=c(min(MSD_sim.qc$Estimates[,i]), 0),      xlim=c(min(MSD_sim.qc$Estimates[,i]), 0),
     col="black",
     main="Narrow Eta Estimates: Multiscale Discontinuous",
     xlab="True Value",
     ylab="Estimate")
lines(seq(-10,10,0.1),seq(-10,10,0.1), col="green", lwd=2)
abline(lm(MSD_sim.qc$Estimates[,'eta']~MSD_sim.qc$trueVals[,'eta']), col="red", lty=2, lwd=2)

plot(MSD_sim.qc$trueVals[,'zeta'], MSD_sim.qc$Estimates[,'zeta'],
     # ylim=c(-0.05, 0.0),      xlim=c(-0.05, 0.0),
     # ylim=c(min(MSD_sim.qc$Estimates[,i]), 0),      xlim=c(min(MSD_sim.qc$Estimates[,i]), 0),
     col="black",
     main="Narrow Zeta Estimates: Multiscale Discontinuous",
     xlab="True Value",
     ylab="Estimate")
lines(seq(-10,10,0.1),seq(-10,10,0.1), col="green", lwd=2)
abline(lm(MSD_sim.qc$Estimates[,'zeta']~MSD_sim.qc$trueVals[,'zeta']), col="red", lty=2, lwd=2)

qc.etaM<-MSD_sim.qc$Estimates$etaM > -0.25 & MSD_sim.qc$Estimates$etaM < 0.5
qc.zetaM<-MSD_sim.qc$Estimates$zetaM > -0.25 & MSD_sim.qc$Estimates$zetaM < 0.5

# qc.etaM<-T
# qc.zetaM<-T

plot(MSD_sim.qc$trueVals[qc.etaM,'etaM'], MSD_sim.qc$Estimates[qc.etaM,'etaM'],
     # ylim=c(-1, 0.0),      xlim=c(-1, 0.0),
     # ylim=c(min(MSD_sim$Estimates[,i]), 0),      xlim=c(min(MSD_sim$Estimates[,i]), 0),
     col="black",
     main="Broad Eta Estimates: Multiscale Discontinuous",
     xlab="True Value",
     ylab="Estimate")
lines(seq(-10,10,0.1),seq(-10,10,0.1), col="green", lwd=2)
abline(lm(MSD_sim.qc$Estimates[qc.etaM,'etaM']~MSD_sim.qc$trueVals[qc.etaM,'etaM']), col="red", lty=2, lwd=2)


plot(MSD_sim.qc$trueVals[qc.zetaM,'zetaM'], MSD_sim.qc$Estimates[qc.zetaM,'zetaM'],
     # ylim=c(-0.05, 0.0),      xlim=c(-0.05, 0.0),
     # ylim=c(min(MSD_sim$Estimates[,i]), 0),      xlim=c(min(MSD_sim$Estimates[,i]), 0),
     col="black",
     main="Broad Zeta Estimates: Multiscale Discontinuous",
     xlab="True Value",
     ylab="Estimate")
lines(seq(-10,10,0.1),seq(-10,10,0.1), col="green", lwd=2)
abline(lm(MSD_sim.qc$Estimates[qc.zetaM,'zetaM']~MSD_sim.qc$trueVals[qc.zetaM,'zetaM']), col="red", lty=2, lwd=2)


dev.off()


lm( MSD_sim$Estimates$eta  ~ MSD_sim$trueVals[,'eta'] )
lm( MSD_sim$Estimates$zeta ~ MSD_sim$trueVals[,'zeta'] )

lm( MSD_sim$Estimates$etaM[qc.etaM]  ~ MSD_sim$trueVals[qc.etaM,'etaM'] )
lm( MSD_sim$Estimates$zetaM[qc.zetaM] ~ MSD_sim$trueVals[qc.zetaM,'zetaM'] )

