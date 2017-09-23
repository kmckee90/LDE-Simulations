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
source("lib/NL/N_simPar.R")
source("lib/MS/MS_simPar.R")
source("lib/MS/comparison/simPar.R")

mxOption(NULL, 'Default optimizer', 'NPSOL')
options(width=80)
options(scipen=999)

# Settings ----------------------------------------------------------------
genPlotting=T
par(mfrow=c(2,1), mai=c(0.4,1,0.4,1))

nTrials<- 50
#Generate true values
randEta1<-matrix(runif(nTrials,  -0.9, -0.2))
randZeta1<-matrix(runif(nTrials, -0.5, -0.1 ))
randEta2<-matrix(runif(nTrials,  -0.04, -0.01))
randZeta2<-matrix(runif(nTrials, -0.03, -0.001))
trueVals<-cbind(randEta1, randZeta1, randEta2, randZeta2)
colnames(trueVals)<-c("eta", "zeta", "etaM", "zetaM")

#Generate noise for equivalent series:
source("lib/MS/MS_defaults.R")
tESD<-1/tSNR
noise<-rnorm(N, mean=0, sd=tESD)


# Compare Model original vs duplicate within-row --------------------------

#With broad model
source("lib/MS/MS_defaults.R")
source("lib/MS/MS_genData.old.R")
model$data<-mxData(tEmbedded, type="raw")
model<-mxRun(model)
model<-mxTryHard(model, extraTries=15, greenOK=T)
summary(model)

## Plot example series
# pdf("results/Multiscale Continuous Series.pdf", height=6, width=9)
# par(mfrow=c(2,1), mai=c(0.4,1,0.4,1))
# for(i in 9:11){
#   initCond.X<<- runif(1, 1.5, 4)*sample(c(-1,1),1)
#   initCond.Y<<- runif(1, 0.5, 1)*sample(c(-1,1),1)
#   etaX1<<-trueVals[i,'eta']
#   etaX2<<-trueVals[i,'etaM']
#   zetaX1<<-trueVals[i,'zeta']
#   zetaX2<<-trueVals[i,'zetaM']
#   gammaX1<<- -etaX1
#   source("lib/MS/MS_genData.old.R")
# }
# dev.off()

MS_sim<-MS_simPar("pEvent", 0, trueVals)

  # save(MS_sim, file="MS_sim_500_3x7_c.Rdata")


#Plot results

load("MS_sim_500_3x7_c.Rdata")

pdf("results/Multiscale_Continuous 3x7.pdf", height=6, width=9)
  par(mfrow=c(1,1), mai=c(1,1,1,1))
  
  #Mult-scale basic sim
  plot(MS_sim$trueVals[,'eta'], MS_sim$Estimates[,'eta'],
       # ylim=c(-0.05, 0.0),      xlim=c(-0.05, 0.0),
       # ylim=c(min(MS_sim$Estimates[,i]), 0),      xlim=c(min(MS_sim$Estimates[,i]), 0),
       col="black",
       main="Narrow Eta Estimates: Multiscale",
       ylab="Estimates",
       xlab="True Values")
  lines(seq(-10,10,0.1),seq(-10,10,0.1), col="green")
  abline(lm(MS_sim$Estimates[,'eta']~MS_sim$trueVals[,'eta']), col="red", lty=2, lwd=2)
  
  plot(MS_sim$trueVals[,'zeta'], MS_sim$Estimates[,'zeta'],
       # ylim=c(-0.05, 0.0),      xlim=c(-0.05, 0.0),
       # ylim=c(min(MS_sim$Estimates[,i]), 0),      xlim=c(min(MS_sim$Estimates[,i]), 0),
       col="black",
       main="Narrow Zeta Estimates: Multiscale",
       ylab="Estimates",
       xlab="True Values")
  lines(seq(-10,10,0.1),seq(-10,10,0.1), col="green")
  abline(lm(MS_sim$Estimates[,'zeta']~MS_sim$trueVals[,'zeta']), col="red", lty=2, lwd=2)
  
  qc.etaM<-MS_sim$Estimates$etaM > -0.07 & MS_sim$Estimates$etaM < -0.005
  qc.etaM<-T
  
  plot(MS_sim$trueVals[qc.etaM,'etaM'], MS_sim$Estimates[qc.etaM,'etaM'],
       ylim=c(-0.05, 0.0),      xlim=c(-0.05, 0.0),
       # ylim=c(min(MS_sim$Estimates[,i]), 0),      xlim=c(min(MS_sim$Estimates[,i]), 0),
       col="black",
       main="Broad Eta Estimates: Multiscale",
       ylab="Estimates",
       xlab="True Values")
  lines(seq(-10,10,0.1),seq(-10,10,0.1), col="green")
  abline(lm(MS_sim$Estimates[qc.etaM,'etaM']~MS_sim$trueVals[qc.etaM,'etaM']), col="red", lty=2, lwd=2)
  
  qc.zetaM<-MS_sim$Estimates$zetaM > -0.07 & MS_sim$Estimates$zetaM < -0.002
  qc.zetaM<-T
  plot(MS_sim$trueVals[qc.zetaM,'zetaM'], MS_sim$Estimates[qc.zetaM,'zetaM'],
       # ylim=c(-0.05, 0.0),      xlim=c(-0.05, 0.0),
       # ylim=c(min(MS_sim$Estimates[,i]), 0),      xlim=c(min(MS_sim$Estimates[,i]), 0),
       col="black",
       main="Broad Zeta Estimates: Multiscale",
       ylab="Estimates",
       xlab="True Values")
  lines(seq(-10,10,0.1),seq(-10,10,0.1), col="green")
  abline(lm(MS_sim$Estimates[qc.zetaM,'zetaM']~MS_sim$trueVals[qc.zetaM,'zetaM']), col="red", lty=2, lwd=2)
  

dev.off()






lm( MS_sim$Estimates$eta  ~ MS_sim$trueVals[,'eta'] )
lm( MS_sim$Estimates$zeta ~ MS_sim$trueVals[,'zeta'] )
lm( MS_sim$Estimates$etaM  ~ MS_sim$trueVals[,'etaM'] )
lm( MS_sim$Estimates$zetaM ~ MS_sim$trueVals[,'zetaM'] )


