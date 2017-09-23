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

nTrials<- 200
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




#Run SimPar for the normal model
# source("lib/NL/N_defaults.R")
# embedRange<-c(1, 2, 3, 4)
# 
# embedD2<<-1
# out_1<-NLvsMS_simPar("embedD2", embedD2, trueVals[,1:2], exper=T, noise=noise)
# 
# embedD2<<-2
# out_2<-NLvsMS_simPar("embedD2", embedD2, trueVals[,1:2], exper=T, noise=noise)
# 
# embedD2<<-3
# out_3<-NLvsMS_simPar("embedD2", embedD2, trueVals[,1:2], exper=T, noise=noise)
# 
# embedD2<<-4
# out_4<-NLvsMS_simPar("embedD2", embedD2, trueVals[,1:2], exper=T, noise=noise)

# results<-list(out_1, out_2, out_3, out_4)
# save(results, file="MultiEmbed_Comparison.RData")


#With broad model:
# 


source("lib/MS/MS_defaults.R")
source("lib/MS/MS_genData.old.R")
model$data<-mxData(tEmbedded, type="raw")
model<-mxRun(model)
model<-mxTryHard(model, extraTries=15, greenOK=T)
summary(model)


# MS_sim<-MS_simPar("pEvent", 0, trueVals)

  # save(MS_sim, file="MS_sim_200_4x3.Rdata")

embedD2<<-1
out_1_B<-NLvsMS_simPar("embedD2", embedD2, trueVals, noise=noise)

embedD2<<-3
out_2_B<-NLvsMS_simPar("embedD2", embedD2, trueVals, noise=noise)

# embedD2<<-4
# out_3_B<-NLvsMS_simPar("embedD2", embedD2, trueVals, noise=noise)

# 
# embedD2<<-3
# out_3<-NLvsMS_simPar("embedD2", embedD2, trueVals[,1:2], exper=T, noise=noise)
# 
# embedD2<<-4
# out_4<-NLvsMS_simPar("embedD2", embedD2, trueVals[,1:2], exper=T, noise=noise)

# results<-list(out_1, out_2, out_3, out_4)
# save(results, file="MultiEmbed_Comparison.RData")




#Plot results

load("MultiEmbed_Comparison.RData")
 