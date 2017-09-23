
# Compile results ---------------------------------------------------------
# load("results/simResults.RData")
# s.tSNR.vals<-rev(c(1/4, 1/3, 1/2, 1, 2, 3, 4))
# sEA.tSNR.vals<-rev(c(1/4, 1/3, 1/2, 1, 2, 3, 4))

source("plotSim.R")
load("results/results.RData")

pdf("results/LDESimulation_EventAnalysis.pdf", height=6, width=9)

# N ---------------------------------------------------------------------
plotSimL(S.N$Mean.Estimates, output.sd=S.N$SD.Estimates, plotVars.S, "topright",trueVals = T)
title(main="Mean Univariate Dynamic Estimates by Series Length", 
      sub="100 runs per condition")

plotSimL(S.N$Median.Estimates, output.sd=S.N$SD.Estimates, plotVars.S, "topright",trueVals = T)
title(main="Median Univariate Dynamic Estimates by Series Length", 
      sub="100 runs per condition")

plotSimL(C.N$Mean.Estimates, output.sd=C.N$SD.Estimates, plotVars.C, "topright",trueVals = T)
title(main="Mean Bivariate Dynamic Estimates by Series Length", 
      sub="100 runs per condition")

plotSimL(C.N$Median.Estimates, output.sd=C.N$SD.Estimates, plotVars.C, "topright",trueVals = T)
title(main="Median Bivariate Dynamic Estimates by Series Length", 
      sub="100 runs per condition")


# pEvent ---------------------------------------------------------------------
plotSimL(S.pEvent$Mean.Estimates, output.sd=S.pEvent$SD.Estimates, plotVars.S, "topright",trueVals = T)
title(main="Mean Univariate Dynamic Estimates by Event Probability", 
      sub="100 runs per condition")

plotSimL(S.pEvent$Median.Estimates, output.sd=S.pEvent$SD.Estimates, plotVars.S, "topright",trueVals = T)
title(main="Median Univariate Dynamic Estimates by Event Probability", 
      sub="100 runs per condition")

plotSimL(C.pEvent$Mean.Estimates, output.sd=C.pEvent$SD.Estimates, plotVars.C, "topright",trueVals = T)
title(main="Mean Bivariate Dynamic Estimates by Event Probability", 
      sub="100 runs per condition")

plotSimL(C.pEvent$Median.Estimates, output.sd=C.pEvent$SD.Estimates, plotVars.C, "topright",trueVals = T)
title(main="Median Bivariate Dynamic Estimates by Event Probability", 
      sub="100 runs per condition")



# SNR ---------------------------------------------------------------------
plotSimL(S.tSNR$Mean.Estimates, output.sd=S.tSNR$SD.Estimates, plotVars.S, "topright",trueVals = T)
title(main="Mean Univariate Dynamic Estimates by Signal/Noise Ratio", 
      sub="100 runs per condition")

plotSimL(S.tSNR$Median.Estimates, output.sd=S.tSNR$SD.Estimates, plotVars.S, "topright",trueVals = T)
title(main="Median Univariate Dynamic Estimates by Signal/Noise Ratio", 
      sub="100 runs per condition")

plotSimL(C.tSNR$Mean.Estimates, output.sd=C.tSNR$SD.Estimates, plotVars.C, "topright",trueVals = T)
title(main="Mean Bivariate Dynamic Estimates by Signal/Noise Ratio", 
      sub="100 runs per condition")

plotSimL(C.tSNR$Median.Estimates, output.sd=C.tSNR$SD.Estimates, plotVars.C, "topright",trueVals = T)
title(main="Median Bivariate Dynamic Estimates by Signal/Noise Ratio", 
      sub="100 runs per condition")


# Univariate Distribution Plots ------------------------------------------------------


for(i in plotVars.S){
  plotSimD(S.N, rev(N.vals), i, ylim=c(0,60), xlim=NULL, skip=1)
  abline(v=eval(parse(text=i)))
  title(main=paste("Single Oscillator: ",i," Distributions",sep=" "), 
        sub="100 runs per condition")
}
for(i in plotVars.S){
  plotSimD(S.pEvent, pEvent.vals, i, ylim=c(0,60), xlim=NULL, skip=3)
  abline(v=eval(parse(text=i)))
  title(main=paste("Single Oscillator: ",i," Distributions",sep=" "), 
        sub="100 runs per condition")
}
for(i in plotVars.S){
  plotSimD(S.tSNR, rev(tSNR.vals), i, ylim=c(0,60), xlim=NULL, skip=1)
  abline(v=eval(parse(text=i)))
  title(main=paste("Single Oscillator: ",i," Distributions",sep=" "), 
        sub="100 runs per condition")
}


for(i in plotVars.C){
  plotSimD(C.N, rev(N.vals), i, ylim=c(0,60), xlim=NULL, skip=1)
  abline(v=eval(parse(text=i)))
  title(main=paste("Coupled Oscillators: ",i," Distributions",sep=" "), 
        sub="100 runs per condition")
}
for(i in plotVars.C){
  plotSimD(C.pEvent, pEvent.vals, "zetaX2", ylim=c(0,60), xlim=NULL, skip=3)
  abline(v=eval(parse(text=i)))
  title(main=paste("Coupled Oscillators: ",i," Distributions",sep=" "), 
        sub="100 runs per condition")
}
for(i in plotVars.C){
  plotSimD(C.tSNR, rev(tSNR.vals), i, ylim=c(0,60), xlim=NULL, skip=1)
  abline(v=eval(parse(text=i)))
  title(main=paste("Coupled Oscillators: ",i," Distributions",sep=" "), 
        sub="100 runs per condition")
}

dev.off()

# Ratios of estimate means to true values:----------------------------------------------
#Single
cbind(N.vals, S.N$Mean.Estimates[,2:3]/matrix(c(etaX1,zetaX1,etaX2,zetaX2,gammaX2,gammaX1),length(N.vals),6,byrow=T))
cbind(tSNR.vals, S.tSNR$Mean.Estimates[,2:3]/matrix(c(etaX1,zetaX1,etaX2,zetaX2,gammaX2,gammaX1),length(tSNR.vals),6,byrow=T))
cbind(pEvent.vals, S.pEvent$Mean.Estimates[,2:3]/matrix(c(etaX1,zetaX1,etaX2,zetaX2,gammaX2,gammaX1),length(pEvent.vals),6,byrow=T))

#Coupled
cbind(N.vals, C.N$Mean.Estimates[,2:7]/matrix(c(etaX1,zetaX1,etaX2,zetaX2,gammaX2,gammaX1),length(N.vals),6,byrow=T))
cbind(tSNR.vals, C.tSNR$Mean.Estimates[,2:7]/matrix(c(etaX1,zetaX1,etaX2,zetaX2,gammaX2,gammaX1),length(tSNR.vals),6,byrow=T))
cbind(pEvent.vals, C.pEvent$Mean.Estimates[,2:7]/matrix(c(etaX1,zetaX1,etaX2,zetaX2,gammaX2,gammaX1),length(pEvent.vals),6,byrow=T))

