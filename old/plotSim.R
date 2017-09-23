
# Plot estimates and std errors from simulation ---------------------------

# Line plot ---------------------------------------------------------------
plotSimL<-function(output.m, output.sd=NULL, plotVars, coords, trueVals=F){
  parname<-colnames(output.m)[1]
  plot(x=output.m[,parname],y=output.m[,2],type="l", lwd=2, ylim=c(-1,1), col=1,
       #main=paste("Mean Parameter Estimates by ",parname,", ",nTrials," trials",sep=""),
       xlab=parname,
       ylab="Estimate"
  )
  for(i in 1:length(plotVars)){
    lines(output.m[,parname], output.m[,plotVars[i]], lwd=2, type="l", col=i)
  }
  
  #Conf Intervals
  if(!is.null(output.sd)){
    for(i in 1:length(plotVars)){
      lines(output.m[,parname], output.m[,plotVars[i]]+1.96*output.sd[,plotVars[i]], lwd=1, type="l", col=i, lty=3)
      lines(output.m[,parname], output.m[,plotVars[i]]-1.96*output.sd[,plotVars[i]], lwd=1, type="l", col=i, lty=3)
      }
  }
  
  if(trueVals==T){
    for(i in 1:nVars){
      abline(h=indDyn[i,1], lty=2, lwd=2, col=i*1)
      abline(h=indDyn[i,2], lty=2, lwd=2, col=i*2)
    }
  # lines(output.m[,parname],output.m[,parname],lty=2, lwd=2, col=7)
  abline(h=0, lty=1 )
  }
  # abline(h=eval(parse(text=paste("c(", paste(plotVars, collapse=","), ")"))), lty=2, col=1:length(plotVars))
  
  legend(coords,
         plotVars,
         lwd=2, col=c(1:length(plotVars)))
}




# Density plots -----------------------------------------------------------


plotSimD<-function(output, vals, plotVar, ylim=NULL, xlim=NULL, skip=1){
  # Zeta --------------------------------------------------------------------
  output.d<-list()
  output.d[[1]]<-density(output$Estimates[output$Estimates[,1]==vals[1],plotVar], na.rm=T)
  parname<-colnames(output$Estimates)[1]
  
  plot(output.d[[1]], type="l", 
       main=" ",
       xlab=plotVar,
       xlim=xlim,
       ylim=ylim)
  cols<-rainbow(length(vals),start=0.5,end=1,alpha=1)
  for(i in 1:length(vals)){
    output.d[[i]]<-density(output$Estimates[output$Estimates[,1]==vals[i],plotVar], na.rm=T)
    # hist(sEA.pEvent$Estimates[sEA.pEvent$Estimates[,1]==vals[i],3], breaks=24)
    lines(output.d[[i]], type="l", lwd=1, col=cols[i])
    # abline(v=sEA.pEvent$Mean.Estimates[i,3], col=cols[i], lty=3)
  }
  # abline(v=eta, lty=1)
  everyOther<-(1:(length(vals)/skip))*skip
  legend("topright",legend=vals[everyOther],col=cols[everyOther], lty=1, lwd=3, 
         title=colnames(output$Estimates)[1])
}