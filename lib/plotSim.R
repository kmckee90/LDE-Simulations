

plotSimRandEffect<-function(simOutput, pars, ylim=NULL){
  
  # Plot estimates ----------------------------------------------------------
  vals<-unique(simOutput$Estimates[,1])
  parName<-colnames(simOutput$Estimates)[1]
  trueVals<-as.data.frame(simOutput$trueVals[,pars])
  estimates<-as.data.frame(simOutput$Estimates[,pars])
  
  Est.resid<- estimates - matrix(1,length(vals),1)%x%simOutput$trueVals[,pars]
  Est.resid<-cbind(simOutput$Estimates[,parName], Est.resid)
  colnames(Est.resid)[1]<-parName
  
  sumStats<-list()
  for(h in pars){
    sumStats[[h]]<-list()
    for(i in vals){
      sumStats[[h]]<-rbind(sumStats[[h]], summary(Est.resid[,h][Est.resid[,parName]==i]))
    }
  }
  
  # Scatter plot ------------------------------------------------------------
  valColors<-rainbow(length(vals), s=1, v=0.75, alpha=1)
  for(i in pars){
    
  for(h in 1:length(vals)){
  
      plot(c(min(trueVals[,i]), max(trueVals[,i])),c(0,0), type="n",
           main=paste0(i," Estimates by True Values, ",parName,"=",vals[h]),
           sub = paste(length(Est.resid[,1])/length(vals),"simulated sets"),
           xlab="True Parameter Values",
           ylab="Estimated Parameters")
      points(trueVals[,i], estimates[simOutput$Estimates[,1]==vals[h], i], pch=1, lwd=1)#, col=valColors[h])
    
    
      abline(a=0,b=1, col="green")
      abline(lm(estimates[simOutput$Estimates[,1]==vals[h], i]~trueVals[,i]), col="red", lty=2)
      corVals<-cor(trueVals[,i],estimates[simOutput$Estimates[,1]==vals[h], i], use="pairwise.complete")
  }    
    }

  for(h in vals){
    corVals<-cor(trueVals[,'etaX1'], trueVals[,'zetaX1'], use="pairwise.complete")

    plot(trueVals[,'etaX1'], trueVals[,'zetaX1'],
         main=paste0("Eta x Zeta Estimates, ",parName,"=",h),
         sub = paste(length(Est.resid[,1])/length(vals),"simulated sets"),
         xlab="Estimated Eta",
         ylab="Estimated Zeta",
         col="darkgray",
         ylim=c(-0.75,0),
         xlim=c(-0.75,0))
    points(estimates[simOutput$Estimates[,1]==h,'etaX1'], estimates[simOutput$Estimates[,1]==h, 'zetaX1'],
           col="black", lwd=1)

  }
  
  plotColors<-rainbow(length(pars), s=1, v=0.75, alpha=1)
  plotColors.muted<-rainbow(length(pars),  s=1, v=0.75,  alpha=0.5)
  
  if(is.null(ylim)){
    bounds<-c(-max(abs(Est.resid[,-1])), max(abs(Est.resid[,-1])))
  }else{
    bounds<-ylim
  }

  plot(c(min(vals), max(vals)), c(0,0), type="l", lwd=1, col="black",
       ylim=bounds,
       main=paste("Parameter Estimate Error by",parName),
       sub = paste(length(Est.resid[,1])/length(vals),"simulated sets each"),
       xlab=paste(parName,"value"),
       ylab="Estimated Parameter - True Parameter")
  abline(h=0,lwd=2)
  for(h in 1:length(sumStats)){
    lines(vals, sumStats[[h]][,'Median'], col=plotColors[h], lwd=2)
    lines(vals, sumStats[[h]][,'1st Qu.'], col=plotColors[h], lwd=1, lty=2)
    lines(vals, sumStats[[h]][,'3rd Qu.'], col=plotColors[h], lwd=1, lty=2)
  }
  
  legend("topright",
         legend=pars,
         col=plotColors, 
         lty=1, lwd=2)
  
}





compareSE<-function(simOutput, pars, ylim=NULL, xlim=NULL){
  vals<-unique(simOutput$Estimates[,1])
  parName<-colnames(simOutput$Estimates)[1]
  
  for(h in pars){
    for(j in vals){
      eta.range<-range(simOutput$Estimates[,h])
      bins.seq<-seq(eta.range[1], 0, 0.05)
      bins<-cbind(bins.seq, c(bins.seq[-1],0))
      
      SE<-NULL
      for(i in 1:(nrow(bins)-1)){
        this.bin<-simOutput$Estimates[,h]>=bins[i,1] & simOutput$Estimates[,h] < bins[i,2]
        this.sd<-sd(simOutput$Estimates[this.bin & simOutput$Estimates[,1]==j, h], na.rm=T)
        this.se<-mean(simOutput$Std.Errs[this.bin & simOutput$Estimates[,1]==j, h], na.rm=T)
        SE<-rbind(SE, c("Bin"=bins[i], "sd"=this.sd, "se"=this.se))
      }
      SE<-as.data.frame(SE)
      
      plot(SE$Bin, SE$se, type="l", ylim=ylim,xlim=xlim,
           main=paste(h, "Empirical vs Estimated Standard Error,",parName,"=",j),
           xlab="Parameter value",
           ylab="Standard Error")
      lines(SE$Bin, SE$sd, type="l", lty=2)
      legend("topright",
             legend=c("Estimated","Empirical"),
             lty=c(1,2), lwd=1)
    }
  }
}

