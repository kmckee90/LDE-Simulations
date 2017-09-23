
#------------------------------------------------------------------------\
# Event Analysis v5
# By Kevin McKee
# 05/27/2017
#
#------------------------------------------------------------------------/
load("lib/ObjectiveWeightsLM.RData")

# Model-based Event Detection ---------------------------------------------
eventAnalysis<-function(tData, model,
                        threshold.sv=2, min.data=8,optim=T, 
                        plotting=T, prediction=T, tries=3, use=NULL
){
  
  tData<-as.matrix(tData)
  embedD<-as.numeric(model$embedD$values)
  nVars<-as.numeric(model$nVars$values)
  model$data<-mxData(tData, type="raw")
  selCol<-round(seq(1, nVars*embedD, embedD))
  
  
  # Threshold-based data exclusion func -------------------------------------
  cutData<-function(thres){
    dat<-tData
    pointInfo.s.exc<-which(pointInfo.s>thres, arr.ind=T)
    if(length(pointInfo.s.exc)>0){
      pointInfo.s.exc[,2]<-selCol[pointInfo.s.exc[,2]]
      for(i in 0:(embedD-1)){
        pointInfo.s.exc.i<-pointInfo.s.exc
        pointInfo.s.exc.i[,2]<-pointInfo.s.exc[,2]+i
        dat[pointInfo.s.exc.i]<-NA
      }
      exclusions<-list()
      for(i in 1:nVars){
        exclusions[[i]]<-pointInfo.s.exc[pointInfo.s.exc[,2]==selCol[i], 1]
      }
    }else{
      exclusions<-rep(list(NULL),nVars)
    }
    model$data<-mxData(dat, type="raw")
    model<-mxRun(model, silent=T)
    if(model$output$status$code >= 6 && tries > 0){
      model<-mxTryHard(model, extraTries = tries, greenOK=T)
    }
    return(list("fit"=model, "exclusions"=exclusions))
  }
  
  
  
  #Multiplotting!
  plotEA<-function(num=NULL){
    
    par(mfrow=c(nVars+1,1), mai=c(0.3, 1, 0.025, 0.5))
    plotColors<-rainbow(nVars, s=1, v=1, alpha=1)
    plotColors.muted<-rainbow(nVars, s=0.5, v=0.25, alpha=1)
    varCols<-round(seq(1, nVars*embedD, embedD))
    
    #Variable graph
    for(i in 1:nVars){
      plot(tData[,varCols[i]], type="l", lwd=1, col=plotColors.muted[i], xaxt='n', xlab=NULL, ylab=paste0("Series X",i))
      abline(v=out$exclusions[[i]], col="orange", lty=3)
      if(prediction){
        predictDat<-fitODE(tData, out)
        lines(1:nrow(tData), predictDat$Prediction[,i], type="l", lwd=1, col=plotColors[i])
        abline(predictDat$Means[i,1], predictDat$Means[i,2], col=plotColors[i], lty=1, lwd=1)
      }
    }
    
    #-2LL graph
    plot(c(0,nrow(tData)), c(0,0), type="n", xlab=NULL, ylab="Marginal -2LL", ylim=c(min(pointInfo.s, na.rm=T), max(pointInfo.s, na.rm=T)))
    for(i in 1:nVars){
      lines(pointInfo.s[,i], type="l", col=plotColors[i])
      abline(v=out$exclusions[[i]], col=plotColors[i], lty=3)
      abline(h=num)
    }
  }
  
  
  # Main algorithm ----------------------------------------------------------
  cat("\nGetting best start values...\r")
  if(is.null(use)){
    # if(plotting){plotEA()}
    model$data<-mxData(tData, type="raw")
    model<-mxRun(model)
    model<-mxTryHard(model, greenOK=T, extraTries=29)
    
    # Get -2LL or AIC of refsample with each point ----------------------------
    pointInfo<-matrix(NA,nrow(tData),nVars)
    
    for(h in 1:nVars){
      for(i in 1:nrow(tData)){
        cat(paste0("\rCalculating likelihood for timepoint ",i," of X",h,"..."))
        tData.drop<-tData
        
        for(j in 0:(embedD-1)){
          tData.drop[i,selCol[h]+j]<-NA
        }
        
        # print(tData.drop)
        model$data<-mxData(tData.drop, type="raw")
        model<-mxRun(model, silent=T)
        
        if(model$output$status$code >= 6 && tries > 0){
          model<-mxTryHard(model, extraTries = tries, greenOK=T)
        }
        
        pointInfo[i,h] <- -1*model$output$Minus2LogLikelihood
      }
      
      #Standardize
      pointInfo[,h]<- pointInfo[,h] - min(pointInfo[,h],na.rm=T)
      pointInfo[,h]<- pointInfo[,h] / sum(pointInfo[,h],na.rm=T)
    }
    pointInfo.s<-as.matrix(pointInfo)
  }else{
    pointInfo.s<-use
  }
  
  
  
  
  
  # Data cleaning and analysis ----------------------------------------------
  # fileid<<-0
  if(optim){
    cat("\nFinding exclusion threshold...\r")
    info.sorted<-sort(as.vector(pointInfo.s), decreasing=T)
    cuts.norm<<-NULL
    cutProp<-length(info.sorted)- min.data
    
    for(i in 1:cutProp){
      out<-cutData(info.sorted[i])
      
      #Norms to use in determining threshold
      i.UX<-mean(omxGetParameters(out$fit)[grepl("Ux",names(omxGetParameters(out$fit)))], na.rm=T)
      i.par.std<-mean(out$fit$output$standardErrors[grepl("zeta",rownames(out$fit$output$standardErrors))], na.rm=T)
      percentExclude<-length(unique(unlist(out$exclusions))) / nrow(tData)
      
      #Mixed norm as sum of norms, each weighted a priori
      norms.row<-matrix(c(i.UX, i.par.std, percentExclude),1,3)
      colnames(norms.row)<-c("UX","SEzeta","fEvents")
      norms.row<-(norms.row-ObjEA$objMeans)/sqrt(ObjEA$objVars)
      this.norm<- norms.row %*% as.matrix(ObjEA$objWeights^2)
      
      #Optional plotting of data cutting process
      if(plotting){plotEA(num = info.sorted[i])}
      
      cuts.norm<<-rbind(cuts.norm, this.norm)
    }
    i.opt<-which(cuts.norm==min(cuts.norm, na.rm=T)[1])
    out<-cutData(info.sorted[i.opt])
    i.thres<-info.sorted[i.opt]
    
  }else{
    out<-cutData(threshold.sv)
    i.thres<-threshold.sv
  }
  
  # Final plot
  if(plotting){
    # if(prediction){fitODE(tData, out)}
    plotEA(num=i.thres)
  }
  
  # Output ------------------------------------------------------------------
  return(list("fit"=out$fit,"exclusions"=out$exclusions, "info"=pointInfo.s))
}














