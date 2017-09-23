
#------------------------------------------------------------------------\
# Event Analysis v5
# By Kevin McKee
# 05/27/2017
#
#------------------------------------------------------------------------/

# Model-based Event Detection ---------------------------------------------
eventAnalysis_trainObj<-function(tData, model,
                        threshold.sv=2, min.data=8,optim=T, 
                        plotting=T, tries=3, use=NULL
                        ){
  tData<-as.matrix(tData)
  nv<-(ncol(tData)-1)/embedD
  model$data<-mxData(tData, type="raw")
  selCol<-round(seq(1, nv*embedD, embedD))
  
  
  # Threshold-based data exclusion func -------------------------------------
  cutData<-function(thres){
    

    dat<-tData
    pointInfo.s.exc<-which(pointInfo.s>thres, arr.ind=T)
    # print(pointInfo.s.exc)
    
    if(length(pointInfo.s.exc)>0){
      pointInfo.s.exc[,2]<-selCol[pointInfo.s.exc[,2]]
      for(i in 0:(embedD-1)){
        pointInfo.s.exc.i<-pointInfo.s.exc
        pointInfo.s.exc.i[,2]<-pointInfo.s.exc[,2]+i
        dat[pointInfo.s.exc.i]<-NA
      }
      exclusions<-list()
      for(i in 1:nv){
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
  plotEA<-function(num=NULL, fitODE=F){
    
    par(mfrow=c(nVars+1,1), mai=c(0.3, 1, 0.025, 0.5))
    
    plotColors<-rainbow(nv, s=1, v=1, alpha=1)
    plotColors.muted<-rainbow(nv, s=0.5, v=0.25, alpha=1)
    varCols<-round(seq(1, nv*embedD, embedD))
    
    #Variable graph
    for(i in 1:nv){
      plot(tData[,varCols[i]], type="l", lwd=1, col=plotColors.muted[i], xaxt='n', xlab=NULL, ylab=paste0("Series X",i))
      abline(v=out$exclusions[[i]], col="orange", lty=3)
      if(fitODE){
        lines(1:nrow(tData), tOscData.ex[,i], type="l", lwd=1, col=plotColors[i])
        abline(means.ex[i,1], means.ex[i,2], col=plotColors[i], lty=1, lwd=1)
      
        # scores<-mxFactorScores(out$fit, type='WeightedML')
        # scores.comb<-as.matrix(as.data.frame(scores)[,1:3])%*%t(LMatrix)
        # lines(scores.comb[,3], type="l", col="blue")
        
        
      }
    }
    # rm(means.ex)
    # rm(tOscData.ex)
    #-2LL graph
    plot(c(0,nrow(tData)), c(0,0), type="n", xlab=NULL, ylab="Marginal -2LL", ylim=c(min(pointInfo.s, na.rm=T), max(pointInfo.s, na.rm=T)))
    for(i in 1:nv){
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
    pointInfo<-matrix(NA,nrow(tData),nv)
    
    for(h in 1:nv){
      for(i in 1:nrow(tData)){
        cat(paste0("\rCalculating likelihood for timepoint ",i," of X",h,"..."))
        tData.drop<-tData
        for(j in 0:(embedD-1)){
        tData.drop[i,selCol[h]+j]<-NA
        }
        model$data<-mxData(tData.drop, type="raw")
        model<-mxRun(model, silent=T)
        if(model$output$status$code >= 6 && tries > 0){
          model<-mxTryHard(model, extraTries = tries, greenOK=T)
        }
        pointInfo[i,h] <- -1 * model$output$Minus2LogLikelihood
      }
      pointInfo[,h]<- pointInfo[,h] - min(pointInfo[,h],na.rm=T)
      pointInfo[,h]<- pointInfo[,h] / sum(pointInfo[,h],na.rm=T)
    }
    pointInfo.s<-as.matrix(pointInfo)
    #Combine event vectors like old version?
    # pointInfo.s<-matrix(rowSums(pointInfo.s), nrow(tData), nv) 
    
  }else{
    pointInfo.s<-use
  }
  
  
  
  
  
  # Data cleaning and analysis ----------------------------------------------
  # fileid<<-0
  
  logErr<-norms<-NULL
  
  if(optim){
    cat("\nFinding exclusion threshold...\r")
    
    info.sorted<-sort(as.vector(pointInfo.s), decreasing=T)
    cuts.norm<<-NULL
    cutProp<-length(info.sorted)- min.data
    
    for(i in 1:cutProp){
      out<-cutData(info.sorted[i])
      
      i.MAE.int<-fitODE(tData, out)
      i.UX<-mean(omxGetParameters(out$fit)[grepl("Ux",names(omxGetParameters(out$fit)))], na.rm=T)
      # i.UX<-out$fit$output$estimate[["Ux1"]]
      # i.zeta.std<-out$fit$output$standardErrors["zetaX1",]
      # i.zeta.std<-sum(out$fit$output$standardErrors[grepl("zetaX",rownames(out$fit$output$standardErrors))])
      i.par.std<-mean(out$fit$output$standardErrors[grepl("zetaX",rownames(out$fit$output$standardErrors))])
      # i.m2ll<-out$fit$output$Minus2LogLikelihood
      i.AIC<-summary(out$fit)$AIC
      i.BIC<-summary(out$fit)$BIC
      
      percentExclude<-length(unique(unlist(out$exclusions))) / nrow(tData)
      norms.row<-as.vector(c(i.MAE.int, i.AIC, i.BIC, i.UX, i.par.std, percentExclude))
      names(norms.row)<-c("MAE","AIC","BIC","UX","SEzeta","fEvents")
                        
      
      # normTest<-glm(normTest.dat[,1]~normTest.dat[,-1])
      this.norm<- i.AIC#norms.row%*%(-1/as.matrix(normTest$coefficients[-1]))
      
      logErr<- -log((omxGetParameters(out$fit)[["zetaX1"]]-indDyn[1,2])^2 + (omxGetParameters(out$fit)[["etaX1"]]-indDyn[1,1])^2)
      norms<-rbind(norms, c(logErr,norms.row))

      if(plotting){
        # fileid<<-fileid+1
        # png(paste0("gfx/",fileid,".png"), width=800, height = 800)
          plotEA(num = info.sorted[i], fitODE=T)
        # dev.off()
        }
      
      cuts.norm<<-rbind(cuts.norm, this.norm)
    }
    i.opt<-which(cuts.norm==min(cuts.norm, na.rm=T)[1])
    out<-cutData(info.sorted[i.opt])
    i.thres<-info.sorted[i.opt]

  }else{
    out<-cutData(threshold.sv)
    i.thres<-threshold.sv
    # best.norm<-summary(out$fit)$AIC
  }
  
  # Final plot
  if(plotting){
    # fileid<<-fileid+1
    # png(paste0("gfx/",fileid,".png"), width=800, height = 800)
    
    fitODE(tData, out)
    plotEA(num=i.thres, fitODE = T)
    
    # dev.off()
  }

  # Output ------------------------------------------------------------------
  
  normTest.dat<<-rbind(normTest.dat, norms)
  
  return(list("fit"=out$fit,"exclusions"=out$exclusions, "info"=pointInfo.s, "norms"=norms))
}














