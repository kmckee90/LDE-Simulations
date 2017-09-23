
N_simPar <- function(parName, pars, trueVals=NULL, useEventAnalysis=F, defaults) {
    out.par.e<- 
      out.par.se <-
      out.par.sd <-
      out.par.e.m <-
      out.par.e.med <-
      out.par.se.m <-
      out.par.e.sd <- NULL
    
    
    for(h in pars) {
      defaults[[parName]]<-h
      
      source("lib/NL/N_matrices.R", local=T)
      source("lib/NL/N_model.R", local=T)
      source("lib/NL/N_constraints.R", local=T)
      
      print(model$L$values)
      # Main simulation code ----------------------------------------------------
      Ests<-Se<-NULL
      for(i in 1:defaults$nTrials){
        cat(paste("\r",h,i))
        #Generate and embed data
        # print(defaults$embedD)
        dat<-N_genData(defaults$N, SNR=defaults$SNR, numVars=defaults$nVars, trueVals[i,], initCond[i,], 
                       events=defaults$pEvent, eventNormal=defaults$eventNormal, eventType="level", eventScale = 1,
                       plotting=defaults$plotting)
        
        tEmbedded<-NULL
        for(i in 1:defaults$nVars){
          tEmbedded <- cbind(tEmbedded, gllaEmbed(dat[,i], embed=defaults$embedD, tau=defaults$theTau, idColumn=FALSE))
        }
        # print(head(tEmbedded))
        colnames(tEmbedded)<-manifestVars
        tEmbedded<-cbind(tEmbedded, "Occasion"=1:nrow(tEmbedded))
        
        
        #Analysis with or without events
        if(useEventAnalysis){
          EA<-eventAnalysis(tEmbedded, model, plotting=defaults$plotting, prediction=F, tries=defaults$tryHard.tries)
          model<-EA$fit
        }else{
          model$data<-mxData(tEmbedded, type="raw")
          model<-mxRun(model, silent=T)
        if(model$output$status$code%in%c(5,6)){
          model<-mxTryHard(model, extraTries = defaults$tryHard.tries, greenOK=T)
        }
      }
      
      #Store results
      Ests<-rbind(Ests, c("par"=h, omxGetParameters(model), "statuscode"=model$output$status$code))
      Se<-rbind(Se, c("par"=h, model$output$standardErrors, "statuscode"=model$output$status$code))
      
    }
    
    colnames(Se)<-colnames(Ests)
    colnames(Se)[1]<-colnames(Ests)[1]<-parName
    
    #calculate means and medians 
    # Ests.m<-colMeans(Ests, na.rm=T)
    # Se.m<-colMeans(Se, na.rm=T)
    # Ests.sd<-c(h, apply(Ests[,-1],2,sd, na.rm=T))
    # Ests.med<-c(h, apply(Ests[,-1],2,median, na.rm=T))
    
    
    #Store total results
    out.par.e <- rbind(out.par.e, Ests)
    out.par.se <- rbind(out.par.se, Se)
    
    # rm(defaultEnv)
    
  }
  
  return( list("Estimates"=out.par.e, 
               "Std.Errs"=out.par.se, 
               "trueVals"=trueVals))
}