
MSD_simPar <- function(parName, pars, trueVals=NULL, useEventAnalysis=F, defaults=list() ) {
  out.par.e<- 
    out.par.se <-
    out.par.sd <-
    out.par.e.m <-
    out.par.e.med <-
    out.par.se.m <-
    out.par.e.sd <- NULL
  
  
  for(h in pars) {
    defaults[[parName]]<-h
    
    source("lib/MS/MSD_matrices.R", local=T)
    source("lib/MS/MS_model.R", local=T)
    source("lib/MS/MS_constraints.R", local=T)
    
    # Main simulation code ----------------------------------------------------
    Ests<-Se<-NULL
    for(i in 1:defaults$nTrials){
      
      #Generate and embed data
      dat<-MSD_genData(defaults$N, SNR=defaults$SNR, trueVals[i,], initCond[i,], defaults,
                     eventScale = 1, plotting=T)
      
      

      
      tEmbedded<- gllaEmbed(dat, embed=defaults$embedD1*defaults$embedD2, tau=defaults$theTau, idColumn=FALSE)
      
      
      tEmbedded<-NULL
      disConts<-seq(1, defaults$N,defaults$embedD1)
      
      for(i in 1:defaults$nVars){
        dat.emb <- gllaEmbed(dat[,i], embed=defaults$embedD, tau=defaults$theTau, idColumn=FALSE)
        dat.emb<-dat.emb[disConts[disConts<nrow(dat.emb)],]
        tEmbedded <- cbind(tEmbedded, dat.emb)
      }
      
      colnames(tEmbedded)<-manifestVars
      tEmbedded<-cbind(tEmbedded, "Occasion"=1:nrow(tEmbedded))
      
      
      #Analysis with or without events
      if(useEventAnalysis){
        EA<-eventAnalysis(tEmbedded, model, plotting=T, prediction=F, tries=defaults$tryHard.tries)
        model<-EA$fit
      }else{
        model$data<-mxData(tEmbedded, type="raw")
        model<-mxRun(model, silent=F)
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

    
    #Store total results
    out.par.e <- rbind(out.par.e, Ests)
    out.par.se <- rbind(out.par.se, Se)
    
  }
  
  return( list("Estimates"=out.par.e, 
               "Std.Errs"=out.par.se, 
               "trueVals"=trueVals))
}