
MS_simPar <- function(parName, pars, trueVals) {
    out.par.e<- 
    out.par.se <-
    out.par.sd <-
    out.par.e.m <-
    out.par.e.med <-
    out.par.se.m <-
    out.par.e.sd <- NULL
    
    
    
  for (h in pars) {
    eval(parse(text = paste(parName, "<<-h", sep = "")))
    source("lib/MS/MS_matrices.R")
    source("lib/MS/MS_model.r")
    source("lib/MS/MS_constraints.R")
    # Main simulation code ----------------------------------------------------
    Ests<-Se<-NULL
    for(i in 1:nTrials){
      cat(paste("\r",i))
        etaX1<<-trueVals[i,'eta']
        zetaX1<<-trueVals[i,'zeta']
        etaX2<<-trueVals[i, 'etaM']
        zetaX2<<-trueVals[i, 'zetaM']
        gammaX1<<- -etaX1

        source("lib/MS/MS_genData.old.R")
        
        
      model$data<-mxData(tEmbedded, type="raw")
      model<-mxTryHard(model, extraTries = 20, greenOK=T)
      
      Ests<-rbind(Ests, cbind(h, t(omxGetParameters(model))))
      Se<-rbind(Se, cbind(h,t(model$output$standardErrors)))
      }
    
    Ests.m<-t(colMeans(Ests, na.rm=T))#; row.names(Ests.m)<-NULL
    Se.m<-t(colMeans(Se, na.rm=T))#; row.names(Se.m)<-NULL
    Ests.sd<-cbind(h, t(apply(Ests[,-1],2,sd, na.rm=T)))
    Ests.med<-cbind(h, t(apply(Ests[,-1],2,median, na.rm=T)))
    
    out.par.e <- rbind(out.par.e, Ests)
    out.par.e.m <- rbind(out.par.e.m, Ests.m)
    out.par.e.med <- rbind(out.par.e.med, Ests.med)
    out.par.e.sd <- rbind(out.par.e.sd, Ests.sd)
    out.par.se <- rbind(out.par.se, Se)
    out.par.se.m <- rbind(out.par.se.m, Se.m)
  }
    
    out.par.e<-as.data.frame(out.par.e); colnames(out.par.e)[1]<-parName
    out.par.e.m <- as.data.frame(out.par.e.m); colnames(out.par.e.m)[1]<-parName
    out.par.e.med <- as.data.frame(out.par.e.med); colnames(out.par.e.med)[1]<-parName
    out.par.e.sd<-as.data.frame(out.par.e.sd); colnames(out.par.e.sd)[1]<-parName
    out.par.se<-as.data.frame(out.par.se); colnames(out.par.se)[1]<-parName
    out.par.se.m<-as.data.frame(out.par.se.m); colnames(out.par.se.m)[1]<-parName
 
    return( list("Estimates"=out.par.e, 
                 "Std.Errs"=out.par.se, 
                 "Mean.Estimates"=out.par.e.m, 
                 "Median.Estimates"=out.par.e.med, 
                 "Mean.Std.Errs"=out.par.se.m,
                 "SD.Estimates"=out.par.e.sd,
                 "trueVals"=trueVals))
}