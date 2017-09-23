# ---------------------------------------------------------------------
#  Use estimates from an LDE model to generate a purelydeterministic prediction
#  Author: Kevin McKee
# ---------------------------------------------------------------------

# Define the damped linear oscillator function.
fitODE<-function(data, model, plot=T){
  
  # data<-tDat
  # model<-out
  # dx<-initCond
  # nVars<-ncol(data)-1/embedD
  nVars<-as.numeric(model$fit$nVars$values)
  embedD<-as.numeric(model$fit$embedD$values)
  
  varNames<-as.vector(rbind(paste("x",1:nVars, sep=""), paste("dx",1:nVars, sep="")))
  
  tCenter<-round(embedD/2)
  selCol<-round(seq(1, nVars*embedD, embedD))
  theTimes  <- seq(tCenter, nrow(data), 1)

  #Get parameters from the model--------------------------------------------
  initCond<-indDyn<-means<-matrix(0, nVars, 2); gammas<-diag(1, nVars, nVars)
  means<-matrix(omxGetParameters(model$fit)[cbind(paste0('intX',1:nVars), paste0('slopeX',1:nVars))], nVars, 2)
  indDyn<-matrix(omxGetParameters(model$fit)[p.l], nVars, 2)
  gammas.L<-gammas.S<-diag(1, nVars, nVars)
  gammaLabs.L<-omxGetParameters(model$fit)[gL.l]
  gammaLabs.L<-gammaLabs.L[complete.cases(gammaLabs.L)]
  gammaLabs.S<-omxGetParameters(model$fit)[gS.l]
  gammaLabs.S<-gammaLabs.S[complete.cases(gammaLabs.S)]
  gammas.L[which(gammas.L==0)]<-gammaLabs.L
  gammas.S[which(gammas.S==0)]<-gammaLabs.S
  parms<-list(indDyn, gammas.L, gammas.S)
  
  
  
  # ODE ---------------------------------------------------------------------
  DLOmodel <- function(t, prevState, parms) {
    dx<-matrix(prevState, nVars, 2, byrow=T)
    states<-matrix(NA, nVars, 2)
    with(list(parms, states, dx),
         {
           states[,1] <- dx[,2]
           states[,2] <- parms[[2]]%*%(dx*parms[[1]])[,1] + parms[[3]]%*%(dx*parms[[1]])[,2]
           vals<-as.vector(t(states))
           names(vals)<-varNames
           res<-vals
           list(res)
         }
    )
  }
  

  
  # Set up event data -------------------------------------------------------
  
  varNames.vec<-NULL
  has.events<-NULL
  eventTimes.bounds<-c(1:(tCenter), (nrow(data)-(embedD)):nrow(data))
  eventTimes<-list()
  
  #Step 1... set up initial conditions
    for(i in 1:nVars){
      eventTimes[i]<-list(model$exclusions[[i]][!model$exclusions[[i]]%in%eventTimes.bounds])
      # print(eventTimes[[i]])
      
      eventTimes[i]<-list(unique(c(eventTimes[[i]], eventTimes[[i]]+1)))
      # print(eventTimes[[i]])
      
      if(!is.null(eventTimes[[i]])){
        has.events<-c(has.events, i)
      }
    } 
  
  # print(eventTimes)
  # print("A")
  
  
    shocks.m<-NULL
    for(i in 1:nVars){
      eNoise<-omxGetParameters(model$fit)[[paste0("Ux",i)]] 
      eSig<-omxGetParameters(model$fit)[[paste0("VX",i)]]
      eNoiseVar<-eNoise/(eSig+eNoise)
      
      #Initial Conditions
      #Option 1: Smoothed level and slope
      ic.level.s<-data[1, selCol[i]:(selCol[i]+(embedD-1)) ]%*%(L1/embedD) - means[i,1] - means[i,2]*tCenter
      ic.slope.s<-data[1, selCol[i]:(selCol[i]+(embedD-1)) ]%*%(L2/(embedD-1)) - means[i,2]*tCenter
      #Option 2: Raw level and slope
      ic.level.r<-data[tCenter, selCol[i]] - means[i,1] - means[i,2]*tCenter
      ic.slope.r<-data[tCenter, selCol[i]+1]-data[tCenter, selCol[i]] - means[i,2]*tCenter
      #Variable ratio between the two depending on estimated noise variance
      initCond[i,1]<-(1-eNoiseVar)*ic.level.r + eNoiseVar*ic.level.s
      initCond[i,2]<-(1-eNoiseVar)*ic.slope.r + eNoiseVar*ic.slope.s
      
      
      if(i%in%has.events){
      varNames.vec<-c(varNames.vec, rep(varNames[(i*2-1)],length(theTimes)), rep(varNames[i*2],length(theTimes)))
        
      #Shocks
      shocks.slopes.r<-shocks.slopes.s<-shocks.levels.r<-shocks.levels.s<-rep(NA, nrow(data))
     
       #Option 1: Smoothed level and slope
      shocks.levels.s[eventTimes[[i]]]<-data[eventTimes[[i]]-tCenter+1, selCol[i]:(selCol[i]+(embedD-1)) ]%*%(L1/embedD)  - means[i,1] - means[i,2]*eventTimes[[i]]
      X.s<-cbind(rep(1,embedD), 0:(embedD-1))
      Y.s<-t(as.matrix(data[eventTimes[[i]], selCol[i]:(selCol[i]+(embedD-1))]))
      XYslope<-(solve(t(X.s)%*%X.s)%*%t(X.s)%*%Y.s)
      shocks.slopes.s[eventTimes[[i]]]<-XYslope[2,]-means[i,2]

      #Option 2: Raw level and slope
      shocks.levels.r[eventTimes[[i]]]<-data[eventTimes[[i]], selCol[i]] - means[i,1] - means[i,2]*eventTimes[[i]]
      shocks.slopes.r[eventTimes[[i]]]<-data[eventTimes[[i]], selCol[i]+1]-data[eventTimes[[i]], selCol[i]] - means[i,2]
      
      #Variable ratio between the two depending on estimated noise variance
      shocks.levels<-(1-eNoiseVar)*shocks.levels.r + eNoiseVar*shocks.levels.s
      shocks.slopes<-(1-eNoiseVar)*shocks.slopes.r + eNoiseVar*shocks.slopes.s
      
      shocks.levels<-shocks.levels[theTimes]      
      shocks.slopes<-shocks.slopes[theTimes] 
      shocks.m<-c(shocks.m, shocks.levels, shocks.slopes)
      }
    }
    
    if(length(has.events)>0){
    eventdat <- data.frame(var = varNames.vec,
                           time = rep(theTimes, length(has.events)*2),
                           value = shocks.m,
                           method = rep("replace", length(theTimes)*length(has.events)*2))
    
    eventdat<-eventdat[which(!is.na(eventdat$value)),]
    # print(eventdat)
    }else{
      eventdat<-NULL
    }
    
    # print(eventdat)
  # Generate data -----------------------------------------------------------
  xstart <- as.vector(t(initCond))
  names(xstart)<-varNames
  pred<-NULL
  
  eventTimes.collapse<-unique(unlist(eventTimes))
  
  # print(eventTimes.collapse)
  # print(varNames)
  # print(theTimes)
  
  if(length(eventTimes.collapse)>0){
    theTimes.micro<- c(eventTimes.collapse+0.000001)
    theTimes.event<-sort(c(theTimes, theTimes.micro))
    
    # print(theTimes.event)
    # print(eventdat)
    
    try(pred <- as.data.frame(ode(xstart, theTimes.event, DLOmodel, parms, 
                                  events = list(data=eventdat), method="lsoda")))
    
    
    pred[which(pred$time%in%eventTimes.collapse),-1]<-pred[which(pred$time%in%theTimes.micro),-1]
    pred<-pred[-which(pred$time%in%theTimes.micro),]
  }else{
    try(pred <- as.data.frame(ode(xstart, theTimes, DLOmodel, parms,
                                  events = list(data=eventdat), method="lsoda")))
  }
  # print(pred)
  
  if(is.null(pred)){return(NA)}
  
  # print(pred)
  # Outputs -----------------------------------------------------------------
  means.adj<- matrix(c(rep(1, length(theTimes)), 1:length(theTimes)),length(theTimes), 2)%*%t(means)
  tOscData<-rbind(matrix(NA,tCenter-1, nVars), as.matrix(pred[,(1:nVars)*2]) + means.adj)
  
  # tOscData.ex<<-tOscData
  # means.ex<<-means
  
  # print(dim(tOscData)); print(dim(data))
  residuals<-as.matrix(tOscData-data[,selCol])
  MSE<-mean(residuals^2, na.rm=T)
  MAE<-mean(abs(residuals),na.rm=T)
  
  return(list("MAE"=MAE,"MSE"=MSE,"Means"=means,"Prediction"=tOscData))
  
}



