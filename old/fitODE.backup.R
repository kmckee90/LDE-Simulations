# ---------------------------------------------------------------------
#  Simulate a single Damped Linear Oscillator
#  Author: Original script by Steve Boker, modified by Kevin McKee
# ---------------------------------------------------------------------
# ----------------------------------
# Define the damped linear oscillator function.
fitODE<-function(data, model, plot=T){
  
  # data<-tDat
  # model<-out
  # dx<-initCond
  theTimes  <- 1:nrow(data)
  varNames<-as.vector(rbind(paste("x",1:nVars, sep=""), paste("dx",1:nVars, sep="")))
  varNames.vec<-NULL
  for(i in varNames){
    varNames.vec<-c(varNames.vec, rep(i, nrow(data)))
  }
  selCol<-round(seq(1, nVars*embedD, embedD))
  
  # Initial conditions and parameters--------------------------------------------
  initCond<-indDyn<-means<-matrix(0, nVars, 2); gammas<-diag(1, nVars, nVars)
  means<-matrix(omxGetParameters(model$fit)[cbind(paste0('intX',1:nVars), paste0('slopeX',1:nVars))], nVars, 2)
  
  for(i in 1:nVars){
    initCond[i,1]<-data[1, selCol[i]] - means[i,1]
    initCond[i,2]<-data[2, selCol[i]]-data[1, selCol[i]] - means[i,2]
  }
  
  indDyn<-matrix(omxGetParameters(model$fit)[p.l], nVars, 2)
  
  gammas.L<-gammas.S<-diag(1, nVars, nVars)
  gammaLabs.L<-omxGetParameters(model$fit)[gL.l]
  gammaLabs.L<-gammaLabs.L[complete.cases(gammaLabs.L)]
  gammaLabs.S<-omxGetParameters(model$fit)[gS.l]
  gammaLabs.S<-gammaLabs.S[complete.cases(gammaLabs.S)]
  gammas.L[which(gammas.L==0)]<-gammaLabs.L
  gammas.S[which(gammas.S==0)]<-gammaLabs.S
  
  parms<-list(indDyn, gammas.L, gammas.S, means)
  
  
  
  
  
  
  
  
  
  DLOmodel <- function(t, prevState, parms) {
    dx<-matrix(prevState, nVars, 2, byrow=T)
    states<-matrix(NA, nVars, 2)
    # meanMat<- matrix(c(1,t),2,1)
    with(list(parms, states, dx),
         {
           states[,1] <- dx[,2]
           states[,2] <- parms[[2]]%*%(dx*parms[[1]])[,1] + parms[[3]]%*%(dx*parms[[1]])[,2] #+ parms[[4]]%*%meanMat
           vals<-as.vector(t(states))
           varNames<-as.vector(rbind(paste("x",1:nVars, sep=""), paste("dx",1:nVars, sep="")))
           names(vals)<-varNames
           res<-vals
           list(res)
         }
    )
  }
  
  

  
  
  # Set up event data -------------------------------------------------------
  if(length(model$exclusions)>0 && !is.null(model$exclusions) && is.numeric(model$exclusions)){
    shocks.m<-NULL
    for(i in 1:nVars){
      shocks.levels<-shocks.slopes<-rep(NA, nrow(data))
      shocks.levels[model$exclusions]<-data[model$exclusions, selCol[i] ]
      shocks.slopes[model$exclusions]<-data[model$exclusions, selCol[i]+1 ]-data[model$exclusions, selCol[i] ]
      # shocks.slopes<-shocks.slopes[1:length(shocks.levels)]
      shocks.m<-c(shocks.m, shocks.levels, shocks.slopes)
    }
    
    eventdat <- data.frame(var = varNames.vec,
                           time = rep(theTimes, nVars*2),
                           value = shocks.m,
                           method = rep("replace", nrow(data)*nVars*2))
    
    # eventdat<-eventdat[eventdat$var=="x1",]
    
    eventdat$value[eventdat$var=="x1"]<-eventdat$value[eventdat$var=="x1"] - means[i,1] - means[i,2]*eventdat$time[eventdat$var=="x1"]
    eventdat$value[eventdat$var=="dx1"]<-eventdat$value[eventdat$var=="dx1"] -  means[i,2] 
    # 
    eventdat<-eventdat[which(!is.na(eventdat$value)),]
    # eventdat<-eventdat[which(eventdat$var=="x1"),]
    
    # eventdat$time<-eventdat$time+1
    print(eventdat)
  }else{
    eventdat<-NULL
  }
  

  
  # Generate data -----------------------------------------------------------
  tOffsets <- c(1:nrow(data))
  xstart <- as.vector(t(initCond))
  names(xstart)<-varNames
  pred<-NULL
  
  if(length(model$exclusions)>0){
    theTimes.micro<- c(model$exclusions+0.000001)
    theTimes.event<-sort(c(theTimes, theTimes.micro))

  # try(pred <- as.data.frame(ode(xstart, theTimes, DLOmodel, parms, 
  #                                 events = list(data=eventdat), method="euler"
  #                                 )))
  
  try(pred <- as.data.frame(ode(xstart, theTimes.event, DLOmodel, parms, 
                                events = list(data=eventdat), method="lsoda")))
      pred[which(pred$time%in%model$exclusions),-1]<-pred[which(pred$time%in%theTimes.micro),-1]
      pred<-pred[-which(pred$time%in%theTimes.micro),]

  }else{
    try(pred <- as.data.frame(ode(xstart, theTimes, DLOmodel, parms, 
                                  events = list(data=eventdat), method="lsoda")))
  }
  
  if(is.null(pred)){
    return(NA)
  }

  # Outputs -----------------------------------------------------------------
  means.adj<- matrix(c(rep(1, nrow(data)), 1:nrow(data)), nrow(data), 2)%*%t(means)
  # means.adj[model$exclusions,]<-0
  # means.adj[1,]<-0
  
  
  tOscData<-as.matrix(pred[,(1:nVars)*2]) + means.adj

  tOscData.ex<<-tOscData
  means.ex<<-means
  residuals<-as.matrix(tOscData[-model$exclusions,]-data[-model$exclusions, selCol])
  MSE<-sum(residuals^2, na.rm=T) / nrow(residuals)
  MSE<- sqrt(MSE)
  
  # print("RAN")
  
  return(MSE)
  
}

