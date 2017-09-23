# ---------------------------------------------------------------------
#  Simulate a single Damped Linear Oscillator
#  Author: Kevin McKee
# ---------------------------------------------------------------------
N_genData<-function(N=100, SNR=2, numVars, parameters, initialConditions,
                    events=0, sigmaEvent=NULL, eventNormal=F, eventScale=1, eventType="level",
                    plotting=F) {
  
  
  
  theTimes  <- 1:N
  varNames<-as.vector(rbind(paste("x",1:numVars, sep=""), paste("dx",1:numVars, sep="")))
  
  
  # Format inputs -----------------------------------------------------------
  indDyn<-cbind(
    matrix(parameters[substr(names(parameters),1,3)=="eta"], ncol=1),
    matrix(parameters[substr(names(parameters),1,4)=="zeta"],  ncol=1))
  
  randMeans<-cbind(
    matrix(parameters[substr(names(parameters),1,3)=="int"], ncol=1),
    matrix(parameters[substr(names(parameters),1,5)=="slope"], ncol=1))
  
  numVars.tot<-nrow(indDyn)
  
  gammas.L<-matrix(parameters[substr(names(parameters),1,5)=="gamma" & substr(names(parameters),9,9)=="L"], numVars.tot, numVars.tot)
  gammas.S<-matrix(parameters[substr(names(parameters),1,5)=="gamma" & substr(names(parameters),9,9)=="S"], numVars.tot, numVars.tot)
  
  
  
  #Trim down to numvars
  indDyn<-matrix(indDyn[1:numVars,],  nrow=numVars, ncol=2)
  randMeans<- matrix(randMeans[1:numVars,], nrow=numVars, ncol=2)
  
  gammas.L<-gammas.L[1:numVars, 1:numVars]
  gammas.S<-gammas.S[1:numVars, 1:numVars]
  xstart <- initialConditions[1:(numVars*2)]
  
  
  # print(indDyn)
  # print(gammas.L)
  # print(gammas.S)
  # print(randMeans)
  
  # Function for lsoda to generate data from a system of diff eqs. ---------
  DLOmodel <- function(t, prevState, parms) {
    dx<-matrix(prevState, numVars, 2, byrow=T)
    states<-matrix(NA, numVars, 2)
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
  
  
  # Generate events ---------------------------------------------------------
  
  #If a probability is entered for events, generate events.
  if(is.numeric(events)){
    
    #generate shocks  
    shocks<-list()
    for(i in 1:numVars){
      shocks[[i]]<-eventScale*rbinom(N,1,events/numVars)*sample(c(-1, 1), N, replace=T)
      if(eventNormal){shocks[[i]]<-shocks[[i]]*rnorm(N)}
    }
    
    shocks.m<-NULL
    for(i in 1:numVars){
      if(eventType=="slope"){
        shocks.m<-c(shocks.m, rep(0,N), shocks[[i]])
      }else if(eventType=="both"){
        shocks.m<-c(shocks.m, rep(shocks[[i]],2))
      }else{
        shocks.m<-c(shocks.m, shocks[[i]], rep(0,N))
      }
    }
    varNames.vec<-NULL
    for(i in varNames){varNames.vec<-c(varNames.vec, rep(i, N))}
    
    eventdat <- data.frame(var = varNames.vec,
                           time = rep(theTimes, numVars*2),
                           value = shocks.m,
                           method = rep("add", N*numVars*2))
    eventdat<-eventdat[eventdat$value!=0,]
  }
  
  # print(eventdat)
  
  #If a list of event times is entered, generate with those events. 
  ###to do
  
  
  
  # Generate data -----------------------------------------------------------
  parms<-list(indDyn, gammas.L, gammas.S)
  names(xstart)<-varNames
  # sim <- as.data.frame(lsoda(xstart, theTimes, DLOmodel, parms, events = list(data=eventdat)))
  
  
  
  eventTimes.collapse<-unique(eventdat$times)
  
  if(length(eventTimes.collapse)>0){
    theTimes.micro<- c(eventTimes.collapse+0.000001)
    theTimes.event<-sort(c(theTimes, theTimes.micro))
    
    try(sim <- as.data.frame(ode(xstart, theTimes.event, DLOmodel, parms,
                                 events = list(data=eventdat), method="lsoda")))
    
    
    sim[which(sim$time%in%eventTimes.collapse),-1]<-sim[which(sim$time%in%theTimes.micro),-1]
    sim<-sim[-which(sim$time%in%theTimes.micro),]
  }else{
    try(sim <- as.data.frame(ode(xstart, theTimes, DLOmodel, parms,
                                 events = list(data=eventdat), method="lsoda")))
  }
  
  
  
  # ----------------------------------
  # Scale error for a chosen signal to noise ratio.
  tESD <- 1 / SNR
  tOscData<-NULL
  # tEmbedded<-NULL
  
  
  for(i in (1:numVars)*2){
    tSD <- sd(sim[,i])
    thisVar<-sim[,i]/tSD + rnorm(N, mean=0, sd=tESD)
    thisVar<- thisVar+ randMeans[i/2,2]*(1:length(thisVar)) + randMeans[i/2,1]
    tOscData <- cbind(tOscData, thisVar)
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  dimnames(tOscData) <- list(NULL, paste("X",1:numVars,sep=""))
  
  
  
  
  if(plotting){
    plotColors<-rainbow(numVars, v=0.75)
    plotColors.muted<-rainbow(numVars, s=1, v=0.75, alpha=0.5)
    plot(c(min(theTimes), max(theTimes[1:N])), c(-5, 5),
         xlab="Time",
         ylab="Score",
         type='n',
    main=paste(numVars,"-Oscillator Network", sep=""))
    for(i in 1:numVars){
      lines(sim$time, tOscData[,i], type='l', lwd=1, col=plotColors.muted[i])
    }
    for(i in 1:numVars){
      lines(shocks[[i]], type="h", lwd=3, col=plotColors[i])
    }
  }
  
  return(tOscData)
}







