# ---------------------------------------------------------------------
#  Simulate a single Damped Linear Oscillator
#  Author: Kevin McKee
# ---------------------------------------------------------------------
MSD_genData<-function(N=100, SNR=2, parameters, initialConditions, defaults,
                     eventScale=1, plotting=F) {
  
  # nVars<-1
  
  theTimes  <- 1:N
  varNames<-as.vector(rbind(paste("x",1:defaults$nVars, sep=""), paste("dx",1:defaults$nVars, sep=""), paste("xM",1:defaults$nVars, sep=""), paste("dxM",1:defaults$nVars, sep="")))
  # varNames<-c(varNames.1, varNames.2)
  
  # Format inputs -----------------------------------------------------------
  indDyn.short<-cbind(
    matrix(parameters[substr(names(parameters),1,4)=="etaX"], ncol=1),
    matrix(parameters[substr(names(parameters),1,5)=="zetaX"],  ncol=1))
  
  indDyn.long<-cbind(
    matrix(parameters[substr(names(parameters),1,5)=="mEtaX"], ncol=1),
    matrix(parameters[substr(names(parameters),1,6)=="mZetaX"],  ncol=1))
  
  randMeans<-cbind(
    matrix(parameters[substr(names(parameters),1,3)=="int"], ncol=1),
    matrix(parameters[substr(names(parameters),1,5)=="slope"], ncol=1))
  
  nVars.tot<-nrow(indDyn.short)
  
  
  gammas.L<-matrix(parameters[substr(names(parameters),1,5)=="gamma" & substr(names(parameters),9,9)=="L"], nVars.tot, nVars.tot)
  gammas.S<-matrix(parameters[substr(names(parameters),1,5)=="gamma" & substr(names(parameters),9,9)=="S"], nVars.tot, nVars.tot)
  
  print(indDyn.short)
  print(indDyn.long)
  print(randMeans)
  print(gammas.L)
  print(gammas.S)
    
  #Trim down to nVars
  indDyn.short<-matrix(indDyn.short,  nrow=defaults$nVars, ncol=2)
  indDyn.long<-matrix(indDyn.long,  nrow=defaults$nVars, ncol=2)
  randMeans<- matrix(randMeans, nrow=defaults$nVars, ncol=2)
  gammas.L<-gammas.L[1:defaults$nVars, 1:defaults$nVars]
  gammas.S<-gammas.S[1:defaults$nVars, 1:defaults$nVars]
  
  xstart <- initialConditions[1:(defaults$nVars*4)]
  
  
  
  print(indDyn.short)
  print(indDyn.long)
  
  DLOmodel <- function(t, prevState, parms) {
    dx<-matrix(prevState, defaults$nVars, 4, byrow=T)
    states<-matrix(NA, defaults$nVars, 4)
    with(list(parms, states, dx),
         {
           states[,1] <- dx[,2]
           states[,2] <- parms[[2]]%*%(dx[,1:2]*parms[[1]])[,1] + parms[[3]]%*%(dx[,1:2]*parms[[1]])[,2]
           
           states[,3] <- dx[,4]
           states[,4] <- parms[[2]]%*%(dx[,3:4]*parms[[4]])[,1] + parms[[3]]%*%(dx[,3:4]*parms[[4]])[,2]
           
           vals<-as.vector(t(states))
           names(vals)<-varNames
           res<-vals
           list(res)
         }
    )
  }
  
  # print(eventdat)
  theTimes  <- 1:N  
  
  disConts<-seq(1, N,defaults$embedD1)
  nEvent<-length(disConts) #round(runif(1, 1, length(disConts)*defaults$pEvent-2))
  # shocks.t<-sort(sample(disConts[-1], nEvent, replace=F))
  shocks.v<-rnorm(nEvent,0,defaults$eventScale)
  
  eventdat <- data.frame(var = rep(varNames[1], nEvent),
                         time = disConts-1,
                         value = shocks.v,
                         method = rep("replace", nEvent))
  
  #If a list of event times is entered, generate with those events. 
  ###to do
  
  
  
  # Generate data -----------------------------------------------------------
  parms<-list(indDyn.short, gammas.L, gammas.S, indDyn.long)
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
  
  
  short<-sim[,"x1"]
  long<-sim[,"xM1"]

  segs<-round(length(short)/defaults$embedD1)
  short.seg<-sort(rep(1:segs, defaults$embedD1))
  short.seg.rand<-runif(segs, 0, 1)[short.seg]
  short.seg.rand<-short.seg.rand + (1:length(short))*0.000001
  short<-short[order(short.seg.rand)]
  
  short.adj<-short
  for(i in disConts){
    selSeg<-i:(i+(defaults$embedD1-1))
    this.seg<-short[selSeg]
    this.seg<-this.seg + long[selSeg[round(defaults$embedD1/2)]] 
    short.adj[selSeg]<-this.seg
  }
  

  
  
  
  # ----------------------------------
  # Scale error for a chosen signal to noise ratio.
  tESD <- 1 / SNR
  tOscData<-NULL
  # tEmbedded<-NULL
  print(head(sim))
  
  selCols<-seq(2, defaults$nVars*4+1, 4)
  for(i in selCols){
    tSD <- sd(short.adj)
    thisVar<-(short.adj)/tSD + rnorm(N, mean=0, sd=tESD)
    thisVar<- thisVar+ randMeans[i/2,2]*(1:length(thisVar)) + randMeans[i/2,1]
    tOscData <- cbind(tOscData, thisVar)
  }
  

  dimnames(tOscData) <- list(NULL, paste("X",1:defaults$nVars,sep=""))
  
  if(plotting){
  #   plotColors<-rainbow(nVars, v=0.75)
  #   plotColors.muted<-rainbow(nVars, s=1, v=0.75, alpha=0.5)
  #   plot(c(min(theTimes), max(theTimes[1:N])), c(-5, 5),
  #        xlab="Time",
  #        ylab="Score",
  #        type='n',
  #        main=paste("Multiscale Oscillator (Discontinuous)", sep=""))
  #   for(i in 1:nVars){
  #     lines(sim$time, tOscData[,i], type='l', lwd=1, col=plotColors.muted[i])
  #   }
  #   abline(v=disConts-1, lty=3, col="orange")
    
    # for(i in 1:nVars){
    #   lines(shocks[[i]], type="h", lwd=3, col=plotColors[i])
    # }
  
  plot( c(0, length(short)),c(0,0), type="l", 
        main="Discontinuous Multi-scale Series",
        xlab="Time",
        ylab="Score",
        col="darkgray",
        ylim=c(-4,4))
  for(i in disConts){
    selSeg<-i:(i+(defaults$embedD1-1))
    lines((0:length(short))[selSeg], short.adj[selSeg], col="blue")
  }
  abline(v=disConts-1, lty=3, col="orange")
  }
  
  return(tOscData)
}








