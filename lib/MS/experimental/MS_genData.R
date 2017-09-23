# ---------------------------------------------------------------------
#  Simulate a single Damped Linear Oscillator
#  Author: Original script by Steve Boker, modified by Kevin McKee
# ---------------------------------------------------------------------

# ----------------------------------
# Define the damped linear oscillator function.
DLOmodel <- function(t, prevState, parms) {
  dx<-matrix(prevState, nVars, 2, byrow=T)
  states<-matrix(NA, nVars, 2)
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



# ---------------------------------------------
theTimes  <- 1:N
varNames<-as.vector(rbind(paste("x",1:nVars, sep=""), paste("dx",1:nVars, sep="")))

#Generate correlated shocks
shocks<-list()
for(i in 1:nVars){
  shocks[[i]]<-eventScale*rbinom(N,1,pEvent)*sample(c(-1, 1), N, replace=T)
  if(eventNormal){
    shocks[[i]]<-shocks[[i]]*rnorm(N)
  }
}
# shocks<-rep(1,N*nVars*2)

shocks.m<-NULL
for(i in 1:nVars){
  if(eventType=="slope"){
    shocks.m<-c(shocks.m, rep(0,N), shocks[[i]])
  }else if(eventType=="both"){
    shocks.m<-c(shocks.m, rep(shocks[[i]],2))
  }else{
    shocks.m<-c(shocks.m, shocks[[i]], rep(0,N))
  }
}


varNames.vec<-NULL
for(i in varNames){
  varNames.vec<-c(varNames.vec, rep(i, N))
}

eventdat <- data.frame(var = varNames.vec,
                       time = rep(theTimes, nVars*2),
                       value = shocks.m,
                       method = rep("add", N*nVars*2))


parms<-list(indDyn, gammas)
# tOffsets <- c(1:N)
tOffsets <- c(1:N)
xstart <- as.vector(t(initCond))
names(xstart)<-varNames
# xstart[c(2,4,6)]<-0.5
out1 <- as.data.frame(lsoda(xstart, theTimes, DLOmodel, parms, events = list(data=eventdat)))#[tOffsets,]

# ----------------------------------
# Scale error for a chosen signal to noise ratio.

tESD <- 1 / tSNR
tOscData<-NULL
tEmbedded<-NULL

if(!stationary){
  randInt<-rnorm(nVars)
  randSlope<-rnorm(nVars,0,0.1)
  randMeans<-cbind(randInt, randSlope)
}else{
  randInt<-rep(0, nVars)
  randSlope<-rep(0, nVars)
  randMeans<-cbind(randInt, randSlope)}


for(i in (1:nVars)*2){
  tSD <- sd(out1[,i])
  thisVar<-out1[,i]/tSD + rnorm(N, mean=0, sd=tESD)
  
  thisVar<- thisVar+ randMeans[i/2,2]*(1:length(thisVar)) + randMeans[i/2,1]
  
  tOscData <- cbind(tOscData, thisVar)
  tEmbedded <- cbind(tEmbedded, gllaEmbed(thisVar, embed=embedD, tau=theTau, idColumn=FALSE))
}

tEmbedded<-cbind(tEmbedded, "Occasion"=1:nrow(tEmbedded))

dimnames(tOscData) <- list(NULL, paste("X",1:nVars,sep=""))
colnames(tEmbedded)[-(nVars*embedD+1)]<-manifestVars


# series<-list("LDE"=out1, "Data"=tOscData, "EventsX"=shocksX1r,"EventsY"=shocksX2r,"EventsR"=shocksR, "Info"=simInfo)


if(genPlotting==T){
  plotColors<-rainbow(nVars, v=0.75)
  plotColors.muted<-rainbow(nVars, s=1, v=0.75, alpha=0.5)
  
  plot(c(min(theTimes), max(theTimes[1:N])), c(-5, 5),
       xlab="Time",
       ylab="Score",
       type='n')
  # main=paste(nVars,"-Oscillator Network", sep=""))
  for(i in 1:nVars){
    lines(out1$time, tOscData[,i], type='l', lwd=1, col=plotColors.muted[i])
  }
  for(i in 1:nVars){  
    lines(shocks[[i]], type="h", lwd=3, col=plotColors[i])
  }
  # title(main="Family Stress Dynamics")
  # legend("topright",
  #        c("Person1", "Person2", "Person3", "Person4", "Person5"),
  #        col=plotColors.muted,
  #        lty=1,
  #        lwd=4
  #        )
  
}

# tData<-tOscData



