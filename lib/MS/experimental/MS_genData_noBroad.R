# ---------------------------------------------------------------------
#  Simulate a single Damped Linear Oscillator
#  Author: Original script by Steve Boker, modified by Kevin McKee
# ---------------------------------------------------------------------

# ----------------------------------
# Define the damped linear oscillator function.
DLOmodel <- function(t, prevState, parms) {
  x <- prevState[1] # x[t]
  y <- prevState[2] # dx[t]
  z <- prevState[3] # z[t]
  w <- prevState[4] # dz[t]
  
  with(as.list(parms), 
       {      
         dx <- y
         dy <- parms[1]*x + parms[2]*y  #+ parms[3]*z
         dz <- w
         dw <- parms[4]*z + parms[5]*w 
         res<-c(dx,dy,dz,dw)
         list(res)
       }
  )
}

# ---------------------------------------------
theTimes  <- 1:N  

#Generate correlated shocks
shocksX<-rbinom(N,1,pEvent)*eventScale*sample(c(-1,1),N, replace=T)
if(eventNormal){ shocksX<-shocksX*rnorm(N)}
shocksY<-shocksX*eventScaleRatio

# shocksX<-rbinom(N,1,pEvent)*eventScaleX1*sample(c(-1,1),N, replace=T)
# shocksY<-rbinom(N,1,pEvent)*eventScaleX2*sample(c(-1,1),N, replace=T)

  # shocksX<-shocksX*rnorm(N)


if(eventTypeX1=="slope"){
  shocks.m<-c(rep(0,N), shocksX)
}else if(eventTypeX1=="both"){
  shocks.m<-c(rep(shocksX,2))
}else{
  shocks.m<-c(shocksX, rep(0,N))
}


if(eventTypeX2=="slope"){
  shocks.m<-c(shocks.m, rep(0,N), shocksY)
}else if(eventTypeX2=="both"){
  shocks.m<-c(shocks.m, rep(shocksY,2))
}else{
  shocks.m<-c(shocks.m, shocksY, rep(0,N))
}



eventdat <- data.frame(var = c(rep("x", N),rep("y", N),rep("z", N),rep("w", N)),
                       time = rep(theTimes, 4),
                       value = shocks.m,
                       method = rep("add", N*4))


# ----------------------------------
# Simulate a damped linear oscillator.
parms <- c(etaX1, zetaX1, gammaX1, etaX2, zetaX2, gammaX2)
# tOffsets <- c(1:N)
tOffsets <- c(1:N)
xstart <- c(x = initCond.X, y = initCond.dX, z = initCond.Y, w = initCond.dY)
out1 <- as.data.frame(lsoda(xstart, theTimes, DLOmodel, parms, events = list(data=eventdat)))



# eventTimes.collapse<-which(shocks!=0)

# print(eventTimes.collapse)
# print(varNames)
# print(theTimes)

# if(length(eventTimes.collapse)>0){
#   theTimes.micro<- c(eventTimes.collapse+0.000001)
#   theTimes.event<-sort(c(theTimes, theTimes.micro))
# 
#   # print(theTimes.event)
#   # print(eventdat)
# 
#   try(pred <- as.data.frame(ode(xstart, theTimes.event, DLOmodel, parms,
#                                 events = list(data=eventdat), method="lsoda")))
# 
# 
#   pred[which(pred$time%in%eventTimes.collapse),-1]<-pred[which(pred$time%in%theTimes.micro),-1]
#   pred<-pred[-which(pred$time%in%theTimes.micro),]
# }else{
#   try(pred <- as.data.frame(ode(xstart, theTimes, DLOmodel, parms,
#                                 events = list(data=eventdat), method="lsoda")))
# }
# # print(pred)
# 



# out1<-pred











# ----------------------------------
# Scale error for a chosen signal to noise ratio.





tSD <- sqrt(var(c(out1$x,out1$z)))
tESD <- 1 / tSNR
tOscData <- out1$x/tSD + rnorm(N, mean=0, sd=tESD)
# tOscDataY <- out1$z/tSD



# series<-list("LDE"=out1, "Data"=tOscData, "EventsX"=shocksX,"EventsY"=shocksY,"EventsR"=shocksR, "Info"=simInfo)

# ----------------------------------
# Plot the true displacement and first derivative plus the noisy signal.


plot(c(min(theTimes), max(theTimes[1:N])), c(-5, 5),
     xlab="Time",
     ylab="Score",
     type='n',
     main=paste("Single Multiscale Oscillator, rShocks=",pEventX1, sep=""))
lines(out1$time, tOscData, type='l', lwd=1, col='lightblue')
lines(c(0,shocksX), type="h",lwd=3,col="blue")
lines(c(0,shocksY), type="h",lwd=3,col="red")

lines(out1$time, out1$x/tSD, type='l', lwd=1, col='blue')
# lines(out1$time, out1$y/tSD, type='l', lwd=2, col='blue')
lines(out1$time, out1$z/tSD, type='l', lwd=1, col='red')
# lines(out1$time, out1$w/tSD, type='l', lwd=2, col='red')

lines(c(min(theTimes), max(theTimes[1:N])), c(-0, 0), type='l', lty=2, col=1)

plot(tOscData, type='l')

tData<-tOscData
tEmbedded<- gllaEmbed(tData, embed=embedD, tau=theTau, idColumn=FALSE)
colnames(tEmbedded)<-manifestVars
Occasion<-1:nrow(tEmbedded) + round(embedD/2)
tEmbedded<-cbind(tEmbedded, Occasion)


