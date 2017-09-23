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
         dy <- parms[1]*x + parms[2]*y
         dz <- w
         dw <- parms[4]*z + parms[5]*w 
         res<-c(dx,dy,dz,dw)
         list(res)
       }
  )
}

# ---------------------------------------------
theTimes  <- 1:N  

disConts<-seq(1, N,embedD1)

nEvent<-round(runif(1, 1, length(disConts)*pEvent-2))
shocks.t<-sort(sample(disConts[-1], nEvent, replace=F))
shocks.v<-rnorm(length(shocks.t),0,eventScale)

eventdat <- data.frame(var = rep("x", nEvent),
                       time = shocks.t-1,
                       value = shocks.v,
                       method = rep("add", nEvent))


# ----------------------------------
# Simulate a damped linear oscillator.
parms <- c(etaX1, zetaX1, gammaX1, etaX2, zetaX2, gammaX2)
# tOffsets <- c(1:N)
tOffsets <- c(1:N)
xstart <- c(x = initCond.X, y = initCond.dX, z = initCond.Y, w = initCond.dY)
out1 <- as.data.frame(lsoda(xstart, theTimes, DLOmodel, parms, events = list(data=eventdat)))



short<-sim$x
long<-sim$z

segs<-round(length(short)/embedD1)
short.seg<-sort(rep(1:segs, embedD1))
short.seg.rand<-runif(segs, 0, 1)[short.seg]
short.seg.rand<-short.seg.rand + (1:length(short))*0.000001
short<-short[order(short.seg.rand)]

short.adj<-short
for(i in disConts){
  selSeg<-i:(i+(embedD1-1))
  this.seg<-short[selSeg]
  this.seg<-this.seg + long[selSeg[round(embedD1/2)]] #- mean(this.seg) 
  short.adj[selSeg]<-this.seg
}


# plot(short, type="l", 
#      main="Discontinuous Multi-scale Series",
#      xlab="Time",
#      ylab="Score",
#      col="darkgray")
# lines(short.adj, col="blue")
# abline(v=disConts, lty=3, col="orange")


tSD <- sqrt(var(c(out1$x,out1$z)))
tESD <- 1 / tSNR
short.adj <- short.adj/tSD + rnorm(N, mean=0, sd=tESD)
# tOscDataY <- out1$z/tSD

if(genPlotting){
  plot( c(0, length(short)),c(0,0), type="l", 
        main="Discontinuous Multi-scale Series",
        xlab="Time",
        ylab="Score",
        col="darkgray",
        ylim=c(-4,4))
  for(i in disConts){
    selSeg<-i:(i+(embedD1-1))
    lines((0:length(short))[selSeg], short.adj[selSeg], col="blue")
  }
  abline(v=disConts-1, lty=3, col="orange")
}


# series<-list("LDE"=out1, "Data"=tOscData, "EventsX"=shocksX,"EventsY"=shocksY,"EventsR"=shocksR, "Info"=simInfo)

# ----------------------------------
# Plot the true displacement and first derivative plus the noisy signal.
# plot(c(min(theTimes), max(theTimes[1:N])), c(-5, 5),
#      xlab="Time",
#      ylab="Score",
#      type='n',
#      main=paste("Single Multiscale Oscillator, rShocks=",pEventX1, sep=""))
# lines(out1$time, tOscData, type='l', lwd=1, col='lightblue')
# lines(c(0,shocksX), type="h",lwd=3,col="blue")
# lines(c(0,shocksY), type="h",lwd=3,col="red")
# 
# lines(out1$time, out1$x/tSD, type='l', lwd=1, col='blue')
# # lines(out1$time, out1$y/tSD, type='l', lwd=2, col='blue')
# lines(out1$time, out1$z/tSD, type='l', lwd=1, col='red')
# # lines(out1$time, out1$w/tSD, type='l', lwd=2, col='red')
# 
# lines(c(min(theTimes), max(theTimes[1:N])), c(-0, 0), type='l', lty=2, col=1)
# 




# plot(tOscData, type='l')

tData<-short.adj

tEmbedded<- gllaEmbed(tData, embed=embedD, tau=theTau, idColumn=FALSE)
# weekly<-seq(1, length(disConts),embedD2)
tEmbedded<-tEmbedded[disConts[disConts<nrow(tEmbedded)],]

colnames(tEmbedded)<-manifestVars
Occasion<-1:nrow(tEmbedded) + round(embedD/2)
tEmbedded<-cbind(tEmbedded, Occasion)


