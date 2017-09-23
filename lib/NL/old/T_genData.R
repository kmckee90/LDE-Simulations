# ---------------------------------------------------------------------
#  Simulate a single Damped Linear Oscillator
#  Author: Original script by Steve Boker, modified by Kevin McKee
# ---------------------------------------------------------------------

# ----------------------------------

#Generate correlated parameters between twin 1 and twin 2.
#Let's say rMZ will be .7
#Going to need to constrain etas and zetas to be below 0
#Gammas can be anything.


sig<-diag(6)
colnames(sig)<-rownames(sig)<-c("T1.eta","T1.zeta","T2.eta","T2.zeta","T1.gamma","T2.gamma")
sig[3,1]<-sig[1,3]<-rEta; sig[4,2]<-sig[2,4]<-rZeta
sig<-sig*0.01; sig[5,5]<-sig[6,6]<-0.05
mu<-c(-0.4, -0.2, -0.4, -0.2, 0, 0)
T.pars<-mvrnorm(nPairs, mu, Sigma = sig, empirical=T)

T.pars[,c(1,3)][which(T.pars[,c(1,3)]>0.1)]<-0.1 #Eta 0 correction
T.pars[,c(2,4)][which(T.pars[,c(2,4)]>0)]<-0 #Zeta 0 correction

# hist(T.pars[,6])
# plot(T.pars[,2],T.pars[,4])


twinData<-list()


for(i in 1:nPairs){
  
# Define the damped linear oscillator function.
  DLOmodel <- function(t, prevState, parms) {
    x <- prevState[1] # x[t]
    y <- prevState[2] # dx[t]
    z <- prevState[3] # z[t]
    w <- prevState[4] # dz[t]
    
    with(as.list(parms), 
         {      
           dx <- y
           dy <- parms[1]*x + parms[2]*y  + parms[3] * (parms[4]*z + parms[5]*w)
           dz <- w
           dw <- parms[4]*z + parms[5]*w  + parms[6] * (parms[1]*x + parms[2]*y)
           res<-c(dx,dy,dz,dw)
           list(res)
         }
    )
  }

# ---------------------------------------------
theTimes  <- 1:N  

#Generate correlated shocks
shocksX<-rbinom(N,1,pEvent/2)*sample(c(-1, 1), N, replace=T)#*rnorm(N, 0, eventScaleX1)
shocksY<-rbinom(N,1,pEvent/2)*sample(c(-1, 1), N, replace=T)#*rnorm(N, 0, eventScaleX2)
shocksR<-rbinom(N,1,pEventShared)*sample(c(-1, 1), N, replace=T)#*rnorm(N, 0, eventScaleShared)
shocksR.ind<-which(shocksR!=0)
shocksX[shocksR.ind]<-shocksR[shocksR.ind]
shocksY[shocksR.ind]<-shocksR[shocksR.ind]

locShocks<-unique(sort(c(which(shocksX!=0), which(shocksY!=0))+1))


if(eventTypeX1=="slope"){
  shocks.m<-c(rep(0,N), shocksX)
}else if(eventTypeX1=="both"){
  shocks.m<-c(rep(shocksX,2))
}else{
  shocks.m<-c(shocksX, rep(0,N))
}


if(eventTypeX2=="slope"){
  shocks.m<-c(shocks.m, rep(0,N), shocksY)
}else if(eventTypeX1=="both"){
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
parms <- c(T.pars[i,1], T.pars[i,2], T.pars[i,5], T.pars[i,3], T.pars[i,4], T.pars[i,6])
# tOffsets <- c(1:N)
tOffsets <- c(1:N)
xstart <- c(x = initCond.X, y = initCond.dX, z = initCond.Y, w = initCond.dY)
out1 <- as.data.frame(lsoda(xstart, theTimes, DLOmodel, parms, events = list(data=eventdat)))[tOffsets,]

# ----------------------------------
# Scale error for a chosen signal to noise ratio.

tSD <- sqrt(var(c(out1$x,out1$z)))
tESD <- 1 / tSNR
tOscDataX <- out1$x/tSD + rnorm(N, mean=0, sd=tESD)
tOscDataY <- out1$z/tSD + rnorm(N, mean=0, sd=tESD)

tOscData <- cbind(tOscDataX, tOscDataY)
dimnames(tOscData) <- list(NULL, c("X1", "X2"))



tData<-tOscData
tEmbeddedX <- gllaEmbed(tData[,1], embed=embedD, tau=theTau, idColumn=FALSE)
tEmbeddedY <- gllaEmbed(tData[,2], embed=embedD, tau=theTau, idColumn=FALSE)
tEmbedded.cur<-data.frame(cbind(tEmbeddedX, tEmbeddedY))
colnames(tEmbedded.cur)<-manifestVars

twinData[[i]]<-tEmbedded.cur
}


