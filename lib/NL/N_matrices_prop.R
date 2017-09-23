
#------------------------------------------------------------------------
#Set up matrices for univariate DLO
#By Kevin McKee
#------------------------------------------------------------------------

numIndicators<-nVars
nLatentVars<-nVars*3+nVars
# L matrix ----------------------------------------------------------------
L1 <- rep(1,embedD)
L2 <- c(1:embedD)*theTau*deltaT-mean(c(1:embedD)*theTau*deltaT)
L3 <-  (L2^2)/2
LMatrix <- cbind(L1,L2,L3)
LMatrix<-cbind(LMatrix, rep(0, embedD))
LMatrix<-diag(nVars) %x% LMatrix 


# A matrix ----------------------------------------------------------------
latentVar.labels<-matrix(c("VX","VdX","Vd2X","NodeX"),nVars,4, byrow=T)
latentVar.labels<-matrix(paste(latentVar.labels,1:nVars,sep=""),nVars,4)
lV.l<-as.vector(t(latentVar.labels))

par.labels<-matrix(c("etaX","zetaX"),nVars,2, byrow=T)
p.l<-matrix(paste(par.labels,1:nVars,sep=""),nVars,2)
# par.labels<-as.vector(t(par.labels))

cm.i<-matrix(1:nVars, nVars,nVars)
cm<-matrix(paste(t(cm.i),cm.i,sep=""),nVars,nVars)
g.l<-matrix(paste(matrix("gammaX", nVars, nVars), cm ,sep=""),nVars, nVars)



# A matrix ----------------------------------------------------------------
Am.labels<-matrix(NA, nLatentVars,nLatentVars)
rownames(Am.labels)<-colnames(Am.labels)<-lV.l

Am.fixed<-NULL
for(i in 0:(nVars-1)){
  Am.labels[4+i*4, 1+i*4]<-p.l[i+1,1]
  Am.labels[4+i*4, 2+i*4]<-p.l[i+1,2]
  Am.labels[3+i*4, ((1:nVars)*4)]<-g.l[i+1,]
  # Am.fixed<-c(Am.fixed, )
}
Am.fixed<-which(Am.labels%in%diag(g.l))
Am.freepars<-which(!is.na(Am.labels))
Am.freepars<-Am.freepars[!Am.freepars%in%Am.fixed]

Am.free<-matrix(F, nLatentVars,nLatentVars)
Am.free[Am.freepars]<-T

Am.sv<-matrix(0, nLatentVars,nLatentVars)
Am.sv[Am.freepars]<-  -0.2
Am.sv[Am.fixed] <- 1

Am.ub<-matrix(NA, nLatentVars,nLatentVars)
Am.ub[grepl("eta", Am.labels)]<- -0.00001
# Am.ub[]<- -0.00001
# Am.ub[]<- 1.2
Am.lb<-matrix(NA, nLatentVars,nLatentVars)
# Am.lb[A.free]<- -1.2


#S
# Sm.ind<-matrix(1:(nLatentVars^2), nLatentVars,nLatentVars)
# Sm.ind[lower.tri(Sm.ind,diag=T)]



Sm.labels<-matrix(NA, nLatentVars,nLatentVars)
rownames(Sm.labels)<-colnames(Sm.labels)<-lV.l
diag(Sm.labels)<-lV.l


for(i in 0:(nVars-1)){
  Sm.labels[2+i*4, 1+i*4]<-paste("r",latentVar.labels[i+1,1],"_",latentVar.labels[i+1,2],sep="")
  Sm.labels[4+i*4, 4+i*4]<-NA
  Sm.labels[1+i*4, (0:(nVars-1)*4+1)[-(i+1)]]<-paste("r",latentVar.labels[-(i+1),1],"_",latentVar.labels[i+1,1],sep="")
}

Sm.freepars<-which(!is.na(Sm.labels))
Sm.freepars.r<-which(!Sm.labels%in%diag(Sm.labels))

Sm.sv<-matrix(0, nLatentVars,nLatentVars)
Sm.sv[Sm.freepars]<-  0.8
Sm.sv[Sm.freepars.r]<- -0.1


Sm.free<-matrix(F, nLatentVars,nLatentVars)
Sm.free[Sm.freepars]<-T

Sm.ub<-matrix(NA, nLatentVars,nLatentVars)
Sm.ub[Sm.freepars.r]<- 1

Sm.lb<-matrix(NA, nLatentVars,nLatentVars)
Sm.lb[Sm.freepars]<-0.00000001
Sm.lb[Sm.freepars.r]<- -1


Sm.labels<-Sm.labels[lower.tri(Sm.labels,diag=T)]
Sm.sv<-Sm.sv[lower.tri(Sm.sv,diag=T)]
Sm.free<-Sm.free[lower.tri(Sm.free,diag=T)]
Sm.ub<-Sm.ub[lower.tri(Sm.ub,diag=T)]
Sm.lb<-Sm.lb[lower.tri(Sm.lb,diag=T)]

# save.image(file="MatricesUnivariate.RData")

LGC1 <- rep(1,embedD)
LGC2 <- (c(1:embedD)*theTau*deltaT-mean(c(1:embedD)*theTau*deltaT))
LGCMatrix <- rbind(LGC1,LGC2) #%x%matrix(1,1,2)
LGCMatrix2 <- rbind(rep(0,embedD), rep(deltaT,embedD)) #%x%matrix(1,1,2)


manifestVars <- NULL
for (i in 1:numIndicators) {
  manifestVars <-
    c(manifestVars, paste("x", i, "_", 0:(embedD - 1), sep = ""))
}

intLabs <- NULL
slopeLabs <- NULL
uLabs <- NULL
for (i in 1:numIndicators) {
  intLabs <-  c(intLabs, rep(paste0("intX", i), 1))
  slopeLabs <-  c(slopeLabs, rep(paste0("slopeX", i), 1))
  uLabs <-  c(uLabs, rep(paste0("uX", i), embedD))
}
meanLabs<-cbind(intLabs, slopeLabs)

# save.image(file = "matrices.RData")
