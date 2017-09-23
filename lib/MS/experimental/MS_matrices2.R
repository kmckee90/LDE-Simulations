
#------------------------------------------------------------------------
#Set up matrices for univariate DLO
#By Kevin McKee
#------------------------------------------------------------------------
# embedD<-embedD1*embedD2
embedD<-embedD1+embedD2-1

manifestVars <- NULL
for (i in 1:numIndicators) {
  manifestVars <-
    c(manifestVars, paste("x", i, "_", 0:(embedD-1), sep = ""))
}


numIndicators<-1
nLatentVars<-3*(embedD2+1)

DLO1.l<-c("X", "dX", "d2X")
DLO2.l<-c("Xm", "dXm", "d2Xm")

lV.l<-c(rep(DLO1.l, embedD2), DLO2.l)
# lV.l<-c(paste0(DLO1.l,1), paste0(DLO1.l,2), DLO2.l)

ind.l<-paste0("X",1:(embedD))
# L matrix ----------------------------------------------------------------
a.L1 <- rep(1,embedD)
a.L2 <- c(1:(embedD))*theTau*deltaT-mean(c(1:(embedD))*theTau*deltaT)
a.L3 <-  (a.L2^2)/2

b.L1 <- rep(1,embedD1)
b.L2 <- c(1:embedD1)*theTau*deltaT-mean(c(1:embedD1)*theTau*deltaT)
b.L3 <-  (b.L2^2)/2

a.LMatrix <- cbind(a.L1,a.L2,a.L3)
b.LMatrix <- cbind(b.L1,b.L2,b.L3)
c.LMatrix <- NULL
for(i in 1:embedD2){
  b.above<- matrix(0, (i-1), 3)
  b.below<- matrix(0, embedD2-(i), 3)
  b.LMatrix.col <- rbind(b.above, b.LMatrix, b.below)
  c.LMatrix <- cbind(c.LMatrix, b.LMatrix.col )
}

LMatrix<-cbind(c.LMatrix, a.LMatrix)
dimnames(LMatrix)<-list( ind.l, lV.l)
# Z<-matrix(0, embedD, 6)
# LMatrix<-cbind(LMatrix, Z)


# A matrix ----------------------------------------------------------------

A.sub<-matrix(0,3,3)
A.sub[3,1]<-1
A.sub[3,2]<-2
Am<-diag(embedD2+1)%x%A.sub
dimnames(Am)<-list( lV.l, lV.l)
Am[seq(3, 3*(embedD2),3),3*(embedD2+1)-1]<-3
# Am[seq(3, 3*(embedD2),3),3*(embedD2+1)-2]<-6 #Level gamma. Unnecessary?



Am[3*(embedD2+1),3*(embedD2+1)-2]<-4
Am[3*(embedD2+1), 3*(embedD2+1)-1]<-5

A.l<-c("eta", "zeta", "gamma", "etaM", "zetaM")
Am.labels<-Am.free<-Am.sv<-Am.ub<-Am.lb<-Am
Am.labels[Am.labels==0]<-NA
for(i in 1:length(A.l)){Am.labels[Am.labels==i]<-A.l[i]}

Am.free<-Am>0
Am.sv<-ifelse(Am.sv==0, 0, -0.2)
Am.lb<-ifelse(Am.lb==0, NA, -10)
Am.ub<-ifelse(Am.ub==0, NA, 1.5)


#S
S.sub<-diag(1,3)
S.sub[2,1]<-2
Sm<-diag(embedD2+1)%x%S.sub
for(i in 0:2){Sm[3*(embedD2+1)-i, 3*(embedD2+1)-i]<-3}
Sm[3*(embedD2+1)-1, 3*(embedD2+1)-2]<-4
dimnames(Sm)<-list( lV.l, lV.l)
Sm.labels<-Sm.free<-Sm.sv<-Sm.ub<-Sm.lb<-Sm


S.l<-c("VX", "VdX", "Vd2X")
S.l.m<-c("VXM", "VdXM", "Vd2XM")

Sm.labels<-ifelse(Sm==1, S.l, NA)
Sm.labels[Sm==2]<-"rVX_VdX"
Sm.labels[Sm==3]<-S.l.m
Sm.labels[Sm==4]<-"rVXM_VdXM"

Sm.free<-Sm>0
Sm.sv<-ifelse(Sm==1 | Sm==3, 0.8, 0)
Sm.sv[Sm==2 | Sm==4]<- -0.1

Sm.lb<-ifelse(Sm==2 | Sm==4, -1, NA)
Sm.lb[Sm==1|Sm==3]<-0.0001
Sm.ub<-ifelse(Sm==2 | Sm==4, 1, NA)

Sm.labels<-Sm.labels[lower.tri(Sm.labels, diag=T)]
Sm.sv<-Sm.sv[lower.tri(Sm.sv, diag=T)]
Sm.free<-Sm.free[lower.tri(Sm.free, diag=T)]
Sm.lb<-Sm.lb[lower.tri(Sm.lb, diag=T)]
Sm.ub<-Sm.ub[lower.tri(Sm.ub, diag=T)]




LGC1 <- rep(1,embedD)
LGC2 <- (c(1:(embedD))*theTau*deltaT-mean(c(1:(embedD))*theTau*deltaT))
LGCMatrix <- rbind(LGC1,LGC2) #%x%matrix(1,1,2)
LGCMatrix2 <- rbind(rep(0,(embedD)), rep(deltaT,embedD)) #%x%matrix(1,1,2)




intLabs <- NULL
slopeLabs <- NULL
uLabs <- NULL
for (i in 1:numIndicators) {
  intLabs <-  c(intLabs, rep(paste0("intX", i), 1))
  slopeLabs <-  c(slopeLabs, rep(paste0("slopeX", i), 1))
  uLabs <-  c(uLabs, rep(paste0("Ux", i), embedD))
}
meanLabs<-cbind(intLabs, slopeLabs)

# save.image(file = "matrices.RData")

# save.image(file = "matrices.RData")
