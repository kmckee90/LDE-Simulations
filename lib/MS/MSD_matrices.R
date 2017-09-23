
#------------------------------------------------------------------------
#Set up matrices for univariate DLO
#By Kevin McKee
#------------------------------------------------------------------------
defaults$embedD<-defaults$embedD1*defaults$embedD2


numIndicators<-1
nLatentVars<-3*(defaults$embedD2+1)


manifestVars <- NULL
for (i in 1:numIndicators) {
  manifestVars <-
    c(manifestVars, paste("x", i, "_", 0:(defaults$embedD-1), sep = ""))
}


DLO1.l<-c("X1", "dX1", "d2X1")
DLO2.l<-c("Xm1", "dXm1", "d2Xm1")

lV.l<-c(rep(DLO1.l, defaults$embedD2), DLO2.l)
# lV.l<-c(paste0(DLO1.l,1), paste0(DLO1.l,2), DLO2.l)

ind.l<-paste0("X",1:(defaults$embedD))
# L matrix ----------------------------------------------------------------
a.L1 <- rep(1,defaults$embedD)#/ 6.44915
a.L2 <- sort(rep(c(1:(defaults$embedD2))*defaults$theTau*defaults$deltaT-mean(c(1:(defaults$embedD2))*defaults$theTau*defaults$deltaT), defaults$embedD1))*defaults$embedD1
a.L3 <-  (a.L2^2)/2
# a.L2 <- a.L2 / 2.73269988

b.L1 <- rep(1,defaults$embedD1) #/ 0.89340
b.L2 <- c(1:defaults$embedD1)*defaults$theTau*defaults$deltaT-mean(c(1:defaults$embedD1)*defaults$theTau*defaults$deltaT)
b.L3 <-  (b.L2^2)/2
# b.L2 <- b.L2 / 0.990665

a.LMatrix <- cbind(a.L1,a.L2,a.L3)

b.LMatrix <- diag(defaults$embedD2)%x%cbind(b.L1,b.L2,b.L3)

LMatrix<-cbind(b.LMatrix, a.LMatrix)
dimnames(LMatrix)<-list( ind.l, lV.l)
# Z<-matrix(0, defaults$embedD, 6)
# LMatrix<-cbind(LMatrix, Z)


# A matrix ----------------------------------------------------------------

A.sub<-matrix(0,3,3)
A.sub[3,1]<-1
A.sub[3,2]<-2
Am<-diag(defaults$embedD2+1)%x%A.sub
dimnames(Am)<-list( lV.l, lV.l)
Am[seq(3, 3*(defaults$embedD2),3)-2, 3*(defaults$embedD2+1)-2]<-0
# Am[seq(3, 3*(defaults$embedD2),3),3*(defaults$embedD2+1)-1]<-3
# Am[seq(3, 3*(defaults$embedD2),3),3*(defaults$embedD2+1)-2]<-3

# Am[seq(3, 3*(defaults$embedD2),3),3*(defaults$embedD2+1)-2]<-6 #Level gamma. Unnecessary?



Am[3*(defaults$embedD2+1),3*(defaults$embedD2+1)-2]<-4
Am[3*(defaults$embedD2+1), 3*(defaults$embedD2+1)-1]<-5

A.l<-c("etaX1", "zetaX1", "gamma", "mEtaX1", "mZetaX")
Am.labels<-Am.free<-Am.sv<-Am.ub<-Am.lb<-Am
Am.labels[Am.labels==0]<-NA
for(i in 1:length(A.l)){Am.labels[Am.labels==i]<-A.l[i]}

Am.free<-Am>0
Am.sv<-ifelse(Am.sv==0, 0, -0.2)
Am.sv[grepl('mEtaX1', Am.labels)]<- -0.02
Am.lb<-ifelse(Am.lb==0, NA, -1)
Am.ub<-ifelse(Am.ub==0, NA, 1)


#S
S.sub<-diag(1,3)
S.sub[2,1]<-2
Sm<-diag(defaults$embedD2+1)%x%S.sub
for(i in 0:2){Sm[3*(defaults$embedD2+1)-i, 3*(defaults$embedD2+1)-i]<-3}
Sm[3*(defaults$embedD2+1)-1, 3*(defaults$embedD2+1)-2]<-4
dimnames(Sm)<-list( lV.l, lV.l)
Sm.labels<-Sm.free<-Sm.sv<-Sm.ub<-Sm.lb<-Sm


S.l<-c("VX1", "VdX1", "Vd2X1")
S.l.m<-c("VXM1", "VdXM1", "Vd2XM1")

Sm.labels<-ifelse(Sm==1, S.l, NA)
Sm.labels[Sm==2]<-"rVX1_VdX1"
Sm.labels[Sm==3]<-S.l.m
Sm.labels[Sm==4]<-"rVXM1_VdXM1"

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




LGC1 <- rep(1,defaults$embedD)
LGC2 <- (c(1:(defaults$embedD))*defaults$theTau*defaults$deltaT-mean(c(1:(defaults$embedD))*defaults$theTau*defaults$deltaT))
LGCMatrix <- rbind(LGC1,LGC2) #%x%matrix(1,1,2)
LGCMatrix2 <- rbind(rep(0,(defaults$embedD)), rep(defaults$deltaT,defaults$embedD)) #%x%matrix(1,1,2)




intLabs <- NULL
slopeLabs <- NULL
uLabs <- NULL
for (i in 1:numIndicators) {
  intLabs <-  c(intLabs, rep(paste0("intX", i), 1))
  slopeLabs <-  c(slopeLabs, rep(paste0("slopeX", i), 1))
  uLabs <-  c(uLabs, rep(paste0("Ux", i), defaults$embedD))
}
meanLabs<-cbind(intLabs, slopeLabs)

# save.image(file = "matrices.RData")

# save.image(file = "matrices.RData")
