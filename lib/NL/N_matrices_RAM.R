

try(dat<-tEmbedded)

# RAM matrices ------------------------------------------------------------
nInd<-nVars*embedD
nLVars<-nVars*3
nTot<-nInd+nLVars
SS<-diag(1, nTot)
AA<-matrix(0, nTot, nTot)
  
ind.labels<-colnames(dat)[1:nInd]
latentVar.labels<-matrix(c("VX","VdX","Vd2X"),nVars,3, byrow=T)
latentVar.labels<-matrix(paste(latentVar.labels,1:nVars,sep=""),nVars,3)
lV.l<-as.vector(t(latentVar.labels))
totLabels<-c(lV.l, ind.labels)

dimnames(SS)<-dimnames(AA)<-list(totLabels, totLabels)

sel.VX<-round(seq(1, nLVars, 3))
sel.VdX<-round(seq(2, nLVars, 3))
sel.Vd2X<-round(seq(3, nLVars, 3))


# S matrix ----------------------------------------------------------------
SS.l<-ifelse(SS==1, totLabels, NA)
for(i in 1:nVars){
  SS.l[sel.VX[-i],sel.VX[i]]<-paste("r",lV.l[sel.VX[i]],"_",lV.l[sel.VX[-i]],sep="")
  SS.l[sel.VdX[i],sel.VX[i]]<-paste("r",lV.l[sel.VX[i]],"_",lV.l[sel.VdX[i]],sep="")
  
  SS.l[grepl(paste0("x",i), SS.l)]<-paste0("Ux",i)
  
}

SS.v<-SS*0.8
SS.v[grepl("Node", SS.l)]<-0

SS.f<-!is.na(SS.l)
SS.f[grepl("Node", SS.l)]<-F

SS.lb<-SS.ub<-matrix(NA, nTot, nTot)
SS.lb<-SS*0.00001
SS.lb[SS.lb==0]<-NA
SS.lb[grepl("Node", SS.l)]<-NA


# SS.lb[grepl("rV", SS.l)]<- -1
# SS.ub[grepl("r", SS.l)]<-  1
SS.v[grepl("rVX", SS.l)]<- -0.2

# A matrix ----------------------------------------------------------------
par.labels<-matrix(c("etaX","zetaX"),nVars,2, byrow=T)
p.l<-matrix(paste(par.labels,1:nVars,sep=""),nVars,2)
# par.labels<-as.vector(t(par.labels))

cm.i<-matrix(1:nVars, nVars,nVars)
cm<-matrix(paste(t(cm.i),cm.i,sep=""),nVars,nVars)
g.l<-matrix(paste(matrix("gammaX", nVars, nVars), cm ,sep=""),nVars, nVars)

cm.i<-matrix(1:nVars, nVars,nVars)
cm<-matrix(paste(t(cm.i),cm.i,sep=""),nVars,nVars)
g.l<-matrix(paste(matrix("gammaX", nVars, nVars), cm ,sep=""),nVars, nVars)

gL.l<-matrix(paste0(g.l,"L"),nVars,nVars)
gS.l<-matrix(paste0(g.l,"S"),nVars,nVars)

AA.l<-AA
AA.l[AA.l==0]<-NA
AA.lb<-AA.ub<-AA.l
  
for(i in 1:nVars){
  AA.l[sel.Vd2X[i],sel.VX[i]]<-p.l[i,1]
  AA.l[sel.Vd2X[i],sel.VdX[i]]<-p.l[i,2]
  
  AA.l[sel.Vd2X[-i],sel.VX[i]]<-gL.l[-i,i]
  AA.l[sel.Vd2X[-i],sel.VdX[i]]<-gS.l[-i,i]
}

# L submatrix ----------------------------------------------------------------
L1 <- rep(1,embedD)
L2 <- c(1:embedD)*theTau*deltaT-mean(c(1:embedD)*theTau*deltaT)
L3 <-  (L2^2)/2
LMatrix <- cbind(L1,L2,L3)
# LMatrix<-cbind(LMatrix, rep(0, embedD))
LMatrix<-diag(nVars) %x% LMatrix 

AA.v<-rbind(matrix(0,nLVars, nTot) , cbind(LMatrix, matrix(0,nInd,nInd)))
dimnames(AA.v)<-list(totLabels, totLabels)
AA.v[grepl("eta", AA.l)]<- -0.2
AA.v[grepl("gamma", AA.l)]<- -0.2

AA.f<-!is.na(AA.l)
AA.lb[AA.f]<- -10
AA.ub[AA.f]<-  10





#  ------------------------------------------------------------------------
SS.l<-SS.l[lower.tri(SS.l,diag=T)]
SS.v<-SS.v[lower.tri(SS.v,diag=T)]
SS.f<-SS.f[lower.tri(SS.f,diag=T)]
SS.ub<-SS.ub[lower.tri(SS.ub,diag=T)]
SS.lb<-SS.lb[lower.tri(SS.lb,diag=T)]



# F Matrix ----------------------------------------------------------------
FF<-cbind(matrix(0,nInd,nLVars), diag(1, nInd, nInd))
dimnames(FF)<-list(ind.labels, totLabels)


# Nonstationary means -----------------------------------------------------
LGC1 <- rep(1,embedD)
LGC2 <- (c(1:embedD)*theTau*deltaT-mean(c(1:embedD)*theTau*deltaT))
LGCMatrix <- rbind(LGC1,LGC2) #%x%matrix(1,1,2)
LGCMatrix2 <- rbind(rep(0,embedD), rep(deltaT,embedD)) #%x%matrix(1,1,2)


manifestVars <- NULL
for (i in 1:nVars) {
  manifestVars <-
    c(manifestVars, paste("x", i, "_", 0:(embedD - 1), sep = ""))
}

intLabs <- NULL
slopeLabs <- NULL
uLabs <- NULL
for (i in 1:nVars) {
  intLabs <-  c(intLabs, rep(paste0("intX", i), 1))
  slopeLabs <-  c(slopeLabs, rep(paste0("slopeX", i), 1))
  uLabs <-  c(uLabs, rep(paste0("uX", i), embedD))
}
meanLabs<-cbind(intLabs, slopeLabs)

# save.image(file = "matrices.RData")
