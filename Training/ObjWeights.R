
source("klmFuncs.R")
# Get Weights for Mixed Objective Function --------------------------------
normTest.dat1<-as.matrix(read.csv("output/normTest.dat.1.csv",header=T))
normTest.dat2<-as.matrix(read.csv("output/normTest.dat.2.csv",header=T))
normTest.dat<-rbind(normTest.dat1, normTest.dat2)
dim(normTest.dat)

normTest.dat[abs(normTest.dat[,1])>10,1]<-NA
normTest.dat[abs(normTest.dat[,2])>2,2]<-NA
normTest.dat[abs(normTest.dat[,5])>0.5,5]<-NA
normTest.dat[abs(normTest.dat[,6])>2,6]<-NA
normTest.dat[abs(normTest.dat[,7])>1,7]<-NA

normTest.dat<-normTest.dat[,-(2:4)]
normTest.dat<-normTest.dat[complete.cases(normTest.dat),]
normTest.dat.var<-diag(var(normTest.dat))[-1]
normTest.dat.mu<-colMeans(normTest.dat)[-1]
normTest.z<-apply(normTest.dat, 2, stdz)
normTest<-glm(normTest.z[,1]~normTest.z[,-1])

ObjEA<-list("objWeights"=as.matrix(normTest$coefficients[-1]), 
                    "objVars"=normTest.dat.var, 
                    "objMeans"=normTest.dat.mu)
save(ObjEA, file="ObjectiveWeightsLM.RData")
