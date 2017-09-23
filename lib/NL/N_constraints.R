
#Model name is always DLO.

# model<-omxSetParameters(model, labels="intX1", values=0, free=F)
# model<-omxSetParameters(model, labels="slopeX1", values=0, free=F)
model<-omxSetParameters(model, labels=p.l, lbound=-1)
model<-omxSetParameters(model, labels=p.l, ubound=1)

if(defaults$nVars>1){
  model<-omxSetParameters(model, labels=gL.l[!diag(defaults$nVars)], lbound= -1)
  model<-omxSetParameters(model, labels=gS.l[!diag(defaults$nVars)], lbound= -1)
  model<-omxSetParameters(model, labels=gL.l[!diag(defaults$nVars)], ubound= 1)
  model<-omxSetParameters(model, labels=gS.l[!diag(defaults$nVars)], ubound= 1)
}