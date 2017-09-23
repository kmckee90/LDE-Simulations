# ---------------------------------------------------------------------
#   N-linked DLO with shocks
#   Author: Kevin McKee
#
# ---------------------------------------------------------------------

emptyDat<-data.frame()
model <- mxModel("model",
              mxMatrix("Full", values=c(defaults$embedD1,defaults$embedD2), nrow=1, ncol=2, name="embedD"),
                 
               mxMatrix("Full",  
                        values=LMatrix, 
                        free=FALSE, 
                        name="L", 
                        byrow=TRUE,
                        dimnames=list( manifestVars, lV.l)
               ),
               
               # mxAlgebra(-1*eta, name="gamma"),
               # mxConstraint(eta < etaM, name="separate_Etas"),
               
               mxMatrix("Full",
                        values=Am.sv, 
                        labels=Am.labels, 
                        free=Am.free, 
                        ubound=Am.ub,
                        lbound=Am.lb, 
                        name="A", 
                        byrow=TRUE,
                        dimnames=list( lV.l, lV.l)
               ),
               
               
               mxMatrix("Symm",nrow = nLatentVars, ncol = nLatentVars,
                        values=Sm.sv, 
                        labels=Sm.labels, 
                        free=Sm.free, 
                        ubound=Sm.ub,
                        lbound=Sm.lb, 
                        name="S", 
                        byrow=FALSE,
                        dimnames=list( lV.l, lV.l)
               ),
               
               mxMatrix("Diag", defaults$embedD*numIndicators, defaults$embedD*numIndicators, 
                        values=.8, 
                        free=TRUE, 
                        labels=uLabs, 
                        name="U",
                        lbound=0.000001
               ),
               
               
               
               mxMatrix("Full", nrow=numIndicators, ncol=2,
                        values=0.2, free=TRUE, labels=meanLabs, name="Mu"
               ),
               
               mxMatrix("Full", values=LGCMatrix, name="C"),
               mxMatrix("Full", values=LGCMatrix2, name="H"),
               mxMatrix("Full", nrow=2, ncol=2,
                        free=FALSE,
                        labels=c(NA, NA,
                                 NA, "data.Occasion"),
                        byrow=TRUE,
                        name="J"),
               
               mxAlgebra(t(rvectorize(Mu %*% ((J %*% H) + C))),
                         dimnames=list(NULL, manifestVars), name="M"
               ),
               
               
               
               mxMatrix("Iden", nLatentVars, name="I"),
               mxAlgebra(L %*% solve(I-A) %*% S %*% t(solve(I-A))  %*% t(L) + U, 
                         name="R", 
                         dimnames = list(manifestVars, manifestVars)
               ),
               
               
               
               mxExpectationNormal(covariance="R", means="M"),
               mxFitFunctionML(),
               mxData(emptyDat, 
                      type="raw"
               )
)