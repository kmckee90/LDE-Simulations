# ---------------------------------------------------------------------
#   N-linked DLO with shocks
#   Author: Kevin McKee
#
# ---------------------------------------------------------------------
emptyDat<-data.frame()
model <- mxModel("model", type = "RAM",
               
               mxMatrix("Full", nrow=nTot, ncol=nTot,
                        values=AA.v, 
                        labels=AA.l, 
                        free=AA.f, 
                        ubound=AA.ub,
                        lbound=AA.lb, 
                        name="A", 
                        byrow=TRUE
               ),
               
               mxMatrix("Symm", nrow=nTot, ncol=nTot,
                        values=SS.v, 
                        labels=SS.l, 
                        free=SS.f, 
                        ubound=SS.ub,
                        lbound=SS.lb, 
                        name="S", 
                        byrow=FALSE
               ),
               
               mxMatrix("Full", nrow=nInd, ncol=nTot,
                        values=FF, 
                        name="F", 
                        byrow=TRUE
               ),
               mxMatrix("Full", nrow=nVars, ncol=2,
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
               mxMatrix("Full", nrow=1, ncol=nLVars, values=0, free=F, name="Z"),
               
               mxAlgebra(cbind(Z,t(rvectorize(Mu %*% ((J %*% H) + C)))),
                         dimnames=list(NULL, totLabels), name="M"
               ),
               
               mxExpectationRAM("A", "S", "F", "M"),
               mxFitFunctionML(),
               mxData(emptyDat, type="raw")
)