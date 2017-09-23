source("lib/klmFuncs.r")

expectedValsML<-function(model){
  # Factor Scores -----------------------------------------------------------
  scores<-as.matrix(as.data.frame(mxFactorScores(model$fit, type='WeightedML', minManifests = 7)))[,1:(3*nVars)] 
  scores.pred<- scores[,1:2]%*%t(LMatrix[,1:2])
  scores.shocks<-scores[,3]
  scores.filled<-insertRow(scores.pred, NA, model$exclusions)
  rownames(scores.filled)<-1:nrow(scores.filled)
  
  #Reverse time-embedding. Time-collapsing?
  scores.final<-NULL
  for(i in 1:nVars){
    
    scores.mean<-NULL
    for(j in 1:embedD){
      this.col<-scores.filled[,j]
      this.col<-c(rep(NA,j-1), this.col)
      this.col<-this.col[1:(length(this.col)-j+1)]
      scores.mean<-cbind(scores.mean, this.col)
    }
    
    intercept<-omxGetParameters(model$fit)[grep("intX",names(omxGetParameters(model$fit)))[i]]
    slope<-omxGetParameters(model$fit)[grep("slopeX",names(omxGetParameters(model$fit)))[i]]
    scores.mean<-rowMeans(scores.mean, na.rm=T)+intercept+(1:length(scores.mean))*slope
    scores.final<-cbind(scores.final, scores.mean)
    scores.final[model$exclusions]<-NA
  }
  return(list("scores"=scores.final,"shocks"=scores.shocks))
}