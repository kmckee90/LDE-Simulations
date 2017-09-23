
stdz<-function(var){
  return((var-mean(var,na.rm=T))/sd(var,na.rm=T))
}

insertRow<-function(data, vals, rowNum){
  data<-as.matrix(data)
  for(i in rowNum){
    if(i<nrow(data) && i>1){ 
      above<-data[1:(i-1),] 
      below<-data[i:nrow(data),]
      data<-rbind(above, vals, below)
    }else if(i==1){
      data<-rbind(vals, data)
    }else{
      data<-rbind(data, vals)
    }
  }
  return(data)
}