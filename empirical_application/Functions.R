bracket<-function(y,Delta,minb=-8, maxb=8,...) {
  yL<-y
  yU<-y
  brackets<-seq(minb,maxb,by=Delta)
  sample_size<-length(y)
  for (i in 1:sample_size) {
    for (j in 1:(length(brackets)-1)) {
      if (brackets[j] <=y[i]  && y[i]< brackets[j+1]) {
        yL[i]<-brackets[j]
        yU[i]<-brackets[j+1]
      }
    }
  }
  if (Delta==1) {
    yL[y>=floor(max(y))]<-floor(max(y))
    yU[y>=floor(max(y))]<-ceiling(max(y))
  }
  
  if (Delta==2) {
    yL[y>=5]<-5
    yU[y>=5]<-7
  }
  return(data.frame(yL=yL,yU=yU,y=y))
}