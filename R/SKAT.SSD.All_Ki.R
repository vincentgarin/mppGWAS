########################
# SKAT.SSD.All (modif) #
########################

# Modification of the SKAT.SSD.All function to allow to subset the SSD.info
# object per chromosome and use is for the K_i variant of the Kernel association
# scan.

SKAT.SSD.All_Ki = function(SSD.INFO, obj, ..., obj.SNPWeight=NULL, verbose){

  N.Set<-SSD.INFO$nSets
  OUT.Pvalue<-rep(NA,N.Set)
  OUT.Marker<-rep(NA,N.Set)
  OUT.Marker.Test<-rep(NA,N.Set)
  OUT.Error<-rep(-1,N.Set)
  OUT.Pvalue.Resampling<-NULL

  Is.Resampling = FALSE
  n.Resampling = 0

  if(class(obj) == "SKAT_NULL_Model"){
    if(obj$n.Resampling > 0){
      Is.Resampling = TRUE
      n.Resampling = obj$n.Resampling

      OUT.Pvalue.Resampling<-matrix(rep(0,n.Resampling*N.Set),ncol=n.Resampling)
    }
  } else if(class(obj) == "SKAT_NULL_Model_ADJ"){
    if(obj$re1$n.Resampling > 0){
      Is.Resampling = TRUE
      n.Resampling = obj$re1$n.Resampling

      OUT.Pvalue.Resampling<-matrix(rep(0,n.Resampling*N.Set),ncol=n.Resampling)
    }
  }

  if(verbose){pb <- txtProgressBar(min=0, max=N.Set, style=3)}

  for(i in 1:N.Set){
    Is.Error<-TRUE
    try1 = try(SKAT::SKAT.SSD.OneSet_SetIndex(SSD.INFO=SSD.INFO,
                                              SetIndex=SSD.INFO$SetInfo$SetIndex[i],
                                              obj=obj, ...,
                                              obj.SNPWeight=obj.SNPWeight))

    if(class(try1) != "try-error"){
      re<-try1
      Is.Error<-FALSE
    } else {

      err.msg<-geterrmessage()
      msg<-sprintf("Error to run SKAT for %s: %s",SSD.INFO$SetInfo$SetID[i], err.msg)
      warning(msg,call.=FALSE)

    }

    if(!Is.Error){

      OUT.Pvalue[i]<-re$p.value
      OUT.Marker[i]<-re$param$n.marker
      OUT.Marker.Test[i]<-re$param$n.marker.test
      if(Is.Resampling){
        OUT.Pvalue.Resampling[i,]<-re$p.value.resampling
      }
    }

   if(verbose){setTxtProgressBar(pb, i)}

  }
  if(verbose){close(pb)}


  out.tbl<-data.frame(SetID=SSD.INFO$SetInfo$SetID, P.value=OUT.Pvalue,
                      N.Marker.All=OUT.Marker, N.Marker.Test=OUT.Marker.Test)
  re<-list(results=out.tbl,P.value.Resampling=OUT.Pvalue.Resampling)
  class(re)<-"SKAT_SSD_ALL"

  return(re)
}
