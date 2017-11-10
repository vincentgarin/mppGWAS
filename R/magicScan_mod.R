#################
# magicScan_mod #
#################

# modification of the function magicScan
#
# Modification of the function magicScan from Wei and Xu (2016) \code{magicQTL}
# package allowing to have marker with different number of alleles at each
# position.
#
# The only element that change is the argument \code{nfounders} that is
# replaced by \code{nalleles}. \code{nalleles} specify for each marker position
# the number of allele. The modified function does not compute the blup and
# return only the lrt and wald statistics.
#
# The \code{magicQTL} package can be found there:
# (\url{http://www.genetics.org/content/202/2/471.supplemental})
#
# @param dataframe input data frame. Row numbers is the number of observation,
# nobs. Column numbers is more than or equal to 2: the first column is the
# response variable; the second column is the intercept, 1 vector; the third...,
# is the other indicator variables, like sex(male or female), group.
#
# @param gen A list of probability matrix. Each element is the probability
# matrix of the chromosome: row numbers is m0*r, m0 denoting the markers in the
# chromosome, r for the founders; column numbers is the nobs, the number of the
# individuals. Each vetcor,r*1, is the probabilities of an allele derived from
# a certain founder, calculated by the R packages happy.hbrem or qtl.
#
# @param map A list of map. Each element is the data frame corresponding to the gen
# parameter: row numbers equals to the number of markers,m0; The first column,
# named markers, can be the marker name or the order of the marker in the
# chromosome; The second is the chromosome, named chr; The Third is the genetic
# map, named cm. Note, the column names must match the markers,chr and cm.
#
# @param kk.eigen A list, with 2 elements. The First element is the kinship matrix, named
# kk. The second element is the object of the R function eigen, named qq.
#
# @param nalleles \code{list} of vector: one element per chromosome, each vector
# element specifies the number of alleles for each markers.
#
# @param model A character, the model you would like to use to implement QTL mapping.
# Currently, included "Random-A","Fixed-A".
#
# @param index A vector, the orders of the response variables, which is used in
# the permutation. Default is NULL, meaning not permutation.
#
# @param window The size of the window, used in the models "Random-B","Fixed-B"
# and "CIM". Release the markers located in the 0.5*window of the left or right
# of the scanned marker.
#
# @param step A number, the space between markers selected as the cofators only
# in the model CIM.
#
# @param bychr Default 1. If the number of the cofactors is less than chromosome,
# need to decide which chromosome the cofactor is distributed in.
#
#
# @return Return:
#
# \item{parr}{A data frame that mainly contains the test statistics of each
# marker, LRT and wald.}
#
# @author Julong Wei and Shizhong Xu (original function), Vincent Garin (modification)
#
# @references
#
# Wei, J., & Xu, S. (2016). A random-model approach to QTL mapping in
# multiparent advanced generation intercross (MAGIC) populations.
# Genetics, 202(2), 471-486.
#
# @examples
#
#\dontrun{
#
# data(Ara)
# names(Ara)
# gen<-Ara[[1]]
# map<-Ara[[2]]
#
# # select only the two first chromosomes
#
# gen <- gen[c(1, 2)]
# map <- map[c(1, 2)]
# nalleles <- list()
#
# nalleles[[1]] <- rep(19, 14)
# nalleles[[2]] <- rep(19, 11)
#
# Ara.phe<-Ara[[3]]
# kk.eigen<-Ara[[4]]
# chrnum<-length(gen)
#
# indi<-nrow(Ara.phe)
# x<-rep(1,indi)
# y<-Ara.phe[,5] # Phenotype,total length, that is height of the Arabdopsis
# d<-data.frame(y=y,x=x)
#
# scan<-magicScan_mod(dataframe=d,gen=gen,map=map,kk.eigen=kk.eigen,
# nalleles=nalleles,model="Fixed-A") #scan the markers
#
# parms<-lapply(1:chrnum, function(i){ return(scans[[i]][[1]]) })
# parms<-do.call(rbind,parms)
#
#}
#

# could add @export there


magicScan_mod <-
function(dataframe,gen,map,kk.eigen,nalleles,model="Fixed-A",index=NULL,
         window=5,step,bychr=1,verbose=FALSE){

  d<-dataframe
  kk<-kk.eigen[[1]]
  qq<-kk.eigen[[2]]
  cc<-3596
  chrnum<-length(gen)

  genold<-gen
  if (!is.null(index)){
    gennew<-list()
    for(i in 1:chrnum){
      gennew[[i]]<-as.matrix(gen[[i]][,index])
    }
    gen<-gennew
  }

  RR<-mixedPar(dataframe=d,qq=qq,optim.speed=FALSE)

  if(verbose){

    cat("lambda:",RR$lambda,"Residual error:",RR$se2,"\n")
    cat("Model:",model,"\n")

  }

###


##function-1,Random-A,call Rpoly
  if (model=="Random-A"){
    scans<-list()
    for ( i in 1:chrnum){
      scans[[i]]<-Rpoly_mod(dataframe=d,gen=gen[[i]],map=map[[i]],qq=qq,
                        r=nalleles[[i]],mixed=RR,rmfactor=FALSE,
                        verbose=verbose)
    }
  }


  ##function-3,Fixed-A, call Fpoly
  if (model=="Fixed-A"){
    scans<-list()
    for ( i in 1:chrnum){
      scans[[i]]<-Fpoly_mod(dataframe=d,gen=gen[[i]],map=map[[i]],qq=qq,
                        r=nalleles[[i]],mixed=RR,rmfactor=FALSE,
                        verbose=verbose)
    }
  }

  return(scans)

# I will only use Fixed-A and Random-A model for the moment

# if (model=="Random-B"|model=="Fixed-B"){
#   if(!is.null(index)){
#     kknew<-kk[index,index]
#     qqnew<-eigen(kknew,symmetric=TRUE)
#     RRnew<-mixedPar(dataframe=d,qq=qqnew,optim.speed=FALSE)
#   }else{
#     RRnew<-RR
#   }
# }


# ##function-2,Random-B, call Rpoly
#   if (model=="Random-B"){
#     scans<-list()
#     for ( i in 1:chrnum){
#       scans[[i]]<-Rpoly(dataframe=d,gen=gen[[i]],map=map[[i]],qq=qq,
#                         r=nalleles[[i]],cc=cc,window=window,mixed=RR,
#                         mixednew=RRnew,rmfactor=TRUE)
#     }
#   }


# ##function-4, Fixed-B, call Fpoly
#   if (model=="Fixed-B"){
#     scans<-list()
#     for ( i in 1:chrnum){
#       scans[[i]]<-Fpoly(dataframe=d,gen=gen[[i]],map=map[[i]],qq=qq,
#                         r=nalleles[[i]], cc=cc,window=window,mixed=RR,
#                         mixednew=RRnew,rmfactor=TRUE)
#     }
#   }


# ##function-5,Interval mapping, call Fixed
#   if (model=="IM"){
#     scans<-list()
#     for ( i in 1:chrnum){
#       scans[[i]]<-Fixed(dataframe=d,gen=gen[[i]],map=map[[i]],r=nalleles[[i]])
#     }
#   }


# ##function-6,composite interval mapping, call Fixedcim
#   if (model=="CIM"){
#
#     map0<-as.data.frame(do.call(rbind,map))
#     gen0<-as.matrix(do.call(rbind,gen))
#     lasso<-effect.fixed(dataframe=d,gen=gen0,map=map0,r=nalleles,step=step,
#                         bychr=bychr)
#     cat("The number of cofatcors:",length(lasso$covar),"\n")
#     cat("cofactor have been completed:",0,"\n")
#
#     scans<-list()
#     for ( i in 1:chrnum){
#       rr<-Fixedcim(dataframe=d,gen=gen[[i]],map=map[[i]],r=nalleles,
#                    window=window,lasso=lasso)
#       scans[[i]]<-rr
#     }
#   }

     # return(scans)

}
