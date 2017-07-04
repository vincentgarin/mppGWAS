##########
# cal.xx #
##########

# Sub-function from package magicQTL of Wei and Xu (2016)

# Reference:

# Wei, J., & Xu, S. (2016). A random-model approach to QTL mapping in
# multiparent advanced generation intercross (MAGIC) populations.
# Genetics, 202(2), 471-486.

# http://www.genetics.org/content/202/2/471.supplemental

cal.xx <-
function(x,y,H=NA){
   x<-as.matrix(x)
   y<-as.matrix(y)
   nr<-ncol(x)
   nc<-ncol(y)
   n0<-nrow(x)
#      xx<-matrix(0,nrow=nr,ncol=nc)
#      for ( i in 1:nr){
#         for (j in 1:nc){
#            xx[i,j]<-x[,i]*H*y[,j]
#         }
#      }
   if (is.na(H)) H<-diag(n0)
   xx<-t(x)%*%H%*%y
   xx<-as.matrix(xx)
   return(xx)
}
