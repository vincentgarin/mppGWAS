#############
# Fpoly_mod #
#############

# Modified sub-function from package magicQTL of Wei and Xu (2016)

# Reference:

# Wei, J., & Xu, S. (2016). A random-model approach to QTL mapping in
# multiparent advanced generation intercross (MAGIC) populations.
# Genetics, 202(2), 471-486.

# http://www.genetics.org/content/202/2/471.supplemental

Fpoly_mod <-
function(dataframe,gen,map,qq,r,cc,window,mixed,mixednew,rmfactor=FALSE){

##Declare the variable
###
   d<-dataframe
   y<-as.matrix(d[,1])
   x<-as.matrix(d[,-1])
   n<-nrow(y)
   s<-ncol(x)
   r<-r

   # partition to split gen object

   r.part <- split(1:sum(r), rep(1:length(r), r))

   gen<-gen
   map<-as.data.frame(map)
   m<-nrow(map)
   ccM<-map$cm
   chr<-map$chr
#
   qq<-qq
   delta<-qq[[1]]
   uu<-qq[[2]]
#
   mixed<-mixed
   lambda<-mixed$lambda
   s2.pre<-mixed$se2
   tau.pre<-lambda*s2.pre
   h<-1/(delta*lambda+1)
#
   ############# since I only use the A model, this part should not be used

   # if (rmfactor){
   #    blupp<-mixedBlup(dataframe=d,gen=gen,map=map,qq=qq,r=r,cc=cc,
   #                     mixed=mixednew)
   #    effect<-blupp$gamma
   # }

##loglike function - I
   loglike<-function(w,y){
       ns<-ncol(w)
       ns<-0
       ww<-cal.xx(x0=w,y0=w)
       wy<-cal.xx(x0=w,y0=y)
       yy<-cal.xx(x0=y,y0=y)
       b<-solve(ww,wy)

       log.sg2<-log(abs(yy-t(wy)%*%b)/(n-ns))
       log.ww<-unlist(determinant(ww))[1]
       log.v<-sum(log(abs(h)))
#       cc0<-(n-ns)*(1-log(n-ns))
       cc0<-n-ns
#       value<--0.5*((n-ns)*log.sg2+log.v+log.ww+cc0)
       value<--0.5*((n-ns)*log.sg2+log.v+cc0)
       return(value)
   }
##xx matrix calculate function - II
###
   cal.xx<-function(x0,y0){
       nr<-ncol(x0)
       nc<-ncol(y0)
       rr<-matrix(0,nr,nc)
       for ( i in 1:nr){
           for ( j in 1:nc){
               rr[i,j]<-sum(x0[,i]*h*y0[,j])
           }
       }
       return(rr)
   }
#

##loop by marker
###
#
   xu<-t(uu)%*%x
   yu<-t(uu)%*%y
   parr<-numeric()
   # blupp<-numeric()
   for(k in 1:m){
#
       # sub<-seq(((k-1)*r+1),((k-1)*r+r))
       # zu<-t(uu)%*%t(gen[sub,])

       zu<-t(uu)%*%t(gen[r.part[[k]], ])
#
##remove the overlap the snp effect
      if (rmfactor){
        c1<-ccM[k]-0.5*window
        c3<-ccM[k]+0.5*window
        k1<-max(which(ccM<=c1))
        k3<-min(which(ccM>=c3))
        if (k1==-Inf) k1<-1
        if (k3==Inf)  k3<-m
#       k1<-k-10
#       k3<-k+10
#       if (k1<1) k1<-1
#       if (k3>m) k3<-m
       remove<-k1:k3
       cover<-length(remove)
      }
#
       rm.factor<-matrix(0,n,1)
       offset.k<-rm.factor

       # again use only models A for the moment

       # if ( rmfactor){
       #     for(i in 1:cover){
       #         m0<-remove[i]
       #         # co.sub<-seq(((m0-1)*r+1),((m0-1)*r+r))
       #           co.sub<-r.part[[m0]]
       #         # rm.factor<-rm.factor+t(gen[co.sub,])%*%as.matrix(effect[m0,],r,1)
       #         rm.factor<-rm.factor+t(gen[co.sub,])%*%as.matrix(effect[m0,],r[m0],1)
       #
       #     }
       #     offset.k<-t(uu)%*%rm.factor
       # }

       yk<-yu+offset.k
       fn0<-loglike(w=xu,y=yk)
       yk<-yk-t(uu)%*%matrix(rep(mean(y),n),n,1)

##yy,zy and zz
       yy<-cal.xx(x0=yk,y0=yk)
       zy<-cal.xx(x0=zu,y0=yk)
       zz<-cal.xx(x0=zu,y0=zu)
       zzi<-solve(zz)
       b<-solve(zz,zy)
       s2<-(yy-t(zy)%*%b)/(n-r[k])
       v<-zzi*drop(s2)

#lrt
       fn1<-loglike(w=zu,y=yk)
       lrt<-2*(fn1-fn0)
       lrt.p<-pchisq(lrt,df=r[k]-s,lower.tail=FALSE)
       lrt.logp<--log10(lrt.p)
#wald
       g<-b
       gamma<-g
#       wald<-sum(gamma*(1/diag(v))*gamma)
       wald<-t(gamma)%*%solve(v)%*%gamma
       wald.p<-pchisq(wald,df=r[k]-1,lower.tail=FALSE,log.p=FALSE)
       wald.logp<--log10(wald.p)
       stderr<-sqrt(diag(v))

       par<-data.frame(Num=k,chr=chr[k],ccM=ccM[k],lrt,lrt.p,lrt.logp,wald,
                       wald.p,wald.logp,sigma2=s2)
       # blup<-c(gamma,stderr)
       parr<-rbind(parr,par)
       # blupp<-rbind(blupp,blup)
   }

   cat("Data of chr have been completed",0,'\n')
   # colnames(blupp)<-c(paste0("Gamma",1:r),paste0("stderr",1:r))
   # result<-list(parr=parr,blupp=blupp)
   return(parr)
}
