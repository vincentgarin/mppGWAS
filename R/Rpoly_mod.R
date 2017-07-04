#############
# Rpoly_mod #
#############

# Modified sub-function from package magicQTL of Wei and Xu (2016)

# Reference:

# Wei, J., & Xu, S. (2016). A random-model approach to QTL mapping in
# multiparent advanced generation intercross (MAGIC) populations.
# Genetics, 202(2), 471-486.

# http://www.genetics.org/content/202/2/471.supplemental

Rpoly_mod <-
function(dataframe,gen,map,qq,r,cc,window,mixed,mixednew,rmfactor=FALSE){

   lower<--50
   upper<-50

##Declare the variable
###
   d<-dataframe
   y<-as.matrix(d[,1])
   x<-as.matrix(d[,-1])
   n<-nrow(y)
   s<-ncol(x)
   r<-r

   ######## add definition of a partition for unequal number of alleles

   r.part <- split(1:sum(r), rep(1:length(r), r))
#
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
   H<-1/(delta*lambda+1)

   ######## Only use models A so this part will never be called

   # if (rmfactor){
   #    blupp<-mixedBlup(dataframe=d,gen=gen,map=map,qq=qq,r=r,cc=cc,mixed=mixednew)
   #    effect<-blupp$gamma
   # }

##function -I, log-likelihood function
###
   loglike<-function(theta, r){
       xi<-exp(theta)
       tmp0<-xi*zz+diag(r)
       tmp<-solve(tmp0,tol=1e-90)
       yHy<-yy-xi*yz%*%tmp%*%t(yz)
       yHx<-yx-xi*yz%*%tmp%*%t(xz)
       xHx<-xx-xi*xz%*%tmp%*%t(xz)
       V<-abs(yHy-yHx%*%solve((1-1e-03)*xHx+1e-03*diag(rep(1,s)),tol=1e-90)%*%t(yHx))
       log.tmp0<-unlist(determinant(tmp0))[1]
       log.xx<-unlist(determinant(xHx))[1]
#loglike<- -0.5*log(det(diag(1/H))*det(tmp0))-0.5*(n-sn)*log(yHy-yHx%*%solve(xHx)%*%t(yHx)-0.5*log(det(xHx))
       loglike<- -0.5*(log.tmp0+(n-s)*log(V)+log.xx)
       if ( !is.finite(loglike)) loglike<--1e+10
       return(-loglike)
   }

##function -II, used to calculate the parameters, including beta, residual variacne
###Gamma and standard deviation
   fixed<-function(xi, r){
       tmp<-solve(xi*zz+diag(r))
       yHy<-yy-xi*yz%*%tmp%*%t(yz)
       yHx<-yx-xi*yz%*%tmp%*%t(xz)
       xHx<-xx-xi*xz%*%tmp%*%t(xz)
       zHy<-t(yz)-xi*zz%*%tmp%*%t(yz)
       zHx<-t(xz)-xi*zz%*%tmp%*%t(xz)
       zHz<-zz-xi*zz%*%tmp%*%zz
##
       Beta<-solve(xHx,t(yHx))# value of beta
       tmp<-solve(xHx)
       sigma2<-(yHy-yHx%*%tmp%*%t(yHx))/(n-s) # value of residual variance
       Gamma<-xi*zHy-xi*zHx%*%tmp%*%t(yHx) # value of random-effect markers
       Var<-abs((xi*diag(r)-xi^2*zHz)*as.numeric(sigma2)) # var(gamama/Y)=
       Stderr<-sqrt(diag(Var))# vector of variance of each allele
       RR<-list(Gamma=Gamma,Var=Var,Stderr=Stderr,Beta=Beta,sigma2=sigma2)
       return(RR)
   }

##compute the matrix,yy,yx,xx,
###
   y<-y-mean(y)
   yu<-t(uu)%*%y
   xu<-t(uu)%*%x

###loop by the marker
###
   parr<-numeric()

   # blupp<-numeric()

   for(k in 1:m){

     ############ a sub-function for a tryCatch structure could start here...

       # sub<-seq(((k-1)*r+1),((k-1)*r+r))
       # zu<-t(uu)%*%t(gen[sub,])

       zu<-t(uu)%*%t(gen[r.part[[k]], ])

##remove the overlap the snp effect

#        if (rmfactor){
#           c1<-ccM[k]-0.5*window
#           c3<-ccM[k]+0.5*window
#           k1<-max(which(ccM<=c1))
#           k3<-min(which(ccM>=c3))
#           if (k1==-Inf) k1<-1
#           if (k3==Inf)  k3<-m
# #          k1<-k-10
# #          k3<-k+10
# #          if (k1<1) k1<-1
# #          if (k3>m) k3<-m
#           remove<-k1:k3
#           cover<-length(remove)
#        }
#
       rm.factor<-matrix(0,n,1)
       offset.k<-rm.factor


       # if ( rmfactor){
       #     for(i in 1:cover){
       #         m0<-remove[i]
       #         co.sub<-seq(((m0-1)*r+1),((m0-1)*r+r))
       #         rm.factor<-rm.factor+t(gen[co.sub,])%*%as.matrix(effect[m0,],r,1)
       #     }
       #     offset.k<-t(uu)%*%rm.factor
       # }

       yk<-yu+offset.k

##
#xz<-t(xu)%*%diag(H)%*%zu; yz<-t(yu)%*%diag(H)%*%zu; zz<-t(zu)%*%diag(H)%*%zu;
       xz<-matrix(0,s,r[k])
       yz<-matrix(0,1,r[k])
       zz<-matrix(0,r[k],r[k])
       for(i in 1:s){
           for(j in 1:r[k]){
               xz[i,j]<-sum(xu[,i]*H*zu[,j])
           }
       }
       for(i in 1:r[k]){
           yz[i]<-sum(yk*H*zu[,i])
           for(j in 1:r[k]){
               zz[i,j]<-sum(zu[,i]*H*zu[,j])
           }
       }
#yy<-t(yu)%*%diag(H)%*%yu;yx<-t(yu)%*%diag(H)%*%xu;xx<-t(xu)%*%diag(H)%*%xu
       yy<-sum(yk*H*yk)
       yx<-matrix(0,1,s)
       xx<-matrix(0,s,s)
       for(i in 1:s){
          yx[i]<-sum(yk*H*xu[,i])
          for(j in 1:s){
             xx[i,j]<-sum(xu[,i]*H*xu[,j])
          }
       }


###
##
       theta<-0
       # parm<-optim(par=theta,fn=loglike,hessian = TRUE,method="L-BFGS-B",lower=lower,upper=upper)
       parm<-optim(par=theta,fn=loglike, r=r[k] ,hessian = TRUE,method="L-BFGS-B",
                   lower=lower,upper=upper)
       xi<-exp(parm$par)
       conv<-parm$convergence
       fn1<-parm$value
       fn0<-loglike(theta=c(-Inf), r = r[k])
       lrt<-2*(fn0-fn1)
       if(TRUE){
       if (lrt<=1e-03){
         lrt.p<-0.5+pchisq(lrt,df=1,lower.tail=FALSE)*0.5
       }else{
         lrt.p<-1-(0.5+pchisq(lrt,df=1)*0.5)
       }
       }
#       lrt.p<-pchisq(lrt,df=1,lower.tail=FALSE)
       lrt.logp<--log10(lrt.p)
       hess<-parm$hessian

#
       parmfix<-fixed(xi=xi, r = r[k])
       Gamma<-parmfix$Gamma
       Var<-parmfix$Var
       Stderr<-parmfix$Stderr
       Beta<-parmfix$Beta
       sigma2<-parmfix$sigma2
       tau_k<-xi*sigma2
       initial<-exp(theta)
#
       gg<-Gamma
#       wald<-sum(gg^2/Stderr^2)
       wald<-t(gg)%*%solve(Var)%*%gg
       wald.p<-pchisq(wald,df=r[k]-1,lower.tail=FALSE)
       wald.logp<--log10(wald.p)
#
       par<-data.frame(Num=k,chr=chr[k],ccM=ccM[k],lrt,lrt.p,lrt.logp,wald,
                       wald.p,wald.logp,tau_k,sigma2,lam_k=xi,conv)

       ############# could be the end of the sub-function evaluated in tryCatch
       ############### the function return par. if there is an error rep(0,...)

       # blup<-c(Gamma-mean(Gamma),Stderr)
       parr<-rbind(parr,par)
       # blupp<-rbind(blupp,blup)
   }

   cat("Data of chr have been completed",0,'\n')
   # colnames(blupp)<-c(paste0("Gamma",1:r),paste0("stderr",1:r))
   # result<-list(parr=parr,blupp=blupp)
   return(parr)
}
