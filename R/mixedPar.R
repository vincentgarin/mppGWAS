############
# mixedPar #
############

# Sub-function from package magicQTL of Wei and Xu (2016)

# Reference:

# Wei, J., & Xu, S. (2016). A random-model approach to QTL mapping in
# multiparent advanced generation intercross (MAGIC) populations.
# Genetics, 202(2), 471-486.

# http://www.genetics.org/content/202/2/471.supplemental

mixedPar <-
function(dataframe,qq,optim.speed=TRUE){
#
   d<-dataframe
   y<-as.matrix(d[,1])
   x<-as.matrix(d[,-1])
   n<-nrow(d)
   s<-ncol(x)
#
   delta<-qq[[1]]
   uu<-qq[[2]]

##loglike() function
   loglike<-function(theta){
      lambda<-exp(theta)
      d0<-sum(log(abs(lambda*delta+1)))
      h<-1/(lambda*delta+1)
      yy<-cal.xx(x=yu,y=yu,H=diag(h))
      xy<-cal.xx(x=xu,y=yu,H=diag(h))
      xx<-cal.xx(x=xu,y=xu,H=diag(h))

      V<-abs(yy-t(xy)%*%solve(xx,tol=1e-50)%*%xy)
      d1<-unlist(determinant(xx))[1]
      loglike<--0.5*(d0+d1+(n-s)*log(V))
      return(-loglike)
   }

##
   yu<-t(uu)%*%y
   xu<-t(uu)%*%x
##solve parms
if(optim.speed){
   theta<-0
   parm<-optim(par=theta,fn=loglike,hessian = TRUE,method="L-BFGS-B",lower=-100,upper=100)
   lambda<-exp(parm$par)
   conv<-parm$convergence
   fn1<-parm$value
}else{
   intervals<-c(-100,-15,-10,-5,0,20,40,80)
   intials<-c(-16,-12,-7,-2,12,30,60)
   ni<-length(intervals)-1
   thetas<-NULL
   for ( i in 1:ni){
       theta<-intials[i]
       parms<-optim(par=theta,fn=loglike,hessian=TRUE,method="L-BFGS-B",lower=intervals[i],upper=intervals[i+1])
       theta<-parms$par
       parms<-optim(par=theta,fn=loglike,hessian=TRUE,method="L-BFGS-B",lower=-100,upper=100)
       fi<-parms$par
       thetas<-c(thetas,fi)
   }
   fns<-sapply(1:ni,function(i){
       theta<-thetas[i]
       fn<-loglike(theta)
       return(fn)
       } )
   ii<-which.min(fns)
   lambda<-exp(thetas[ii])
   fn1<-fns[ii]
}
   fn0<-loglike(-Inf)
   lrt<-2*(fn0-fn1)

##
   h<-1/(lambda*delta+1)
   xx<-cal.xx(x=xu,y=xu,H=diag(h))
   xy<-cal.xx(x=xu,y=yu,H=diag(h))
   beta<-solve(xx,xy,tol=1e-50)
   yt<-yu-xu%*%beta
   yPy<-sum(yt*h*yt)
   sigma2<-yPy/(n-s)
   v.beta<-solve(xx,tol=1e-50)*sigma2
   str.beta<-sqrt(diag(v.beta))
   tau<-lambda*sigma2
#
   random<-(delta*lambda)*h*yt
   rr<-list(beta=beta,v.beta=v.beta,random=random,se2=sigma2,lambda=lambda,sg2=tau)
   return(rr)
}
