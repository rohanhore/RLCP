#---------------------------------------------------
#---------------simulation settings-----------------
#---------------------------------------------------
simulation=function(n,d,setting){
  X=rmvnorm(n,mean=rep(0,d),sigma=diag(d))
  
  ##setting 1
  if(setting==1){
    Y=0.5*apply(X,1,mean)+apply(abs(sin(X)),1,sum)*rnorm(n,0,1)
  }
  
  ##setting 2
  if(setting==2){
    Y=0.5*apply(X,1,mean)+apply(2*dnorm(X,0,1.5),1,sum)*rnorm(n,0,1)
  }
  
  ##setting 3 : P_X uniform on cube
  if(setting==3){
    X=matrix(runif(n*d,-3,3),nrow=n,ncol=d)
    Y=0.5*apply(X,1,mean)+apply(abs(sin(X)),1,sum)*rnorm(n,0,1)
  }
  
  data=as.data.frame(cbind(Y,X))
  colnames(data)=c("Y",paste0("X",1:d))
  return(data)
}
#--------------------------------------------------
#-------Simulating Y|X on a grid of ---------------
#--------feature points in (-3,3)-----------------
#--------------------------------------------------
conditional_simulation=function(n=100,setting){
  X=as.matrix(seq(-3,3,by=0.01))
  N=rnorm(length(X),0,1)
  
  ##setting 1
  if(setting==1){Y=0.5*apply(X,1,mean)+abs(sin(X))*N}
  ##setting 2
  if(setting==2){Y=0.5*apply(X,1,mean)+2*dnorm(X,0,1.5)*N}
  
  data=as.data.frame(cbind(Y,X));d=1
  colnames(data)=c("Y",paste0("X",1:d))
  return(data)
}
