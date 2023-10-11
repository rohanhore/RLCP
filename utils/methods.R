#--------------------------------------------------
#--------generating X_tilde for box kernel---------
#--------------------------------------------------
runifball=function(n,center,radius){
  d=length(center)
  data=matrix(0,nrow=n,ncol=d)
  
  U=runif(n,min=0,max=1)
  Z=matrix(rnorm(n*d),nrow=n,ncol=d)
  for (i in 1:n){
    data[i,]=center+radius*U[i]^(1/d)*Z[i,]/sqrt(sum(Z[i,]^2))}
  return(data)
}

#--------------------------------------------------
#--------localization kernel for box kernel--------
#--------------------------------------------------
euclid_distance=function(x,y){return(sqrt(sum((x-y)^2)))}

#--------------------------------------------------
#--------computing vanilla weighted quantile-------
#--------------------------------------------------
weighted_quantile=function(x,q,w){
  w=w/sum(w)
  ordering=order(x)
  emp_cdf=cumsum(w[ordering])
  quantile=x[ordering][min(which(emp_cdf>=q))]
  return(quantile)
}

#---------------------------------------------------
#----computing smoothed weighted quantile-----------
#---------------------------------------------------
##----here v consists of the unique scores, while indices contain the data-indices that
##----has led to same score, w is the weight vector.
smoothed_weighted_quantile=function(v,alpha,w,indices){
  w=w/sum(w)
  U=runif(1,min=0,max=1)
  
  #finding the weights corresponding to each unique score.
  v_tilde=v
  w_tilde=rep(0,length(v))
  for(i in 1:length(v)){
    w_tilde[i]=sum(w[indices[[i]]])
  }
  
  #computing p-value at points in between the calibration scores.
  p_values=rep(0,length(v_tilde))
  for(i in 1:length(v_tilde)-1){
    p_values[i]=sum(w_tilde[i:(length(v_tilde)-1)])+U*(tail(w_tilde,1))
  }
  #computing p-value at a point higher than all calibration scores.
  p_values[length(v_tilde)]=U*(tail(w_tilde,1))
  
  #if pvalue is never greater than alpha, we output the empty set.
  if(sum((p_values>alpha))>0){
    id=max(which(p_values>alpha))
    #now we check, whether the prediction interval will be a closed interval or open.
    quantile=v_tilde[id]
    if(id<length(v_tilde)-1){closed=(sum(w_tilde[(id+1):(length(v_tilde)-1)])+U*(w_tilde[id]+tail(w_tilde,1))>alpha)}
    if(id==length(v_tilde)){closed=FALSE}
    if(id==length(v_tilde)-1){closed=(U*(w_tilde[id]+tail(w_tilde,1))>alpha)}}
  else{quantile=-Inf;closed=FALSE}
  return(c(quantile,closed))
}

#---------------------------------------------------
#---------optimum bandwidth for calLCP with---------
#---------------prefixed effective size-------------
#---------------------------------------------------
##----given a pre-training dataset 'Xtrain', it aims to find the optimum bandwidth 
##----so that the effective sample size for calLCP, matches some prerfixed quantity 'eff_size'
optimum_calLCP_h=function(Xtrain,kernel,h_min,eff_size){
  ntrain=dim(Xtrain)[1];d=dim(Xtrain)[2]
  
  #estimating the effective sample size of the kernel at bandwidth 'h' for data 'Xtrain'.
  eff.size=function(h){
    H=matrix(0,nrow=ntrain,ncol=ntrain)
    for(i in 1:ntrain){
      if(kernel=="gaussian"){H[i,]=dmvnorm(Xtrain,mean=Xtrain[i,],sigma=diag(d)*h^2)}
      if(kernel=="box"){H[i,]=apply(Xtrain,1,FUN=function(x){(euclid_distance(x,Xtrain[i,])<=h)+0})}
      H[i,]=H[i,]/sum(H[i,])
    }
    eff.size=ntrain/norm(H,type="F")^2-1
    return(eff.size)
  }
  
  #candidate bandwidth choice grid, starting from the user-define minimum value.
  candidate_bandwidths=seq(h_min,6,by=0.02)
  
  #finding optimum bandwidth choice
  i=1;optimizer=0
  while((optimizer<eff_size)|(i>length(candidate_bandwidths))){
    optimizer=eff.size(candidate_bandwidths[i])
    h_opt=candidate_bandwidths[i]
    i=i+1
  }
  return(h_opt)
}

#---------------------------------------------------
#---------optimum bandwidth for RLCP with-----------
#---------------prefixed effective size-------------
#---------------------------------------------------
##----given a pre-training dataset 'Xtrain', it aims to find the optimum bandwidth 
##----so that the effective sample size for RLCP, matches some prerfixed quantity 'eff_size'
optimum_RLCP_h=function(Xtrain,kernel,h_min,eff_size){
  ntrain=dim(Xtrain)[1];d=dim(Xtrain)[2]
  
  #estimating the effective sample size of the kernel at bandwidth 'h' for data 'Xtrain'.
  eff.size=function(h){
    H=matrix(0,nrow=ntrain,ncol=ntrain)
    for(i in 1:ntrain){
      if(kernel=="gaussian"){
        xtilde_train=rmvnorm(1,mean=Xtrain[i,],sigma=diag(d)*h^2)
        H[i,]=dmvnorm(Xtrain,mean=xtilde_train,sigma=diag(d)*h^2)
      }
      if(kernel=="box"){
        xtilde_train=runifball(1,Xtrain[i,],h)
        H[i,]=apply(Xtrain,1,FUN=function(x){(euclid_distance(x,xtilde_train)<=h)+0})
      }
      H[i,]=H[i,]/sum(H[i,])
    }
    eff.size=ntrain/norm(H,type="F")^2-1
    return(eff.size)
  }
  
  #candidate bandwidth choice grid, starting from the user-define minimum value.
  candidate_bandwidths=seq(h_min,6,by=0.02)
  
  #finding optimum bandwidth choice
  i=1;optimizer=0
  while((optimizer<eff_size)|(i>length(candidate_bandwidths))){
    optimizer=eff.size(candidate_bandwidths[i])
    h_opt=candidate_bandwidths[i]
    i=i+1
  }
  return(h_opt)
}

####--------------------------------------------------------------------#####
####---------------Locally Weighted CP Methods--------------------------#####
####--------------------------------------------------------------------#####

#--------------------------------------------------
#-----------------------RLCP-----------------------
#--------------------------------------------------
RLCP=function(Xcalib,scores_calib,Xtest,scores_test,kernel,h,alpha){
  ntest=dim(Xtest)[1];d=dim(Xtest)[2]
  coverage=score_threshold=rep(0,ntest)
  
  #sorting with respect to the order of calibration scores.
  Xcalib=as.matrix(Xcalib[order(scores_calib),]);scores_calib=sort(scores_calib)
  
  #finding unique scores and the indices where each of these unique scores have been repeated.
  scores=c(scores_calib,Inf)
  indices=list();j=1;i=1
  scores_unique=vector()
  while(i<=length(scores)){
    scores_unique=c(scores_unique,scores[i])
    indices[[j]]=which(scores==scores[i])
    i=i+length(indices[[j]]);j=j+1
  }
  
  for(i in 1:ntest){
    xtest=Xtest[i,];test_score=scores_test[i]
    cov_data=rbind(Xcalib,xtest)
    
    #finding the weights and the score threshold
    if(kernel=="gaussian"){
      xtilde_test=rmvnorm(1,mean=xtest,sigma=diag(d)*h^2)
      weights=dmvnorm(cov_data,mean=xtilde_test,sigma=diag(d)*h^2)
      result=smoothed_weighted_quantile(scores_unique,alpha,weights,indices)
    }
    if(kernel=="box"){
      xtilde_test=runifball(1,xtest,h)
      weights=apply(cov_data,1,FUN=function(x){(euclid_distance(x,xtilde_test)<=h)+0})
      result=smoothed_weighted_quantile(scores_unique,alpha,weights,indices)
    }
    
    score_threshold[i]=result[1] #score_threshold
    closed=result[2]   #whether it's a closed interval
    
    #coverage
    coverage[i]=(test_score<score_threshold[i])+0
    if(closed==TRUE){coverage[i]=(test_score<=score_threshold[i])+0}
  }
  return(cbind(coverage,score_threshold))
}

#---------------------------------------------------
#---------------------baseLCP-----------------------
#---------------------------------------------------
baseLCP=function(Xcalib,scores_calib,Xtest,scores_test,kernel,h,alpha){
  ntest=dim(Xtest)[1];d=dim(Xtest)[2]
  coverage=score_threshold=rep(0,ntest)
  
  #sorting with respect to the order of calibration scores.
  Xcalib=as.matrix(Xcalib[order(scores_calib),]);scores_calib=sort(scores_calib)
  
  #finding unique scores and the indices where each of these unique scores have been repeated.
  scores=c(scores_calib,Inf)
  indices=list();j=1;i=1
  scores_unique=vector()
  while(i<=length(scores)){
    scores_unique=c(scores_unique,scores[i])
    indices[[j]]=which(scores==scores[i])
    i=i+sum(scores==scores[i]);j=j+1
  }
  
  for(i in 1:ntest){
    xtest=Xtest[i,];test_score=scores_test[i]
    cov_data=rbind(Xcalib,xtest)
    
    #finding the weights and the score threshold
    if(kernel=="gaussian"){
      weights=dmvnorm(cov_data,mean=xtest,sigma=diag(d)*h^2)
      result=smoothed_weighted_quantile(scores_unique,alpha,weights,indices)
    }
    if(kernel=="box"){
      weights=apply(cov_data,1,FUN=function(x){(euclid_distance(x,xtest)<=h)+0})
      result=smoothed_weighted_quantile(scores_unique,alpha,weights,indices)
    }
    score_threshold[i]=result[1] #score_threshold
    closed=result[2]   #whether it's a closed interval
    
    #coverage
    coverage[i]=(test_score<score_threshold[i])+0
    if(closed==TRUE){coverage[i]=(test_score<=score_threshold[i])+0}
  }
  return(cbind(coverage,score_threshold))
}

#---------------------------------------------------
#---------------------calLCP------------------------
#---------------------------------------------------
calLCP=function(Xcalib,scores_calib,Xtest,scores_test,kernel,h,alpha){
  ntest=dim(Xtest)[1];d=dim(Xtest)[2];ncalib=dim(Xcalib)[1]
  coverage=score_threshold=rep(0,ntest)
  
  #sorting with respect to the order of calibration scores.
  Xcalib=as.matrix(Xcalib[order(scores_calib),])
  scores_calib=(sort(scores_calib))
  
  #finding the un-normalized weight matrix on the calibration set.
  H=matrix(0,nrow=ncalib,ncol=ncalib)
  for(i in 1:ncalib){
    if(kernel=="gaussian"){H[i,]=dmvnorm(Xcalib,mean=Xcalib[i,],sigma=diag(d)*h^2)}
    if(kernel=="box"){H[i,]=apply(Xcalib,1,FUN=function(x){(euclid_distance(x,Xcalib[i,])<=h)+0})}
  }
  
  #denominator of transformed score if the original test score is Infinite.
  denom_calib=apply(H,1,sum)
  #numerator of transformed score if the original test score is Infinite.
  num_calib=rep(0,ncalib)
  for(i in 1:ncalib){num_calib[i]=sum(H[i,]*(scores_calib<scores_calib[i]))}

  for(i in 1:ntest){
    U=runif(1,0,1);p_values=rep(0,ncalib+1)
    
    #defining smoothed p-value for calLCP in terms of transformed scores.
    smoothed_p_value=function(x){return(sum(x>tail(x,1))/length(x)+(U*sum(x==tail(x,1)))/length(x))}
    
    xtest=Xtest[i,];test_score=scores_test[i]
    cov_data=rbind(Xcalib,xtest)
    scores=c(scores_calib,Inf)
    
    #weights of the test point on itself and all the calibration points.
    if(kernel=="gaussian"){weights=dmvnorm(cov_data,mean=xtest,sigma=diag(d)*h^2)}
    if(kernel=="box"){weights=apply(cov_data,1,FUN=function(x){(euclid_distance(x,xtest)<=h)+0})}
    
    #computing the transformed scores when v is smaller than all calibration points.
    T_values=rep(0,ncalib+1)
    for(j in 1:ncalib){
      T_values[j]=(num_calib[j]+weights[j])/(denom_calib[j]+weights[j])
    }
    T_values[ncalib+1]=sum(weights*(scores<scores[1]))/sum(weights)
    #the corresponding smoothed p-value.
    p_values[1]=smoothed_p_value(T_values)
    
    #computing the transformed scores when v is in between two calibration points.
    for(j in 2:(ncalib+1)){
      #making suitable changes to transformed score T
      T_values[ncalib+1]=sum(weights*(scores<scores[j]))/sum(weights)
      T_values[j-1]=(num_calib[j-1])/(denom_calib[j-1]+weights[j-1])
      #corresponding p-value
      p_values[j]=smoothed_p_value(T_values)
    }
    
    #if pvalue is never greater than alpha, we output the empty set.
    if(sum((p_values>alpha))>0){
      id=max(which(p_values>alpha))
      #checking whether it's a closed interval or not.
      if(id==ncalib+1){closed=FALSE}
      if(id<=ncalib){
        T_values[ncalib+1]=sum(weights*(scores<scores[id]))/sum(weights)
        for(j in 1:id){T_values[j]=(num_calib[j])/(denom_calib[j]+weights[j])}
        if(id<ncalib){
          for(j in (id+1):ncalib){T_values[j]=(num_calib[j])/(denom_calib[j]+weights[j])}}
        closed=(smoothed_p_value(T_values)>alpha)
      }
      score_threshold[i]=scores[id]
      
      #coverage
      if(closed==TRUE){coverage[i]=(test_score<=score_threshold[i])+0}
      else{coverage[i]=(test_score<score_threshold[i])+0}
      }
    else{coverage[i]=0;score_threshold[i]=-Inf}
  }
  return(cbind(coverage,score_threshold))
}

#---------------------------------------------------
#-----------------------mRLCP-----------------------
#---------------------------------------------------
##-----finding coverage for mrLCP for m=10,20,30,.. largest multiple of 10 before m.
mRLCP=function(Xcalib,scores_calib,Xtest,scores_test,kernel,h,alpha,m){
  ntest=dim(Xtest)[1];d=dim(Xtest)[2]
  coverage=matrix(0,nrow=ntest,ncol=floor(m/10))
  
  for(i in 1:ntest){
    xtest=Xtest[i,];test_score=scores_test[i]
    cov_data=rbind(Xcalib,xtest)
    scores=c(scores_calib,test_score)
    
    #finding p values for m runs of RLCP.
    pval=rep(0,m)
    for(k in 1:m){
      if(kernel=="gaussian"){
        xtilde_test=rmvnorm(1,mean=xtest,sigma=diag(d)*h^2)
        weights=dmvnorm(cov_data,mean=xtilde_test,sigma=diag(d)*h^2)
        weights=weights/sum(weights)
        pval[k]=sum(weights*(scores>test_score))+sum(weights*(scores==test_score))*runif(1)
      }
      if(kernel=="box"){
        xtilde_test=runifball(1,xtest,h)
        weights=apply(cov_data,1,FUN=function(x){(euclid_distance(x,xtilde_test)<=h)+0})
        weights=weights/sum(weights)
        pval[k]=sum(weights*(scores>test_score))+sum(weights*(scores==test_score))*runif(1)
      }
    }
    
    #p-value aggregation.
    p_value=cumsum(apply(matrix(pval,ncol=10),1,mean))/(1:(m/10))
    coverage[i,]=(p_value>alpha)+0
  }
  return(coverage)
}
