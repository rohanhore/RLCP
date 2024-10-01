#---------------------------------------------------
#-------------------source code---------------------
#---------------------------------------------------
source("../utils/simu_sett.R")
source("../utils/methods.R")

#---------------------------------------------------
#------------------load libraries-------------------
#---------------------------------------------------
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(mvtnorm))
suppressPackageStartupMessages(library(randomForest))
suppressPackageStartupMessages(library(neuralnet))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(ggplot2))
#--------------------------------------------------------
#--------data loading & pre-processing-------------------
#--------------------------------------------------------
data=read.table("abalone.data",sep=",")
data=as.data.frame(data)
data=data[,-c(6,7,8)]
colnames(data)=c("sex","length","diameter","height","whole_weight","rings")
data$sex=as.factor(data$sex)
data$sex=dplyr::recode_factor(data$sex,"M" = "1","F" = "0", "I"="2")
n=dim(data)[1]

# scaling data for neural nets
scaled_data=data
scaled_data[,-1]=as.matrix(sweep(data[,-1],2,apply(data[,-1],2,min))) %*% 
  diag(1/(apply(as.matrix(data[,-1]),2,max)-apply(as.matrix(data[,-1]),2,min)))
max_rings=max(data$rings);min_rings=min(data$rings)
scaled_data=as.data.frame(scaled_data)
colnames(scaled_data)=c("sex","length","diameter","height","whole_weight","rings")
dmy=dummyVars(" ~ .", data = scaled_data)
scaled_data <- data.frame(predict(dmy, newdata = scaled_data))


###-----computing local coverage across diameter & length bins-----------
coverage_diam=function(coverage,test_data){
  quantile_points=seq(0.125,0.875,by=0.01)
  coverage_diam=rep(0,length(quantile_points))
  for(i in 1:length(quantile_points)){
    id=which((test_data$diameter>=quantile(data$diameter,probs=quantile_points[i]-0.025)) & 
               (test_data$diameter<=quantile(data$diameter,probs=quantile_points[i]+0.025)))
    coverage_diam[i]=mean(coverage[id])
  }
  return(coverage_diam)
}

coverage_length=function(coverage,test_data){
  quantile_points=seq(0.125,0.875,by=0.01)
  coverage_length=rep(0,length(quantile_points))
  for(i in 1:length(quantile_points)){
    id=which((test_data$length>=quantile(data$length,probs=quantile_points[i]-0.025)) & 
               (test_data$length<=quantile(data$length,probs=quantile_points[i]+0.025)))
    coverage_length[i]=mean(coverage[id])
  }
  return(coverage_length)
}

###-----local coverage across sex groups-----------
coverage_sex=function(coverage,test_data){
  test_data$sex=as.factor(test_data$sex)
  coverage_sex=rep(0,3)
  for(i in 1:3){
    id=which(test_data$sex==c(0,1,2)[i])
    coverage_sex[i]=mean(coverage[id])
  }
  return(coverage_sex)
}

#----------------------------------------------------------
#-----------------RLCP for real data-----------------------
#----------------------------------------------------------

RLCP_real=function(Xcalib,scores_calib,Xtest,scores_test,h,alpha){
  ntest=dim(Xtest)[1];d=dim(Xtest)[2]
  coverage=threshold=rep(0,ntest)
  Xcalib=as.matrix(Xcalib[order(scores_calib),]);scores_calib=sort(scores_calib)
  
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
    xtilde_test=xtest
    xtilde_test[2:5]=xtest[2:5]+runif(4,min=-h,max=h)
    
    cov_data=rbind(Xcalib,xtest)
    
    weights=apply(abs(sweep(cov_data,2,as.numeric(xtilde_test),"-")),1,FUN=function(x) all(x<=c(0,rep(h,4)))+0)
    result=smoothed_weighted_quantile(scores_unique,alpha,weights,indices)
    
    threshold[i]=result[1]
    closed=result[2]
    coverage[i]=(test_score<threshold[i])+0
    if(closed==TRUE){coverage[i]=(test_score<=threshold[i])+0}
  }
  return(cbind(coverage,threshold))
}

#-----------------------------------------------------------
#-----------------baseLCP for real data---------------------
#-----------------------------------------------------------

baseLCP_real=function(Xcalib,scores_calib,Xtest,scores_test,h,alpha){    
  ntest=dim(Xtest)[1];d=dim(Xtest)[2]
  coverage=threshold=rep(0,ntest)
  Xcalib=as.matrix(Xcalib[order(scores_calib),]);scores_calib=sort(scores_calib)
  
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
    
    weights=apply(sweep(cov_data,2,as.numeric(xtest)),1,FUN=function(x) all(x<=c(0,rep(h,4)))+0)
    result=smoothed_weighted_quantile(scores_unique,alpha,weights,indices)
    threshold[i]=result[1]
    closed=result[2]
    coverage[i]=(test_score<threshold[i])+0
    if(closed==TRUE){coverage[i]=(test_score<=threshold[i])+0}
  }
  return(cbind(coverage,threshold))
}

#------------------------------------------------------------
#-----------------calLCP for real data-----------------------
#------------------------------------------------------------

calLCP_real=function(Xcalib,scores_calib,Xtest,scores_test,h,alpha){
  ntest=dim(Xtest)[1];d=dim(Xtest)[2];ncalib=dim(Xcalib)[1]
  coverage=threshold=rep(0,ntest)
  Xcalib=as.matrix(Xcalib[order(scores_calib),])
  scores_calib=(sort(scores_calib))
  
  H=matrix(0,nrow=ncalib,ncol=ncalib)
  for(i in 1:ncalib){
    H[i,]=apply(sweep(Xcalib,2,as.numeric(Xcalib[i,])),1,FUN=function(x) all(x<=c(0,rep(h,4)))+0)
  }
  denom_calib=apply(H,1,sum)
  num_calib=rep(0,ncalib)
  for(i in 1:ncalib){num_calib[i]=sum(H[i,]*(scores_calib<scores_calib[i]))}
  
  p_values=rep(0,ncalib+1)
  for(i in 1:ntest){
    U=runif(1,0,1)
    
    smoothed_p_value=function(x){
      return(sum(x>tail(x,1))/length(x)+(U*sum(x==tail(x,1)))/length(x))
    }
    
    xtest=Xtest[i,];test_score=scores_test[i]
    cov_data=rbind(Xcalib,xtest)
    scores=c(scores_calib,Inf)
    
    weights=apply(sweep(cov_data,2,as.numeric(xtest)),1,FUN=function(x) all(x<=c(0,rep(h,4)))+0)
    
    T_values=rep(0,ncalib+1)
    for(j in 1:ncalib){
      T_values[j]=(num_calib[j]+weights[j])/(denom_calib[j]+weights[j])
    }
    T_values[ncalib+1]=sum(weights*(scores<scores[1]))/sum(weights)
    
    p_values[1]=smoothed_p_value(T_values)
    for(j in 2:(ncalib+1)){
      T_values[ncalib+1]=sum(weights*(scores<scores[j]))/sum(weights)
      T_values[j-1]=(num_calib[j-1])/(denom_calib[j-1]+weights[j-1])
      p_values[j]=smoothed_p_value(T_values)
    }
    
    #if pvalue is never greater than alpha, we output the empty set.
    if(sum((p_values>alpha))>0){
      id=max(which(p_values>alpha));closed=FALSE
      if(id<=ncalib){
        T_values[ncalib+1]=sum(weights*(scores<scores[id]))/sum(weights)
        for(j in 1:id){T_values[j]=(num_calib[j])/(denom_calib[j]+weights[j])}
        if(id<ncalib){
          for(j in (id+1):ncalib){T_values[j]=(num_calib[j])/(denom_calib[j]+weights[j])}}
        closed=(smoothed_p_value(T_values)>alpha)
      }
      threshold[i]=scores[id]
      
      if(closed==TRUE){coverage[i]=(test_score<=threshold[i])+0}
      else{coverage[i]=(test_score<threshold[i])+0}
    }
    else{coverage[i]=0;score_threshold[i]=-Inf}
  }
  return(cbind(coverage,threshold))
}

#-------------------------------------------------------------------------------
#--------------analyzing marginal,local coverages for real data-----------------
#-------------------------------------------------------------------------------
real_analysis=function(h,k){
  split=sample(1:3,n,replace=TRUE,prob=c(1/3,1/3,1/3))
  train_data=data[split==1,];ntrain=dim(train_data)[1]
  calib_data=data[split==2,];ncalib=dim(calib_data)[1]
  test_data=data[split==3,];ntest=dim(test_data)[1]
  
  train_scaled_data=scaled_data[split==1,]
  calib_scaled_data=scaled_data[split==2,]
  test_scaled_data=scaled_data[split==3,]
  
  #------------learning score on train split-------------------------
  Xcalib=calib_data[,1:5]
  Xcalib[,1]=as.numeric(levels(calib_data[,1]))[calib_data[,1]]
  
  Xtest=test_data[,1:5]
  Xtest[,1]=as.numeric(levels(test_data[,1]))[test_data[,1]]
  
  #---linear model---------
  model_lm=lm(rings~.,data=train_data)
  predict_lm_test=predict(model_lm,newdata=test_data)
  scores_lm_calib=abs(calib_data$rings-predict.lm(model_lm,calib_data))
  scores_lm_test=abs(test_data$rings-predict.lm(model_lm,test_data))
  
  result_lm_calLCP=calLCP_real(Xcalib, scores_lm_calib,Xtest,scores_lm_test,h,alpha)
  result_lm_baseLCP=baseLCP_real(Xcalib, scores_lm_calib,Xtest,scores_lm_test,h,alpha)
  result_lm_RLCP=RLCP_real(Xcalib, scores_lm_calib,Xtest,scores_lm_test,h,alpha)
  
  coverage_lm_RLCP=mean(result_lm_RLCP[,1])
  coverage_lm_calLCP=mean(result_lm_calLCP[,1])
  coverage_lm_baseLCP=mean(result_lm_baseLCP[,1])
  
  width_lm_RLCP=2*median(result_lm_RLCP[,2])
  width_lm_calLCP=2*median(result_lm_calLCP[,2])
  width_lm_baseLCP=2*median(result_lm_baseLCP[,2])
  
  coverage_diam_lm_RLCP=coverage_diam(result_lm_RLCP[,1],Xtest)
  coverage_diam_lm_calLCP=coverage_diam(result_lm_calLCP[,1],Xtest)
  coverage_diam_lm_baseLCP=coverage_diam(result_lm_baseLCP[,1],Xtest)
  
  coverage_length_lm_RLCP=coverage_length(result_lm_RLCP[,1],Xtest)
  coverage_length_lm_calLCP=coverage_length(result_lm_calLCP[,1],Xtest)
  coverage_length_lm_baseLCP=coverage_length(result_lm_baseLCP[,1],Xtest)

  coverage_sex_lm_RLCP=coverage_sex(result_lm_RLCP,Xtest)
  coverage_sex_lm_calLCP=coverage_sex(result_lm_calLCP,Xtest)
  coverage_sex_lm_baseLCP=coverage_sex(result_lm_baseLCP,Xtest)

  
  #----random forest----------
  model_rf=randomForest(rings ~ .,data=train_data)
  predict_rf_test=predict(model_rf,newdata=test_data)
  scores_rf_calib=abs(calib_data$rings-predict(model_rf,calib_data))
  scores_rf_test=abs(test_data$rings-predict(model_rf,test_data))
  
  result_rf_calLCP=calLCP_real(Xcalib, scores_rf_calib,Xtest,scores_rf_test,h,alpha)
  result_rf_baseLCP=baseLCP_real(Xcalib, scores_rf_calib,Xtest,scores_rf_test,h,alpha)
  result_rf_RLCP=RLCP_real(Xcalib, scores_rf_calib,Xtest,scores_rf_test,h,alpha)
  
  coverage_rf_RLCP=mean(result_rf_RLCP[,1])
  coverage_rf_calLCP=mean(result_rf_calLCP[,1])
  coverage_rf_baseLCP=mean(result_rf_baseLCP[,1])
  
  width_rf_RLCP=2*median(result_rf_RLCP[,2])
  width_rf_calLCP=2*median(result_rf_calLCP[,2])
  width_rf_baseLCP=2*median(result_rf_baseLCP[,2])
  
  coverage_diam_rf_RLCP=coverage_diam(result_rf_RLCP[,1],Xtest)
  coverage_diam_rf_calLCP=coverage_diam(result_rf_calLCP[,1],Xtest)
  coverage_diam_rf_baseLCP=coverage_diam(result_rf_baseLCP[,1],Xtest)
  
  coverage_length_rf_RLCP=coverage_length(result_rf_RLCP[,1],Xtest)
  coverage_length_rf_calLCP=coverage_length(result_rf_calLCP[,1],Xtest)
  coverage_length_rf_baseLCP=coverage_length(result_rf_baseLCP[,1],Xtest)

  coverage_sex_rf_RLCP=coverage_sex(result_rf_RLCP,Xtest)
  coverage_sex_rf_calLCP=coverage_sex(result_rf_calLCP,Xtest)
  coverage_sex_rf_baseLCP=coverage_sex(result_rf_baseLCP,Xtest)
  
  #----------neural net--------------------
  model_nn=neuralnet(rings~.,data = train_scaled_data, hidden = c(5, 3),
                     threshold=0.05,linear.output = TRUE)
  predict_nn_calib=neuralnet::compute(model_nn,calib_scaled_data[,1:7])$net.result*
    (max_rings-min_rings)+min_rings
  predict_nn_test=neuralnet::compute(model_nn,test_scaled_data[,1:7])$net.result*
    (max_rings-min_rings)+min_rings
  scores_nn_calib=abs(calib_data$rings-predict_nn_calib)
  scores_nn_test=abs(test_data$rings-predict_nn_test)
  
  result_nn_calLCP=calLCP_real(Xcalib, scores_nn_calib,Xtest,scores_nn_test,h,alpha)
  result_nn_baseLCP=baseLCP_real(Xcalib, scores_nn_calib,Xtest,scores_nn_test,h,alpha)
  result_nn_RLCP=RLCP_real(Xcalib, scores_nn_calib,Xtest,scores_nn_test,h,alpha)
  
  coverage_nn_RLCP=mean(result_nn_RLCP[,1])
  coverage_nn_calLCP=mean(result_nn_calLCP[,1])
  coverage_nn_baseLCP=mean(result_nn_baseLCP[,1])
  
  width_nn_RLCP=2*median(result_nn_RLCP[,2])
  width_nn_calLCP=2*median(result_nn_calLCP[,2])
  width_nn_baseLCP=2*median(result_nn_baseLCP[,2])
  
  coverage_diam_nn_RLCP=coverage_diam(result_nn_RLCP[,1],Xtest)
  coverage_diam_nn_calLCP=coverage_diam(result_nn_calLCP[,1],Xtest)
  coverage_diam_nn_baseLCP=coverage_diam(result_nn_baseLCP[,1],Xtest)
  
  coverage_length_nn_RLCP=coverage_length(result_nn_RLCP[,1],Xtest)
  coverage_length_nn_calLCP=coverage_length(result_nn_calLCP[,1],Xtest)
  coverage_length_nn_baseLCP=coverage_length(result_nn_baseLCP[,1],Xtest)

  coverage_sex_nn_RLCP=coverage_sex(result_nn_RLCP,Xtest)
  coverage_sex_nn_calLCP=coverage_sex(result_nn_calLCP,Xtest)
  coverage_sex_nn_baseLCP=coverage_sex(result_nn_baseLCP,Xtest)
  
  return(list(coverage=c(coverage_lm_RLCP,coverage_lm_calLCP,coverage_lm_baseLCP,
                         coverage_rf_RLCP,coverage_rf_calLCP,coverage_rf_baseLCP,
                         coverage_nn_RLCP,coverage_nn_calLCP,coverage_nn_baseLCP),
              width=c(width_lm_RLCP,width_lm_calLCP,width_lm_baseLCP,
                      width_rf_RLCP,width_rf_calLCP,width_rf_baseLCP,
                      width_nn_RLCP,width_nn_calLCP,width_nn_baseLCP),
              sex_coverage=cbind(coverage_sex_lm_RLCP,coverage_sex_lm_calLCP,coverage_sex_lm_baseLCP,
                                    coverage_sex_rf_RLCP,coverage_sex_rf_calLCP,coverage_sex_rf_baseLCP,
                                    coverage_sex_nn_RLCP,coverage_sex_nn_calLCP,coverage_sex_nn_baseLCP),
              diam_coverage=cbind(coverage_diam_lm_RLCP,coverage_diam_lm_calLCP,coverage_diam_lm_baseLCP,
                                 coverage_diam_rf_RLCP,coverage_diam_rf_calLCP,coverage_diam_rf_baseLCP,
                                 coverage_diam_nn_RLCP,coverage_diam_nn_calLCP,coverage_diam_nn_baseLCP),
              length_coverage=cbind(coverage_length_lm_RLCP,coverage_length_lm_calLCP,coverage_length_lm_baseLCP,
                                 coverage_length_rf_RLCP,coverage_length_rf_calLCP,coverage_length_rf_baseLCP,
                                 coverage_length_nn_RLCP,coverage_length_nn_calLCP,coverage_length_nn_baseLCP)
              
  ))
}

library(doParallel)
numcores=detectCores()-1
cl=makeCluster(numcores)
registerDoParallel(cl)
hseq=1:6 *0.05

real_result=vector('list',length=length(hseq))
mc=width=matrix(0,nrow=length(hseq),ncol=9)


comb=function(x,y){return(list(rbind(x[[1]],y[[1]]),x[[2]]+y[[2]],x[[3]]+y[[3]],
                               x[[4]]+y[[4]],x[[5]]+y[[5]]))}
nrep=50;alpha=0.1
for(i in 1:length(hseq)){
  h=hseq[i]
  result_h=foreach(k=1:nrep,.combine=comb,
                   .packages = c("MASS","parallel","doParallel",
                                 "foreach","mvtnorm","randomForest",
                                 "neuralnet")) %dopar% {
                                   real_analysis(h,k) 
                                 }
  result_h=lapply(result_h,FUN=function(x) x/nrep)
  real_result[[i]]=result_h
  width[i,]=real_result[[i]][[2]]
}
stopCluster(cl)

#----------------------------------------------------------
#---------------Visualization------------------------------
#----------------------------------------------------------


#----------------marginal coverage------------------
plot_result=matrix(0,ncol=5,nrow=9*length(hseq))
for(i in 1:length(hseq)){
  id=(9*(i-1)+1):(9*i)
  res=stack(as.data.frame(apply(real_result[[i]][[1]]*nrep,2,mean)))
  plot_result[id,1]=res[,1]
  plot_result[id,2]=rep(c("RLCP","calLCP","baseLCP"),3)
  plot_result[id,3]=rep(c("Linear Model","Random Forest","Neural Net"),each=3)
  plot_result[id,4]=hseq[i]
  plot_result[id,5]=stack(as.data.frame(apply(real_result[[i]][[1]]*nrep,2,sd)))[,1]/sqrt(nrep)
}
plot_result=as.data.frame(plot_result)
colnames(plot_result)=c("coverage","CP_method","base_method","h","se")
plot_result$coverage=as.numeric(plot_result$coverage)
plot_result$h=as.numeric(plot_result$h)
plot_result$se=as.numeric(plot_result$se)
plot_result$CP_method=factor(plot_result$CP_method,level=c('baseLCP','calLCP','RLCP'))

plot1=ggplot(plot_result, aes(x=h, y = coverage, color=CP_method)) +
  geom_line() + 
  geom_point()+
  ylim(0.885,0.915)+
  geom_hline(yintercept=0.9,linetype="dashed")+
  scale_color_manual(values=c("RLCP"="maroon","calLCP"="gold3","baseLCP"="blue"))+
  facet_grid(.~base_method)+
  scale_x_continuous(n.breaks=3)+
  theme(legend.position = "none")+
  labs(color="Method",x="Bandwidth h",y="Coverage")+
  theme(axis.title = element_text(size = 14),
        axis.text=element_text(size=13),
        plot.title = element_text(hjust = 0.5,size=15), 
        legend.position = "bottom",
        legend.text = element_text(size=15),
        strip.text = element_text(size=13))


#--------------------PI widths------------------------------
write.csv(round(t(width),2),"../results/real_data_width.csv")

plot_result=matrix(0,ncol=4,nrow=9*length(hseq))
for(i in 1:length(hseq)){
  id=(9*(i-1)+1):(9*i)
  res=stack(as.data.frame(real_result[[i]][[2]]))
  plot_result[id,1]=res[,1]
  plot_result[id,2]=rep(c("RLCP","calLCP","baseLCP"),3)
  plot_result[id,3]=rep(c("Linear Model","Random Forest","Neural Net"),each=3)
  plot_result[id,4]=hseq[i]
}
plot_result=as.data.frame(plot_result)
colnames(plot_result)=c("width","CP_method","base_method","h")
plot_result$width=as.numeric(plot_result$width)
plot_result$h=as.numeric(plot_result$h)
plot_result$CP_method=factor(plot_result$CP_method,level=c('baseLCP','calLCP','RLCP'))

plot2=ggplot(plot_result, aes(x=h, y = width, color=CP_method)) +
  geom_line() + 
  geom_point()+
  ylim(7.2,8.5)+
  scale_x_continuous(n.breaks=3)+
  scale_color_manual(values=c("RLCP"="maroon","calLCP"="gold3","baseLCP"="blue"))+
  facet_grid(.~base_method)+
  labs(color="Method",x="Bandwidth h",y="Prediction interval width")+
  theme(axis.title = element_text(size = 14),
        axis.text=element_text(size=13),
        plot.title = element_text(hjust = 0.5,size=15), 
        legend.position = "bottom",
        legend.text = element_text(size=15),
        strip.text = element_text(size=13))



pdf(file = "../results/figures/realdata_marginal_results.pdf",width = 10, height = 4.5) 
grid.arrange(plot1,plot2,ncol=2,widths=c(3.3,3.3))
dev.off()



#-------------------------------------------------------------
#-----------------------Local coverage results----------------
#-------------------------------------------------------------

sex_cov=matrix(0,nrow=3*length(hseq),ncol=9)
for(i in 1:length(hseq)){
  sex_cov[(3*(i-1)+1):(3*i),]=matrix(unlist(real_result[[i]][3]),nrow=3,byrow=F)
}
write.csv(sex_cov,"../results/real_data_sex_cov.csv")

#---------------across smoking groups------------------
plot_result=matrix(0,ncol=5,nrow=3*9*length(hseq))
for(i in 1:length(hseq)){
  id=(3*9*(i-1)+1):(3*9*i)
  res=stack(as.data.frame(real_result[[i]][[3]]))
  plot_result[id,1]=res[,1]
  plot_result[id,2]=rep(rep(c("RLCP","calLCP","baseLCP"),each=3),3)
  plot_result[id,3]=rep(c("Linear Model","Random Forest","Neural Net"),each=3*3)
  plot_result[id,4]=rep(c("F","M","I"),9)
  plot_result[id,5]=hseq[i]
}
plot_result=as.data.frame(plot_result)
colnames(plot_result)=c("coverage","CP_method","base method","sex","h")
plot_result$coverage=as.numeric(plot_result$coverage)
plot_result$h=as.numeric(plot_result$h)
plot_result$CP_method=factor(plot_result$CP_method,level=c('baseLCP','calLCP','RLCP'))

plot3=ggplot(plot_result[plot_result$h %in% c(0.05,0.25),], aes(x = sex, y = coverage,group=CP_method,
                                                          fill=CP_method)) +
  geom_col(position='dodge',width=0.5) + 
  coord_cartesian(ylim=c(0.88,0.93))+
  xlab("LCP method")+
  scale_fill_manual(values=c("RLCP"="maroon","calLCP"="gold3","baseLCP"="blue"))+
  geom_hline(yintercept=0.9,linetype="dashed")+
  facet_grid(h~`base method`,labeller = label_bquote(cols = .(`base method`),rows= h ==.(h)))+
  theme(axis.title = element_text(size = 16),
        axis.text=element_text(size=15),
        plot.title = element_text(hjust = 0.5,size=17), 
        legend.position = "bottom",
        legend.text = element_text(size=17),
        strip.text = element_text(size=15))+
  labs(fill="Method",x="Sex",y="Coverage")

#-------------across diameter bins-------------------------
quantiles=seq(0.125,0.875,by=0.01);len_q=length(quantiles)
plot_result=matrix(0,ncol=5,nrow=len_q*9*length(hseq))
for(i in 1:length(hseq)){
  id=(len_q*9*(i-1)+1):(len_q*9*i)
  res=stack(as.data.frame(real_result[[i]][[4]]))
  plot_result[id,1]=res[,1]
  plot_result[id,2]=rep(rep(c("RLCP","calLCP","baseLCP"),each=len_q),3)
  plot_result[id,3]=rep(c("Linear Model","Random Forest","Neural Net"),each=len_q*3)
  plot_result[id,4]=rep(quantiles,9)
  plot_result[id,5]=hseq[i]
}
plot_result=as.data.frame(plot_result)
colnames(plot_result)=c("coverage","CP_method","base method","diam","h")
plot_result$coverage=as.numeric(plot_result$coverage)
plot_result$diam=as.numeric(plot_result$diam)
plot_result$h=as.numeric(plot_result$h)
plot_result$CP_method=factor(plot_result$CP_method,level=c('baseLCP','calLCP','RLCP'))

plot4=ggplot(plot_result[plot_result$h %in% c(0.05,0.25),], aes(x = diam, y = coverage,color=CP_method)) +
  geom_line() + 
  geom_hline(yintercept=0.9,linetype="dashed")+
  xlab("Diameter quantile")+
  scale_color_manual(values=c("RLCP"="maroon","calLCP"="gold3","baseLCP"="blue"))+
  facet_grid(h~ `base method` ,labeller=label_bquote(cols = .(`base method`),rows= h ==.(h)))+
  theme(axis.title = element_text(size = 16),
        axis.text=element_text(size=15),
        plot.title = element_text(hjust = 0.5,size=17), 
        legend.position = "bottom",
        legend.text = element_text(size=17),
        strip.text = element_text(size=15))+
  labs(color="Method",x="Diameter Quantile",y="Local Coverage")

#-------------across length bins-------------------------
quantiles=seq(0.125,0.875,by=0.01);len_q=length(quantiles)
plot_result=matrix(0,ncol=5,nrow=len_q*9*length(hseq))
for(i in 1:length(hseq)){
  id=(len_q*9*(i-1)+1):(len_q*9*i)
  res=stack(as.data.frame(real_result[[i]][[5]]))
  plot_result[id,1]=res[,1]
  plot_result[id,2]=rep(rep(c("RLCP","calLCP","baseLCP"),each=len_q),3)
  plot_result[id,3]=rep(c("Linear Model","Random Forest","Neural Net"),each=len_q*3)
  plot_result[id,4]=rep(quantiles,9)
  plot_result[id,5]=hseq[i]
}
plot_result=as.data.frame(plot_result)
colnames(plot_result)=c("coverage","CP_method","base method","length","h")
plot_result$coverage=as.numeric(plot_result$coverage)
plot_result$length=as.numeric(plot_result$length)
plot_result$h=as.numeric(plot_result$h)
plot_result$CP_method=factor(plot_result$CP_method,level=c('baseLCP','calLCP','RLCP'))

plot5=ggplot(plot_result[plot_result$h %in% c(0.05,0.25),], aes(x = length, y = coverage,color=CP_method)) +
  geom_line() + 
  geom_hline(yintercept=0.9,linetype="dashed")+
  xlab("Length quantile")+
  scale_color_manual(values=c("RLCP"="maroon","calLCP"="gold3","baseLCP"="blue"))+
  facet_grid(h~ `base method` ,labeller=label_bquote(cols = .(`base method`),rows= h ==.(h)))+
  theme(axis.title = element_text(size = 16),
        axis.text=element_text(size=15),
        plot.title = element_text(hjust = 0.5,size=17), 
        legend.position = "bottom",
        legend.text = element_text(size=17),
        strip.text = element_text(size=15))+
  labs(color="Method",x="Length Quantile",y="Local Coverage")

pdf(file = "../results/figures/realdata_conditional_results.pdf",width = 14, height = 5) 
grid.arrange(plot3,plot5,ncol=2,widths=c(3,3))
dev.off()


