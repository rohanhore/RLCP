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

#--------------------------------------------------------
#--------data loading & pre-processing-------------------
#--------------------------------------------------------
data=read.csv("Train_Data.csv",header=T)
data=as.data.frame(data)
data=data %>% distinct()
#disregarding the children and region information
data=data[,-c(5,6)]
#log transformation on insurance charges
data[,5]=log(data[,5])

n=dim(data)[1]
data$sex=as.factor(data$sex);data$smoker=as.factor(data$smoker)
data$smoker=dplyr::recode_factor(data$smoker,"yes" = "1","no" = "0")
data$sex=dplyr::recode_factor(data$sex,"male" = "1","female" = "0")
attach(data);sex=as.factor(sex);smoker=as.factor(smoker)

#scaling data for neural nets
scaled_data=data
min_age=min(scaled_data[,1]);max_age=max(scaled_data[,1])
min_bmi=min(scaled_data[,3]);max_bmi=max(scaled_data[,3])
scaled_data[,1]=(scaled_data[,1]-min_age)/(max_age-min_age)
scaled_data[,3]=(scaled_data[,3]-min_bmi)/(max_bmi-min_bmi)
min_charge=min(scaled_data[,5]);max_charge=max(scaled_data[,5])
scaled_data[,5]=(scaled_data[,5]-min_charge)/(max_charge-min_charge)
dmy=dummyVars(" ~ .", data = scaled_data)
scaled_data <- data.frame(predict(dmy, newdata = scaled_data))


###-----computing local coverage across BMI bins-----------
coverage_bmi=function(coverage,test_data){
  quantile_points=seq(0.125,0.875,by=0.01)
  coverage_bmi=rep(0,length(quantile_points))
  for(i in 1:length(quantile_points)){
    id=which((test_data$bmi>=quantile(data$bmi,probs=quantile_points[i]-0.025)) & 
               (test_data$bmi<=quantile(data$bmi,probs=quantile_points[i]+0.025)))
    coverage_bmi[i]=mean(coverage[id])
  }
  return(coverage_bmi)
}

###-----local coverage across smoking groups-----------
coverage_smoker=function(coverage,test_data){
  test_data$smoker=as.factor(test_data$smoker)
  coverage_smoker=rep(0,2)
  for(i in 1:2){
    id=which(test_data$smoker==c(0,1)[i])
    coverage_smoker[i]=mean(coverage[id])
  }
  return(coverage_smoker)
}

#----------------------------------------------------------
#-----------------RLCP for real data-----------------------
#----------------------------------------------------------

RLCP_real=function(Xcalib,Vcalib,Xtest,Vtest,h,alpha){
  ntest=dim(Xtest)[1];d=dim(Xtest)[2]
  coverage=cutoff=rep(0,ntest)
  Xcalib=as.matrix(Xcalib[order(Vcalib),]);Vcalib=sort(Vcalib)
  
  scores=c(Vcalib,Inf)
  indices=list();j=1;i=1
  scores_unique=vector()
  while(i<=length(scores)){
    scores_unique=c(scores_unique,scores[i])
    indices[[j]]=which(scores==scores[i])
    i=i+length(indices[[j]]);j=j+1
  }
  
  for(i in 1:ntest){
    xtest=Xtest[i,];vtest=Vtest[i]
    xtilde_test=xtest
    xtilde_test[c(1,3)]=xtest[c(1,3)]+runif(2,min=-h,max=h)
    
    cov_data=rbind(Xcalib,xtest)
    
    weights=apply(abs(sweep(cov_data,2,as.numeric(xtilde_test),"-")),1,FUN=function(x) all(x<=c(h,0,h,0))+0)
    result=smoothed_weighted_quantile(scores_unique,alpha,weights,indices)
    
    cutoff[i]=result[1]
    closed=result[2]
    coverage[i]=(vtest<cutoff[i])+0
    if(closed==TRUE){coverage[i]=(vtest<=cutoff[i])+0}
  }
  return(cbind(coverage,cutoff))
}

#-----------------------------------------------------------
#-----------------baseLCP for real data---------------------
#-----------------------------------------------------------

baseLCP_real=function(Xcalib,Vcalib,Xtest,Vtest,h,alpha){    
  ntest=dim(Xtest)[1];d=dim(Xtest)[2]
  coverage=cutoff=rep(0,ntest)
  Xcalib=as.matrix(Xcalib[order(Vcalib),]);Vcalib=sort(Vcalib)
  
  scores=c(Vcalib,Inf)
  indices=list();j=1;i=1
  scores_unique=vector()
  while(i<=length(scores)){
    scores_unique=c(scores_unique,scores[i])
    indices[[j]]=which(scores==scores[i])
    i=i+sum(scores==scores[i]);j=j+1
  }
  
  for(i in 1:ntest){
    xtest=Xtest[i,];vtest=Vtest[i]
    cov_data=rbind(Xcalib,xtest)
    
    weights=apply(sweep(cov_data,2,as.numeric(xtest)),1,FUN=function(x) all(x<=c(h,0,h,0))+0)
    result=smoothed_weighted_quantile(scores_unique,alpha,weights,indices)
    cutoff[i]=result[1]
    closed=result[2]
    coverage[i]=(vtest<cutoff[i])+0
    if(closed==TRUE){coverage[i]=(vtest<=cutoff[i])+0}
  }
  return(cbind(coverage,cutoff))
}

#------------------------------------------------------------
#-----------------calLCP for real data-----------------------
#------------------------------------------------------------

calLCP_real=function(Xcalib,Vcalib,Xtest,Vtest,h,alpha){
  ntest=dim(Xtest)[1];d=dim(Xtest)[2];ncalib=dim(Xcalib)[1]
  coverage=cutoff=rep(0,ntest)
  Xcalib=as.matrix(Xcalib[order(Vcalib),])
  Vcalib=(sort(Vcalib))
  
  H=matrix(0,nrow=ncalib,ncol=ncalib)
  for(i in 1:ncalib){
    H[i,]=apply(sweep(Xcalib,2,as.numeric(Xcalib[i,])),1,FUN=function(x) all(x<=c(h,0,h,0))+0)
  }
  denom_calib=apply(H,1,sum)
  num_calib=rep(0,ncalib)
  for(i in 1:ncalib){num_calib[i]=sum(H[i,]*(Vcalib<Vcalib[i]))}

  p_values=rep(0,ncalib+1)
  for(i in 1:ntest){
    U=runif(1,0,1)
    
    smoothed_p_value=function(x){
      return(sum(x>tail(x,1))/length(x)+(U*sum(x==tail(x,1)))/length(x))
    }
    
    xtest=Xtest[i,];vtest=Vtest[i]
    cov_data=rbind(Xcalib,xtest)
    scores=c(Vcalib,Inf)
    
    weights=apply(sweep(cov_data,2,as.numeric(xtest)),1,FUN=function(x) all(x<=c(h,0,h,0))+0)
    
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
    id=max(which(p_values>alpha));closed=TRUE
    if(id<=ncalib){
      T_values[ncalib+1]=sum(weights*(scores<scores[id]))/sum(weights)
      for(j in 1:id){T_values[j]=(num_calib[j])/(denom_calib[j]+weights[j])}
      if(id<ncalib){
        for(j in (id+1):ncalib){T_values[j]=(num_calib[j])/(denom_calib[j]+weights[j])}}
      closed=(smoothed_p_value(T_values)>alpha)
    }
      cutoff[i]=scores[id]
    
    if(closed==TRUE){coverage[i]=(vtest<=cutoff[i])+0}
    else{coverage[i]=(vtest<cutoff[i])+0}
  }
  return(cbind(coverage,cutoff))
}

#-------------------------------------------------------------------------------
#--------------analyzing marginal,local coverages for real data-----------------
#-------------------------------------------------------------------------------
real_analysis=function(h,k){
  split=sample(1:3,n,replace=TRUE,prob=c(1/5,2/5,2/5))
  train_data=data[split==1,];ntrain=dim(train_data)[1]
  calib_data=data[split==2,];ncalib=dim(calib_data)[1]
  test_data=data[split==3,];ntest=dim(test_data)[1]
  
  train_scaled_data=scaled_data[split==1,]
  calib_scaled_data=scaled_data[split==2,]
  test_scaled_data=scaled_data[split==3,]
  
  #------------learning score on train split-------------------------
  Xcalib=calib_data[,1:4]
  Xcalib[,2]=as.numeric(levels(calib_data[,2]))[calib_data[,2]]
  Xcalib[,4]=as.numeric(levels(calib_data[,4]))[calib_data[,4]]
  
  Xtest=test_data[,1:4]
  Xtest[,2]=as.numeric(levels(test_data[,2]))[test_data[,2]]
  Xtest[,4]=as.numeric(levels(test_data[,4]))[test_data[,4]]
  
  #---linear model---------
  model_lm=lm(charges~.,data=train_data)
  predict_lm_test=predict(model_lm,newdata=test_data)
  scores_lm_calib=abs(calib_data$charges-predict.lm(model_lm,calib_data))
  scores_lm_test=abs(test_data$charges-predict.lm(model_lm,test_data))
  
  result_lm_calLCP=calLCP_real(Xcalib, scores_lm_calib,Xtest,scores_lm_test,h,0.1)
  result_lm_baseLCP=baseLCP_real(Xcalib, scores_lm_calib,Xtest,scores_lm_test,h,0.1)
  result_lm_RLCP=RLCP_real(Xcalib, scores_lm_calib,Xtest,scores_lm_test,h,0.1)
  
  coverage_lm_RLCP=mean(result_lm_RLCP[,1])
  coverage_lm_calLCP=mean(result_lm_calLCP[,1])
  coverage_lm_baseLCP=mean(result_lm_baseLCP[,1])
  
  width_lm_RLCP=2*median(result_lm_RLCP[,2])
  width_lm_calLCP=2*median(result_lm_calLCP[,2])
  width_lm_baseLCP=2*median(result_lm_baseLCP[,2])
  
  coverage_bmi_lm_RLCP=coverage_bmi(result_lm_RLCP[,1],Xtest)
  coverage_bmi_lm_calLCP=coverage_bmi(result_lm_calLCP[,1],Xtest)
  coverage_bmi_lm_baseLCP=coverage_bmi(result_lm_baseLCP[,1],Xtest)
  
  coverage_smoker_lm_RLCP=coverage_smoker(result_lm_RLCP,Xtest)
  coverage_smoker_lm_calLCP=coverage_smoker(result_lm_calLCP,Xtest)
  coverage_smoker_lm_baseLCP=coverage_smoker(result_lm_baseLCP,Xtest)
  
  
  #----random forest----------
  model_rf=randomForest(charges ~ .,data=train_data)
  predict_rf_test=predict(model_rf,newdata=test_data)
  scores_rf_calib=abs(calib_data$charges-predict(model_rf,calib_data))
  scores_rf_test=abs(test_data$charges-predict(model_rf,test_data))
  
  result_rf_calLCP=calLCP_real(Xcalib, scores_rf_calib,Xtest,scores_rf_test,h,0.1)
  result_rf_baseLCP=baseLCP_real(Xcalib, scores_rf_calib,Xtest,scores_rf_test,h,0.1)
  result_rf_RLCP=RLCP_real(Xcalib, scores_rf_calib,Xtest,scores_rf_test,h,0.1)
  
  coverage_rf_RLCP=mean(result_rf_RLCP[,1])
  coverage_rf_calLCP=mean(result_rf_calLCP[,1])
  coverage_rf_baseLCP=mean(result_rf_baseLCP[,1])
  
  width_rf_RLCP=2*median(result_rf_RLCP[,2])
  width_rf_calLCP=2*median(result_rf_calLCP[,2])
  width_rf_baseLCP=2*median(result_rf_baseLCP[,2])
  
  coverage_bmi_rf_RLCP=coverage_bmi(result_rf_RLCP[,1],Xtest)
  coverage_bmi_rf_calLCP=coverage_bmi(result_rf_calLCP[,1],Xtest)
  coverage_bmi_rf_baseLCP=coverage_bmi(result_rf_baseLCP[,1],Xtest)
  
  coverage_smoker_rf_RLCP=coverage_smoker(result_rf_RLCP,Xtest)
  coverage_smoker_rf_calLCP=coverage_smoker(result_rf_calLCP,Xtest)
  coverage_smoker_rf_baseLCP=coverage_smoker(result_rf_baseLCP,Xtest)
  
  #----------neural net--------------------
  model_nn=neuralnet(charges~.,data = train_scaled_data, hidden = c(5, 3),
                     threshold=0.05,linear.output = TRUE)
  predict_nn_calib=neuralnet::compute(model_nn,calib_scaled_data[,1:6])$net.result*
    (max_charge-min_charge)+min_charge
  predict_nn_test=neuralnet::compute(model_nn,test_scaled_data[,1:6])$net.result*
    (max_charge-min_charge)+min_charge
  scores_nn_calib=abs(calib_data$charges-predict_nn_calib)
  scores_nn_test=abs(test_data$charges-predict_nn_test)
  
  result_nn_calLCP=calLCP_real(Xcalib, scores_nn_calib,Xtest,scores_nn_test,h,0.1)
  result_nn_baseLCP=baseLCP_real(Xcalib, scores_nn_calib,Xtest,scores_nn_test,h,0.1)
  result_nn_RLCP=RLCP_real(Xcalib, scores_nn_calib,Xtest,scores_nn_test,h,0.1)
  
  coverage_nn_RLCP=mean(result_nn_RLCP[,1])
  coverage_nn_calLCP=mean(result_nn_calLCP[,1])
  coverage_nn_baseLCP=mean(result_nn_baseLCP[,1])
  
  width_nn_RLCP=2*median(result_nn_RLCP[,2])
  width_nn_calLCP=2*median(result_nn_calLCP[,2])
  width_nn_baseLCP=2*median(result_nn_baseLCP[,2])
  
  coverage_bmi_nn_RLCP=coverage_bmi(result_nn_RLCP[,1],Xtest)
  coverage_bmi_nn_calLCP=coverage_bmi(result_nn_calLCP[,1],Xtest)
  coverage_bmi_nn_baseLCP=coverage_bmi(result_nn_baseLCP[,1],Xtest)
  
  coverage_smoker_nn_RLCP=coverage_smoker(result_nn_RLCP,Xtest)
  coverage_smoker_nn_calLCP=coverage_smoker(result_nn_calLCP,Xtest)
  coverage_smoker_nn_baseLCP=coverage_smoker(result_nn_baseLCP,Xtest)
  
  return(list(coverage=c(coverage_lm_RLCP,coverage_lm_calLCP,coverage_lm_baseLCP,
                         coverage_rf_RLCP,coverage_rf_calLCP,coverage_rf_baseLCP,
                         coverage_nn_RLCP,coverage_nn_calLCP,coverage_nn_baseLCP),
              width=c(width_lm_RLCP,width_lm_calLCP,width_lm_baseLCP,
                      width_rf_RLCP,width_rf_calLCP,width_rf_baseLCP,
                      width_nn_RLCP,width_nn_calLCP,width_nn_baseLCP),
              bmi_coverage=cbind(coverage_bmi_lm_RLCP,coverage_bmi_lm_calLCP,coverage_bmi_lm_baseLCP,
                                 coverage_bmi_rf_RLCP,coverage_bmi_rf_calLCP,coverage_bmi_rf_baseLCP,
                                 coverage_bmi_nn_RLCP,coverage_bmi_nn_calLCP,coverage_bmi_nn_baseLCP),
              smoker_coverage=cbind(coverage_smoker_lm_RLCP,coverage_smoker_lm_calLCP,coverage_smoker_lm_baseLCP,
                                    coverage_smoker_rf_RLCP,coverage_smoker_rf_calLCP,coverage_smoker_rf_baseLCP,
                                    coverage_smoker_nn_RLCP,coverage_smoker_nn_calLCP,coverage_smoker_nn_baseLCP)
  ))
}

library(doParallel)
numcores=detectCores()-1
cl=makeCluster(numcores)
registerDoParallel(cl)
#bandwidth choices
hseq=2:8

real_result=vector('list',length=length(hseq))
mc=width=matrix(0,nrow=length(hseq),ncol=9)

comb=function(x,y){return(list(rbind(x[[1]],y[[1]]),x[[2]]+y[[2]],x[[3]]+y[[3]],
                               x[[4]]+y[[4]]))}
nrep=100
for(i in 1:length(hseq)){
  h=hseq[i]
  print(h)
  result_h=foreach(k=1:nrep,.combine=comb,
                   .packages = c("MASS","parallel","doParallel",
                                 "foreach","mvtnorm","randomForest",
                                 "neuralnet")) %dopar% {
                                   real_analysis(h,k)
                                 }
  result_h=lapply(result_h,FUN=function(x) x/nrep)
  real_result[[i]]=result_h
  print(colMeans(real_result[[i]][[1]]*nrep))
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
  ylim(0.88,0.98)+
  geom_errorbar(aes(ymin=coverage-se, ymax=coverage+se), width=.02,
                position=position_dodge(0.01))+
  geom_hline(yintercept=0.9,linetype="dashed")+
  scale_color_manual(values=c("RLCP"="maroon","calLCP"="gold3","baseLCP"="blue"))+
  facet_grid(.~base_method)+
  scale_x_continuous(n.breaks=3)+
  theme(legend.position = "none")+
  labs(color="Method",x="Bandwidth h",y="Coverage")


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
  ylim(0.5,2.5)+
  scale_x_continuous(n.breaks=3)+
  scale_color_manual(values=c("RLCP"="maroon","calLCP"="gold3","baseLCP"="blue"))+
  facet_grid(.~base_method)+
  labs(color="Method",x="Bandwidth h",y="Prediction interval width")



pdf(file = "../results/realdata_marginal_results.pdf",width = 9, height = 4) 
grid.arrange(plot1,plot2,ncol=2,widths=c(2.6,3))
dev.off()

#-------------------------------------------------------------
#-----------------------Local coverage results----------------
#-------------------------------------------------------------

smoker_cov=matrix(0,nrow=2*length(hseq),ncol=9)
for(i in 1:length(hseq)){
  smoker_cov[(2*(i-1)+1):(2*i),]=matrix(unlist(real_result[[i]][5]),nrow=2,byrow=F)
}
write.csv(smoker_cov,"../results/real_data_smoker_cov.csv")

#---------------across smoking groups------------------
plot_result=matrix(0,ncol=5,nrow=2*9*length(hseq))
for(i in 1:length(hseq)){
  id=(2*9*(i-1)+1):(2*9*i)
  res=stack(as.data.frame(real_result[[i]][[4]]))
  plot_result[id,1]=res[,1]
  plot_result[id,2]=rep(rep(c("RLCP","calLCP","baseLCP"),each=2),3)
  plot_result[id,3]=rep(c("Linear Model","Random Forest","Neural Net"),each=2*3)
  plot_result[id,4]=rep(c(0,1),9)
  plot_result[id,5]=hseq[i]
}
plot_result=as.data.frame(plot_result)
colnames(plot_result)=c("coverage","CP_method","base method","smoker","h")
plot_result$coverage=as.numeric(plot_result$coverage)
plot_result$h=as.numeric(plot_result$h)
plot_result$CP_method=factor(plot_result$CP_method,level=c('baseLCP','calLCP','RLCP'))

plot3=ggplot(plot_result[plot_result$h %in% c(2,6),], aes(x = CP_method, y = coverage,group=smoker,
                                                          fill=smoker)) +
  geom_col(position='dodge',width=0.5) + 
  coord_cartesian(ylim=c(0.8,1))+
  xlab("LCP method")+
  scale_fill_discrete(labels=c("0"="no","1"="yes"))+
  geom_hline(yintercept=0.9,linetype="dashed")+
  facet_grid(h~`base method`,labeller = label_bquote(cols = .(`base method`),rows= h ==.(h)))+
  theme(legend.position = "bottom")+
  labs(fill="Smoker",x="Method",y="Coverage")

#-------------across BMI bins-------------------------
bmi_quantiles=seq(0.125,0.875,by=0.01);len_q=length(bmi_quantiles)
plot_result=matrix(0,ncol=5,nrow=len_q*9*length(hseq))
for(i in 1:length(hseq)){
  id=(len_q*9*(i-1)+1):(len_q*9*i)
  res=stack(as.data.frame(real_result[[i]][[3]]))
  plot_result[id,1]=res[,1]
  plot_result[id,2]=rep(rep(c("RLCP","calLCP","baseLCP"),each=len_q),3)
  plot_result[id,3]=rep(c("Linear Model","Random Forest","Neural Net"),each=len_q*3)
  plot_result[id,4]=rep(bmi_quantiles,9)
  plot_result[id,5]=hseq[i]
}
plot_result=as.data.frame(plot_result)
colnames(plot_result)=c("coverage","CP_method","base method","bmi","h")
plot_result$coverage=as.numeric(plot_result$coverage)
plot_result$bmi=as.numeric(plot_result$bmi)
plot_result$h=as.numeric(plot_result$h)
plot_result$CP_method=factor(plot_result$CP_method,level=c('baseLCP','calLCP','RLCP'))

plot4=ggplot(plot_result[plot_result$h %in% c(2,6),], aes(x = bmi, y = coverage,color=CP_method)) +
  geom_line() + 
  geom_hline(yintercept=0.9,linetype="dashed")+
  xlab("BMI quantile")+
  scale_color_manual(values=c("RLCP"="maroon","calLCP"="gold3","baseLCP"="blue"))+
  facet_grid(h~ `base method` ,labeller=label_bquote(cols = .(`base method`),rows= h ==.(h)))+
  theme(legend.position = "bottom")+
  labs(color="Method",x="BMI Quantile",y="Local Coverage")


pdf(file = "/Users/rohanhore/Dropbox/My projects/rLCP/Results/realdata_conditional_results.pdf",width = 13, height = 5) 
grid.arrange(plot3,plot4,ncol=2,widths=c(2.6,3))
dev.off()

