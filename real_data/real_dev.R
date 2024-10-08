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

#--------------------------------------------------------------
#---------------computing deviation of RLCP PI-----------------
#--------------------------------------------------------------
real_RLCP_deviation=function(h,k,split){
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
  
  result_lm_RLCP=RLCP_real(Xcalib, scores_lm_calib,Xtest,scores_lm_test,h,alpha)
  width_lm_RLCP=2*abs(result_lm_RLCP[,2])
  
  #----random forest----------
  model_rf=randomForest(rings ~ .,data=train_data)
  predict_rf_test=predict(model_rf,newdata=test_data)
  scores_rf_calib=abs(calib_data$rings-predict(model_rf,calib_data))
  scores_rf_test=abs(test_data$rings-predict(model_rf,test_data))
  
  result_rf_RLCP=RLCP_real(Xcalib, scores_rf_calib,Xtest,scores_rf_test,h,alpha)
  width_rf_RLCP=2*abs(result_rf_RLCP[,2])
  
  #----------neural net--------------------
  model_nn=neuralnet(rings~.,data = train_scaled_data, hidden = c(5, 3),
                      threshold=0.05,linear.output = TRUE)
  predict_nn_calib=neuralnet::compute(model_nn,calib_scaled_data[,1:7])$net.result*
    (max_rings-min_rings)+min_rings
  predict_nn_test=neuralnet::compute(model_nn,test_scaled_data[,1:7])$net.result*
    (max_rings-min_rings)+min_rings
  scores_nn_calib=abs(calib_data$rings-predict_nn_calib)
  scores_nn_test=abs(test_data$rings-predict_nn_test)
  
  result_nn_RLCP=RLCP_real(Xcalib, scores_nn_calib,Xtest,scores_nn_test,h,alpha)
  width_nn_RLCP=as.vector(2*abs(result_nn_RLCP[,2]))
  
  return(list(width_lm_RLCP,width_rf_RLCP,width_nn_RLCP))
}
#--------------------------------------------------------------
#----------------------experimental results--------------------
#--------------------------------------------------------------
library(doParallel)
numcores=detectCores()
cl=makeCluster(numcores)
registerDoParallel(cl)
#bandwidth choices
hseq=1:6 *0.05
alpha=0.1

comb_rand=function(x,y){return(list(rbind(x[[1]],y[[1]]),rbind(x[[2]],y[[2]]),rbind(x[[3]],y[[3]])))}

result_rand=matrix(0,ncol=length(hseq),nrow=3)
for(i in 1:length(hseq)){
  h=hseq[i]
  print(h)
  result_h=rep(0,3)
  start=Sys.time()
  for(k in 1:20){
    split=sample(1:3,n,replace=TRUE,prob=c(1/3,1/3,1/3))
    result_split=foreach(l=1:30,.combine=comb_rand,
                         .packages = c("MASS","parallel","doParallel",
                                       "foreach","mvtnorm","randomForest",
                                       "neuralnet")) %dopar% {
                                         real_RLCP_deviation(h,l,split)
                                       }
    #computing deviation of the PI across several random seeds
    deviation=function(x){median(abs(x-median(x)))/median(x)}
    
    width_lm_RLCP=result_split[[1]]
    width_rf_RLCP=result_split[[2]]
    width_nn_RLCP=result_split[[3]]
    dev_lm=mean(apply(width_lm_RLCP[,apply(width_lm_RLCP,2,median)!=Inf],2,deviation))
    dev_rf=mean(apply(width_rf_RLCP[,apply(width_rf_RLCP,2,median)!=Inf],2,deviation))
    dev_nn=mean(apply(width_nn_RLCP[,apply(width_nn_RLCP,2,median)!=Inf],2,deviation))
    
    result_h=result_h+c(dev_lm,dev_rf,dev_nn)
  }
  result_rand[,i]=result_h/20
  end=Sys.time()
  print(end-start)
}
stopCluster(cl)
#--------------------------------------------------------------
#------------------------visualization-------------------------
#--------------------------------------------------------------
deviation_real_df=c(result_rand[1,],result_rand[2,],result_rand[3,])
deviation_real_df=cbind(deviation_real_df,rep(hseq,3))
deviation_real_df=cbind(deviation_real_df,rep(c("Linear Model","Random Forest","Neural Net"),each=length(hseq)))
deviation_real_df=as.data.frame(deviation_real_df)
colnames(deviation_real_df)=c("deviation","h","setting")
deviation_real_df$setting=as.factor(deviation_real_df$setting)
deviation_real_df$deviation=as.numeric(deviation_real_df$deviation)
deviation_real_df$h=as.numeric(deviation_real_df$h)

deviation_real=ggplot(deviation_real_df,aes(x=h,y=deviation,linetype=setting))+
  geom_line()+geom_point()+
  scale_linetype_manual(values=c("Linear Model"="dashed","Random Forest"="dotted","Neural Net"="dotdash"))+
  labs(x="Bandwidth h",y=expression(D(h)),linetype="Base Method")+
  ggtitle("Deviation D(h) of RLCP in real data setting")+
  theme(axis.title = element_text(size = 16),
        axis.text=element_text(size=16),
        plot.title = element_text(hjust = 0.5,size=16), 
        legend.position = "bottom",
        legend.text = element_text(size=14))

pdf(file = "../results/figures/RLCP_deviation.pdf",width = 13.5, height = 5)
gridExtra::grid.arrange(deviation_simul,deviation_real,ncol=2)
dev.off()

