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
suppressPackageStartupMessages(library(ggplot2))
#------------------------------------------------------------------
#-----------visual analysis on effect of randomization-------------
#-------------------------------------------------------------------
randomization_effect=function(setting,nrep,h){
  ntrain=2000;ncalib=2000;d=1
  
  #----------training data and learning score------------
  train_data=simulation(ntrain,d,setting)
  Xnames=paste("X", 1:d, sep="")
  formula=as.formula(paste("Y ~ ", paste(Xnames, collapse= "+")))
  model_lm=lm(formula,data=train_data)
  
  #----------calibration data--------------------------
  calib_data=simulation(ncalib,d,setting)
  Xcalib=as.matrix(calib_data[,-1])
  scores_calib=abs(calib_data$Y-predict(model_lm,newdata=calib_data))
  
  #-------------test data------------------------------
  test_data=conditional_simulation(100,setting)
  Xtest=as.matrix(test_data[,-1])
  scores_test=abs(test_data$Y-predict(model_lm,newdata=test_data))
  
  #-----------evaluating the competing methods-----------
  calLCP_res=calLCP(Xcalib,scores_calib,Xtest,scores_test,"gaussian",h,alpha)
  
  scores_threshold_calLCP=calLCP_res[,2]
  scores_threshold_calLCP[scores_threshold_calLCP==Inf]=max(scores_calib,scores_test)+0.2
  pred_test=predict(model_lm,test_data)
  #calLCP prediction interval
  upper_threshold_calLCP=pred_test+scores_threshold_calLCP
  lower_threshold_calLCP=pred_test-scores_threshold_calLCP
  
  scores_threshold_RLCP=matrix(0,nrow=dim(Xtest)[1],ncol=nrep)
  for(k in 1:nrep){
    RLCP_res=RLCP(Xcalib,scores_calib,Xtest,scores_test,"gaussian",h,alpha)
    scores_threshold_RLCP[,k]=RLCP_res[,2]
    scores_threshold_RLCP[scores_threshold_RLCP[,k]==Inf,k]=max(scores_calib,scores_test)+0.2
  }
  #RLCP prediction intervals
  upper_thresholds_RLCP=pred_test+scores_threshold_RLCP
  lower_thresholds_RLCP=pred_test-scores_threshold_RLCP
  
  #95% and 5% quantiles of the RLCP prediction interval ends across different random seeds
  up_random1=smooth.spline(Xtest,apply(upper_thresholds_RLCP,1,FUN=function(x) {quantile(x,probs=0.95)}),spar=0.3)$y
  up_random2=smooth.spline(Xtest,apply(upper_thresholds_RLCP,1,FUN=function(x) {quantile(x,probs=0.05)}),spar=0.3)$y
  
  low_random1=smooth.spline(Xtest,apply(lower_thresholds_RLCP,1,FUN=function(x) {quantile(x,probs=0.95)}),spar=0.3)$y
  low_random2=smooth.spline(Xtest,apply(lower_thresholds_RLCP,1,FUN=function(x) {quantile(x,probs=0.05)}),spar=0.3)$y
  
  #oracle prediction interval
  if(setting==1){oracle=abs(sin(xseq))}
  if(setting==2){oracle=2*dnorm(xseq,0,1.5)}
  
  return(cbind(0.5*xseq+oracle*qnorm(0.95,0,1),0.5*xseq+oracle*qnorm(0.05,0,1),
               upper_threshold_calLCP,lower_threshold_calLCP,
               up_random1,up_random2,low_random1,low_random2))
}

#----------------------------------------------------------
#-----------------experimental results---------------------
#----------------------------------------------------------
#bandwidth choices
hseq=c(0.5,1,1.5)
#grid on feature space
xseq=seq(-3,3,by=0.01)

alpha=0.1
plot_result=data.frame()
for(i in 1:2){
  for(j in 1:3){
    result_band=randomization_effect(i,nrep=100,hseq[j])
    result_band=cbind(result_band,i)
    result_band=cbind(result_band,hseq[j])
    result_band=cbind(result_band,xseq)
    plot_result=rbind(plot_result,result_band)
  }
}
plot_result=cbind(plot_result,"threshold")
colnames(plot_result)=c("oracle_upper","oracle_lower","LCP_upper","LCP_lower","RLCP_UU","RLCP_UL","RLCP_LU","RLCP_LL","setting","h","covariate","ID")

data=simulation(4000,d=1,1)
result_scatter=matrix(0,ncol=12,nrow=4000)
for(i in 1:3){
  h=hseq[i]
  result_scatter[,1]=data$Y
  result_scatter[,9]=1
  result_scatter[,10]=h
  result_scatter[,11]=data$X1
  result_scatter[,12]="scatter"
  colnames(result_scatter)=c("oracle_upper","oracle_lower","LCP_upper","LCP_lower","RLCP_UU","RLCP_UL","RLCP_LU","RLCP_LL","setting","h","covariate","ID")
  plot_result=rbind(plot_result,result_scatter)
}

data=simulation(4000,d=1,2)
result_scatter=matrix(0,ncol=12,nrow=4000)
for(i in 1:3){
  h=hseq[i]
  result_scatter[,1]=data$Y
  result_scatter[,9]=2
  result_scatter[,10]=h
  result_scatter[,11]=data$X1
  result_scatter[,12]="scatter"
  colnames(result_scatter)=c("oracle_upper","oracle_lower","LCP_upper","LCP_lower","RLCP_UU","RLCP_UL","RLCP_LU","RLCP_LL","setting","h","covariate","ID")
  plot_result=rbind(plot_result,result_scatter)
}
write.csv(plot_result,"../results/simul_rlcp_band.csv")

#----------------------------------------------------------
#--------------------visualization-------------------------
#----------------------------------------------------------
plot_result=read.csv("../results/simul_rlcp_band.csv")[,-1]
colnames(plot_result)=c("oracle_upper","oracle_lower","LCP_upper","LCP_lower","RLCP_UU","RLCP_UL","RLCP_LU","RLCP_LL","setting","h","covariate","ID")
plot_result$oracle_upper=as.numeric(plot_result$oracle_upper)
plot_result$oracle_lower=as.numeric(plot_result$oracle_lower)
plot_result$LCP_upper=as.numeric(plot_result$LCP_upper)
plot_result$LCP_lower=as.numeric(plot_result$LCP_lower)
plot_result$RLCP_UU=as.numeric(plot_result$RLCP_UU)
plot_result$RLCP_UL=as.numeric(plot_result$RLCP_UL)
plot_result$RLCP_LU=as.numeric(plot_result$RLCP_LU)
plot_result$RLCP_LL=as.numeric(plot_result$RLCP_LL)
plot_result$covariate=as.numeric(plot_result$covariate)


pdf(file = "../results/figures/simul_randomization_effect.pdf",width = 12,height = 8)

cols<-c("RLCP"="maroon","oracle"="black","calLCP"="gold3")
linetypes<-c("RLCP"="solid",oracle="dotted","calLCP"="solid")
ggplot(plot_result[plot_result$ID %in% c("threshold"),]) +
  geom_point(data=plot_result[plot_result$ID %in% c("scatter"),],aes(x=covariate,y=oracle_upper),shape=1,col="gray",alpha=0.3)+
  facet_grid(setting~h ,labeller=label_bquote(rows = paste("Setting ", .(setting)),cols= h ==.(h)))+
  geom_line(aes(x=covariate,y=oracle_upper,color="oracle",linetype="oracle"),lwd=1)+
  geom_line(aes(x=covariate,y=oracle_lower,color="oracle",linetype="oracle"),lwd=1)+
  geom_line(aes(x=covariate,y=LCP_upper,color="calLCP",linetype="calLCP"),lwd=1.05)+
  geom_line(aes(x=covariate,y=LCP_lower,color="calLCP",linetype="calLCP"),lwd=1.05)+
  geom_line(aes(x=covariate,y=RLCP_UL,color="RLCP",linetype="RLCP"),lwd=0.5)+
  geom_line(aes(x=covariate,y=RLCP_UU,color="RLCP",linetype="RLCP"),lwd=0.5)+
  geom_line(aes(x=covariate,y=RLCP_LL,color="RLCP",linetype="RLCP"),lwd=0.5)+
  geom_line(aes(x=covariate,y=RLCP_LU,color="RLCP",linetype="RLCP"),lwd=0.5)+
  xlim(-3.5,3.5)+
  scale_color_manual(name="Method",values=cols,breaks=c("oracle","calLCP","RLCP"))+
  scale_linetype_manual(name="Method",values=linetypes,breaks=c("oracle","calLCP","RLCP"))+
  scale_y_continuous(n.breaks=6)+
  labs(color="Method",x="Feature X",y=" ")+
  geom_ribbon(aes(x=covariate,ymin=RLCP_LL,ymax=RLCP_LU),color="maroon",alpha=0.4)+
  geom_ribbon(aes(x=covariate,ymin=RLCP_UL,ymax=RLCP_UU),color="maroon",alpha=0.4)+
  theme_bw()+
  theme(strip.text = element_text(size = 18),
        legend.text=element_text(size=18),
        axis.text=element_text(size=18),
        axis.title=element_text(size=18),
        legend.title=element_text(size=18),
        legend.spacing.y = unit(0.3,'cm'))+
  guides(color = guide_legend(byrow = TRUE))
dev.off()

#---------------------------------------------------------------------------
#------------computing deviation of RLCP for simulation settings------------
#---------------------------------------------------------------------------
RLCP_deviation=function(setting,nrep,h,k){
  ntrain=ncalib=2000;d=1

  #----------training data and learning score------------
  train_data=simulation(ntrain,d,setting)
  Xnames=paste("X", 1:d, sep="")
  formula=as.formula(paste("Y ~ ", paste(Xnames, collapse= "+")))
  model_lm=lm(formula,data=train_data)
  
  #----------calibration data--------------------------
  calib_data=simulation(ncalib,d,setting)
  Xcalib=as.matrix(calib_data[,-1])
  scores_calib=abs(calib_data$Y-predict(model_lm,newdata=calib_data))
  
  #-------------test data------------------------------
  test_data=conditional_simulation(100,setting)
  Xtest=as.matrix(test_data[,-1])
  scores_test=abs(test_data$Y-predict(model_lm,newdata=test_data))
  
  #-----------evaluating the competing methods----------
  set.seed(NULL)
  scores_threshold_RLCP=matrix(0,nrow=dim(Xtest)[1],ncol=nrep)
  for(k in 1:nrep){
    RLCP_res=RLCP(Xcalib,scores_calib,Xtest,scores_test,"gaussian",h,alpha)
    scores_threshold_RLCP[,k]=RLCP_res[,2]
    scores_threshold_RLCP[scores_threshold_RLCP[,k]==Inf,k]=max(scores_calib,scores_test)+0.2
  }
  width=2*scores_threshold_RLCP
  deviation=function(x){median(abs(x-median(x)))/median(x)}
  
  return(mean(apply(width,1,deviation)))
}

library(doParallel)
numcores=detectCores()-1
cl=makeCluster(numcores)
registerDoParallel(cl)
hseq=seq(0.1,2.1,by=0.2)


deviation_est=matrix(0,ncol=length(hseq),nrow=2)
for(j in 1:2){
  for(i in 1:length(hseq)){
    h=hseq[i]
    print(h)
    result_h=foreach(k=1:30,.combine="+",
                     .packages = c("MASS","parallel","doParallel",
                                   "foreach","mvtnorm","randomForest",
                                   "neuralnet")) %dopar% {
                                     RLCP_deviation(setting=j,50,h,k)
                                   }
    deviation_est[j,i]=result_h/30
  }
}
stopCluster(cl)

write.csv(deviation_est,"../results/simul_deviation_estimates.csv")


#--------------------------------------------------------------
#---------------------------visualization----------------------
#--------------------------------------------------------------
deviation_est=read.csv("../results/simul_deviation_estimates.csv")[,-1]
deviation_simul_df=c(unlist(deviation_est[1,]),unlist(deviation_est[2,]))
deviation_simul_df=cbind(deviation_simul_df,rep(hseq,2))
deviation_simul_df=cbind(deviation_simul_df,rep(c(1,2),each=length(hseq)))
deviation_simul_df=as.data.frame(deviation_simul_df)
colnames(deviation_simul_df)=c("deviation","h","setting")
deviation_simul_df$setting=as.factor(deviation_simul_df$setting)
deviation_simul_df$deviation=as.numeric(deviation_simul_df$deviation)
deviation_simul_df$h=as.numeric(deviation_simul_df$h)


deviation_simul=ggplot(deviation_simul_df,aes(x=h,y=deviation,linetype=setting))+
  geom_line()+geom_point()+
  scale_linetype_manual(values=c("1"="dotted","2"="dashed"))+
  labs(x="Bandwidth h",y=expression(D(h)),linetype="Setting")+
  ggtitle("Deviation D(h) of RLCP in simulation settings")+
  theme(axis.title = element_text(size = 16), 
        axis.text = element_text(size=16),
        plot.title = element_text(hjust = 0.5,size=16), 
        legend.position = "bottom",
        legend.text = element_text(size=14))

deviation_simul

