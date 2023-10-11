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
#---------------------------------------------------
#-----------computing average-case PI---------------
#---------------------------------------------------
cc_cutoff_onerep=function(k,d,setting){
  
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
  test_data=conditional_simulation(n=100,setting)
  Xtest=as.matrix(test_data[,-1])
  scores_test=abs(test_data$Y-predict(model_lm,newdata=test_data))
  
  #-----------evaluating the competing methods-----------
  RLCP_res=RLCP(Xcalib,scores_calib,Xtest,scores_test,"gaussian",h,alpha)
  calLCP_res=calLCP(Xcalib,scores_calib,Xtest,scores_test,"gaussian",h,alpha)
  baseLCP_res=baseLCP(Xcalib,scores_calib,Xtest,scores_test,"gaussian",h,alpha)
  
  #--------------cutoffs for score-------------------
  score_threshold_RLCP=RLCP_res[,2]
  score_threshold_RLCP[score_threshold_RLCP==Inf]=max(scores_calib,scores_test)+0.2
  
  score_threshold_calLCP=calLCP_res[,2]
  score_threshold_calLCP[score_threshold_calLCP==Inf]=max(scores_calib,scores_test)+0.2
  
  score_threshold_baseLCP=baseLCP_res[,2]
  score_threshold_baseLCP[score_threshold_baseLCP==Inf]=max(scores_calib,scores_test)+0.2
  
  score_thresholds=cbind(score_threshold_RLCP,score_threshold_calLCP,score_threshold_baseLCP)
  
  #------------upper and lower bounds on response-----------
  
  pred_test=cbind(predict(model_lm,newdata=test_data),predict(model_lm,newdata=test_data),
                  predict(model_lm,newdata=test_data))
  upper_thresholds=pred_test+score_thresholds
  lower_thresholds=pred_test-score_thresholds
  
  return(cbind(upper_thresholds,lower_thresholds))
}
#--------------------------------------------------------------------
#-----------------------uni-variate experiments----------------------
#--------------------------------------------------------------------

d=1

numcores=detectCores()-1
cl=makeCluster(numcores)
registerDoParallel(cl)

#bandwidth choices
hseq=c(0.5,1,1.5)
xseq=seq(-3,3,by=0.01)

#setting 1
setting=1;ntrain=2000;ncalib=2000
nrep=100
for(j in 1:3){
  h=hseq[j]
  result_h=foreach(k=1:nrep,.combine="+",
                   .packages = c("MASS","parallel","doParallel",
                                 "foreach","mvtnorm","quantreg")) %dopar% {
                                   cc_cutoff_onerep(k,d,setting)
                                 }
  result_h=result_h/nrep
  result_h=apply(result_h,2,FUN=function(x) smooth.spline(xseq,x,spar=0.3)$y)
  write.csv(result_h,paste0("../results/setting_1_cc_h_",h,".csv"))
}

#setting 2
setting=2;ntrain=2000;ncalib=2000
nrep=100
for(j in 1:3){
  h=hseq[j]
  result_h=foreach(k=1:nrep,.combine="+",
                   .packages = c("MASS","parallel","doParallel",
                                 "foreach","mvtnorm","quantreg")) %dopar% {
                                   cc_cutoff_onerep(k,d,setting)
                                 }
  result_h=result_h/nrep
  result_h=apply(result_h,2,FUN=function(x) smooth.spline(xseq,x,spar=0.3)$y)
  write.csv(result_h,paste0("../results/setting_2_cc_h_",h,".csv"))
}
stopCluster(cl)

#-----------------------------------------------------
#----------------Visualization------------------------
#-----------------------------------------------------

plot_result=data.frame()
xseq=seq(-3,3,by=0.01)
for(i in 1:3){
  h=hseq[i]
  oracle=abs(sin(xseq))
  h_result=read.csv(paste0("../results/setting_1_cc_h_",h,".csv"),header=T)[,-1]
  h_result=cbind(h_result,0.5*xseq+oracle*qnorm(0.95,0,1))
  h_result=cbind(h_result,0.5*xseq+oracle*qnorm(0.05,0,1))
  h_result=cbind(h_result,1)
  h_result=cbind(h_result,hseq[i])
  h_result=cbind(h_result,xseq)
  h_result=cbind(h_result,"threshold")
  colnames(h_result)=c("RLCP_upper","calLCP_upper","baseLCP_upper",
                       "RLCP_lower","calLCP_lower","baseLCP_lower",
                       "oracle_upper","oracle_lower",
                       "setting","h","covariate","ID")
  plot_result=rbind(plot_result,h_result)
}
for(i in 1:3){
  h=hseq[i]
  oracle=2*dnorm(xseq,0,1.5)
  h_result=read.csv(paste0("../results/setting_2_cc_h_",h,".csv"),header=T)[,-1]
  h_result=cbind(h_result,0.5*xseq+oracle*qnorm(0.95,0,1))
  h_result=cbind(h_result,0.5*xseq+oracle*qnorm(0.05,0,1))
  h_result=cbind(h_result,2)
  h_result=cbind(h_result,hseq[i])
  h_result=cbind(h_result,xseq)
  h_result=cbind(h_result,"threshold")
  colnames(h_result)=c("RLCP_upper","calLCP_upper","baseLCP_upper",
                       "RLCP_lower","calLCP_lower","baseLCP_lower",
                       "oracle_upper","oracle_lower",
                       "setting","h","covariate","ID")
  plot_result=rbind(plot_result,h_result)
}
colnames(plot_result)=c("RLCP_upper","calLCP_upper","baseLCP_upper",
                        "RLCP_lower","calLCP_lower","baseLCP_lower",
                        "oracle_upper","oracle_lower",
                        "setting","h","covariate","ID")

data=simulation(4000,d=1,1)
result_scatter=matrix(0,ncol=12,nrow=4000)
for(i in 1:3){
  result_scatter[,1]=data$Y
  result_scatter[,9]=1
  result_scatter[,10]=hseq[i]
  result_scatter[,11]=data$X1
  result_scatter[,12]="scatter"
  colnames(result_scatter)=c("RLCP_upper","calLCP_upper","baseLCP_upper",
                             "RLCP_lower","calLCP_lower","baseLCP_lower",
                             "oracle_upper","oracle_lower",
                             "setting","h","covariate","ID")
  plot_result=rbind(plot_result,result_scatter)
}

data=simulation(4000,d=1,2)
result_scatter=matrix(0,ncol=12,nrow=4000)
for(i in 1:3){
  result_scatter[,1]=data$Y
  result_scatter[,9]=2
  result_scatter[,10]=hseq[i]
  result_scatter[,11]=data$X1
  result_scatter[,12]="scatter"
  colnames(result_scatter)=c("RLCP_upper","calLCP_upper","baseLCP_upper",
                             "RLCP_lower","calLCP_lower","baseLCP_lower",
                             "oracle_upper","oracle_lower",
                             "setting","h","covariate","ID")
  plot_result=rbind(plot_result,result_scatter)
}

colnames(plot_result)=c("RLCP_upper","calLCP_upper","baseLCP_upper",
                        "RLCP_lower","calLCP_lower","baseLCP_lower",
                        "oracle_upper","oracle_lower",
                        "setting","h","covariate","ID")

plot_result$oracle_upper=as.numeric(plot_result$oracle_upper);plot_result$oracle_lower=as.numeric(plot_result$oracle_lower)
plot_result$calLCP_upper=as.numeric(plot_result$calLCP_upper);plot_result$calLCP_lower=as.numeric(plot_result$calLCP_lower)
plot_result$baseLCP_upper=as.numeric(plot_result$baseLCP_upper);plot_result$baseLCP_lower=as.numeric(plot_result$baseLCP_lower)
plot_result$RLCP_upper=as.numeric(plot_result$RLCP_upper);plot_result$RLCP_lower=as.numeric(plot_result$RLCP_lower)
plot_result$covariate=as.numeric(plot_result$covariate)


pdf(file = "../results/figures/simul_marginal_cutoffs.pdf",width = 12,height = 8)

cols<-c("RLCP"="maroon","oracle"="black","calLCP"="gold3","baseLCP"="blue")
linetypes<-c("RLCP"="solid",oracle="dotted","calLCP"="solid","baseLCP"="solid")
ggplot(plot_result[plot_result$ID %in% c("threshold"),]) +
  geom_point(data=plot_result[plot_result$ID %in% c("scatter"),],aes(x=covariate,y=RLCP_upper),shape=1,col="gray",alpha=0.3)+
  facet_grid(setting~h ,labeller=label_bquote(rows = paste("Setting ", .(setting)),cols= h ==.(h)))+
  geom_line(aes(x=covariate,y=oracle_upper,color="oracle",linetype="oracle"),lwd=1)+
  geom_line(aes(x=covariate,y=oracle_lower,color="oracle",linetype="oracle"),lwd=1)+
  geom_line(aes(x=covariate,y=calLCP_upper,color="calLCP",linetype="calLCP"),lwd=1)+
  geom_line(aes(x=covariate,y=calLCP_lower,color="calLCP",linetype="calLCP"),lwd=1)+
  geom_line(aes(x=covariate,y=baseLCP_upper,color="baseLCP",linetype="baseLCP"),lwd=1)+
  geom_line(aes(x=covariate,y=baseLCP_lower,color="baseLCP",linetype="baseLCP"),lwd=1)+
  geom_line(aes(x=covariate,y=RLCP_upper,color="RLCP",linetype="RLCP"),lwd=1)+
  geom_line(aes(x=covariate,y=RLCP_lower,color="RLCP",linetype="RLCP"),lwd=1)+
  xlim(-3.5,3.5)+
  scale_color_manual(name="Method",values=cols,breaks=c("oracle","baseLCP","calLCP","RLCP"))+
  scale_linetype_manual(name="Method",values=linetypes,breaks=c("oracle","baseLCP","calLCP","RLCP"))+
  scale_y_continuous(n.breaks=6)+
  labs(color="Method",x="Feature X",y=" ")+
  theme_bw()+
  theme(strip.text = element_text(size = 18),legend.text=element_text(size=18),
        axis.text=element_text(size=18),axis.title=element_text(size=18),
        legend.title=element_text(size=18),legend.spacing.y = unit(0.3,'cm'))+
  guides(color = guide_legend(byrow = TRUE))
dev.off()

