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
  
  pred_test=matrix(rep(predict(model_lm,newdata=test_data),3),nrow=dim(test_data)[1])
  
  upper_thresholds=pred_test+score_thresholds
  lower_thresholds=pred_test-score_thresholds
  
  return(cbind(upper_thresholds,lower_thresholds))
}

#----------------------------------------------------------
#-----------computing empirical local coverage-------------
#----------------------------------------------------------
local_ball_cov_onerep=function(k,d,setting){
  
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
  test_data=simulation(ntest,d,setting)
  Xtest=as.matrix(test_data[,-1])
  scores_test=abs(test_data$Y-predict(model_lm,newdata=test_data))
  
  #-----------computing the PI for competing methods-----------
  RLCP_res=RLCP(Xcalib,scores_calib,Xtest,scores_test,"gaussian",h,alpha)
  calLCP_res=calLCP(Xcalib,scores_calib,Xtest,scores_test,"gaussian",h,alpha)
  baseLCP_res=baseLCP(Xcalib,scores_calib,Xtest,scores_test,"gaussian",h,alpha)
  
  #----------local coverage-----------------------
  centers=seq(-2.5,2.5,by=0.25);radii=0.4
  local_ball_id=(abs(matrix(rep(Xtest,length(centers)),nrow=ntest)
                     -matrix(rep(centers,each=ntest),nrow=ntest))<=radii)
  ball_size=apply(local_ball_id,2,sum)
  coverage1=as.vector(t(RLCP_res[,1])%*% local_ball_id/ball_size)
  coverage2=as.vector(t(calLCP_res[,1])%*%local_ball_id/ball_size)
  coverage3=as.vector(t(baseLCP_res[,1])%*%local_ball_id/ball_size)
  
  
  return(list(RLCP_cov=coverage1,
              calLCP_cov=coverage2,
              baseLCP_cov=coverage3))
}
#--------------------------------------------------------------------
#-----------------------uni-variate experiments----------------------
#--------------------------------------------------------------------

d=1;alpha=0.1

numcores=detectCores()-1
cl=makeCluster(numcores)
registerDoParallel(cl)

#bandwidth choices
hseq=c(0.1,0.2,0.4,0.8,1.6)
xseq=seq(-3,3,by=0.01)

#setting 1
setting=1;ntrain=2000;ncalib=2000;ntest=2000
nrep=50
for(j in 1:length(hseq)){
  h=hseq[j]
  print(j)
  result_h=foreach(k=1:nrep,.combine="+",
                   .packages = c("MASS","parallel","doParallel",
                                 "foreach","mvtnorm","quantreg")) %dopar% {
                                   cc_cutoff_onerep(k,d,setting)
                                 }
  result_h=result_h/nrep
  result_h=apply(result_h,2,FUN=function(x) smooth.spline(xseq,x,spar=0.3)$y)
  write.csv(result_h,paste0("../results/setting_1_cc_h_",h,".csv"))
}

comb=function(x,y){return(list(rbind(x[[1]],y[[1]]),rbind(x[[2]],y[[2]]),rbind(x[[3]],y[[3]])))}
RLCP_local_cov.1=calLCP_local_cov.1=baseLCP_local_cov.1=matrix(0,nrow=length(hseq)*2,ncol=21)
for(j in 1:length(hseq)){
  h=hseq[j]
  print(h)
  result_h=foreach(k=1:nrep,.combine=comb,
                   .packages = c("MASS","parallel","doParallel",
                                 "foreach","mvtnorm")) %dopar% {
                                   local_ball_cov_onerep(k,d,setting)
                                 }
  RLCP_local_cov.1[j,]=colMeans(result_h[[1]])
  RLCP_local_cov.1[length(hseq)+j,]=apply(result_h[[1]],2,sd)/sqrt(nrep)
  
  calLCP_local_cov.1[j,]=colMeans(result_h[[2]])
  calLCP_local_cov.1[length(hseq)+j,]=apply(result_h[[2]],2,sd)/sqrt(nrep)
  
  baseLCP_local_cov.1[j,]=colMeans(result_h[[3]])
  baseLCP_local_cov.1[length(hseq)+j,]=apply(result_h[[3]],2,sd)/sqrt(nrep)
}

#setting 2
setting=2;ntrain=2000;ncalib=2000;ntest=2000
nrep=50
for(j in 1:length(hseq)){
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

ntrain=ncalib=ntest=1000;nrep=30
RLCP_local_cov.2=calLCP_local_cov.2=baseLCP_local_cov.2=matrix(0,nrow=length(hseq)*2,ncol=21)
for(j in 1:length(hseq)){
  h=hseq[j]
  h=0.1
  print(h)
  result_h=foreach(k=1:nrep,.combine=comb,
                   .packages = c("MASS","parallel","doParallel",
                                 "foreach","mvtnorm")) %dopar% {
                                   local_ball_cov_onerep(k,d,setting)
                                 }
  RLCP_local_cov.2[j,]=colMeans(result_h[[1]])
  RLCP_local_cov.2[length(hseq)+j,]=apply(result_h[[1]],2,sd)/sqrt(nrep)
  
  calLCP_local_cov.2[j,]=colMeans(result_h[[2]])
  calLCP_local_cov.2[length(hseq)+j,]=apply(result_h[[2]],2,sd)/sqrt(nrep)
  
  baseLCP_local_cov.2[j,]=colMeans(result_h[[3]])
  baseLCP_local_cov.2[length(hseq)+j,]=apply(result_h[[3]],2,sd)/sqrt(nrep)
}
stopCluster(cl)

#-----------------------------------------------------
#----------------Visualization------------------------
#-----------------------------------------------------

##plotting average intervals

plot_result=data.frame()
xseq=seq(-3,3,by=0.01)
for(i in 1:length(hseq)){
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
for(i in 1:length(hseq)){
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
for(i in 1:length(hseq)){
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
for(i in 1:length(hseq)){
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

##plotting local coverage

result.1=cbind(t(rbind(RLCP_local_cov.1[1:5,],calLCP_local_cov.1[1:5,],baseLCP_local_cov.1[1:5,])),
             t(rbind(RLCP_local_cov.1[6:10,],calLCP_local_cov.1[6:10,],baseLCP_local_cov.1[6:10,])))
result.1=as.data.frame(result.1)

plot_result.1=matrix(0,nrow=3*21*length(hseq),ncol=5)
plot_result.1[,1]=stack(result.1[,1:15])$values
plot_result.1[,2]=rep(seq(-2.5,2.5,by=0.25),3*length(hseq))
plot_result.1[,3]=rep(rep(hseq,each=21),3)
plot_result.1[,4]=rep(c("RLCP","calLCP","baseLCP"),each=21*length(hseq))
plot_result.1[,5]=stack(result.1[,16:30])$values

plot_result.1=as.data.frame(plot_result.1)
colnames(plot_result.1)=c("coverage","center","h","method","se")
plot_result.1$coverage=as.numeric(plot_result.1$coverage)
plot_result.1$se=as.numeric(plot_result.1$se)
plot_result.1$h=as.numeric(plot_result.1$h)
plot_result.1$center=as.numeric(plot_result.1$center)
level_order=c('baseLCP','calLCP','RLCP')

result.2=cbind(t(rbind(RLCP_local_cov.2[1:5,],calLCP_local_cov.2[1:5,],baseLCP_local_cov.2[1:5,])),
               t(rbind(RLCP_local_cov.2[6:10,],calLCP_local_cov.2[6:10,],baseLCP_local_cov.2[6:10,])))
result.2=as.data.frame(result.2)

plot_result.2=matrix(0,nrow=3*21*length(hseq),ncol=5)
plot_result.2[,1]=stack(result.2[,1:15])$values
plot_result.2[,2]=rep(seq(-2.5,2.5,by=0.25),3*length(hseq))
plot_result.2[,3]=rep(rep(hseq,each=21),3)
plot_result.2[,4]=rep(c("RLCP","calLCP","baseLCP"),each=21*length(hseq))
plot_result.2[,5]=stack(result.2[,16:30])$values

plot_result.2=as.data.frame(plot_result.2)
colnames(plot_result.2)=c("coverage","center","h","method","se")
plot_result.2$coverage=as.numeric(plot_result.2$coverage)
plot_result.2$se=as.numeric(plot_result.2$se)
plot_result.2$h=as.numeric(plot_result.2$h)
plot_result.2$center=as.numeric(plot_result.2$center)
level_order=c('baseLCP','calLCP','RLCP')





cols<-c("RLCP"="maroon","oracle"="black","calLCP"="gold3","baseLCP"="blue")
linetypes<-c("RLCP"="solid",oracle="dotted","calLCP"="solid","baseLCP"="solid")

plot_1.1=ggplot(plot_result[plot_result$ID %in% c("threshold")&plot_result$setting==1,]) +
  geom_point(data=plot_result[plot_result$ID %in% c("scatter")&plot_result$setting==1,],aes(x=covariate,y=RLCP_upper),shape=1,col="gray",alpha=0.3)+
  facet_grid(.~h ,labeller=label_bquote(rows = paste("Setting ", .(setting)),cols= h ==.(h)))+
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

plot_1.2=ggplot(plot_result[plot_result$ID %in% c("threshold")&plot_result$setting==2,]) +
  geom_point(data=plot_result[plot_result$ID %in% c("scatter")&plot_result$setting==2,],aes(x=covariate,y=RLCP_upper),shape=1,col="gray",alpha=0.3)+
  facet_grid(.~h ,labeller=label_bquote(rows = paste("Setting ", .(setting)),cols= h ==.(h)))+
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

plot_2.1=ggplot(plot_result.1, aes(x = center, y = coverage, col = method)) +
  geom_line() + 
  geom_point()+
  #geom_errorbar(aes(ymin=coverage-se, ymax=coverage+se), width=.05)+
  geom_hline(yintercept=0.9,linetype="dashed")+
  scale_color_manual(values=c("RLCP"="maroon","calLCP"="gold3","baseLCP"="blue"))+
  facet_grid(.~h,labeller=label_bquote(cols=h ==.(h)))+
  labs(x=expression(x[0]))+
  scale_x_continuous(breaks = c(-2,0,2))+
  xlab("Feature X")+
  theme(#legend.position = "none",
        legend.text = element_text(size=16),
        panel.spacing = unit(0.5,"cm",data=NULL),
        strip.text = element_text(size = 18),
        axis.text=element_text(size=18),
        axis.title=element_text(size=18))+
  guides(SET = guide_legend(byrow = TRUE))


plot_2.2=ggplot(plot_result.2, aes(x = center, y = coverage, col = method)) +
  geom_line() + 
  geom_point()+
  #geom_errorbar(aes(ymin=coverage-se, ymax=coverage+se), width=.05)+
  geom_hline(yintercept=0.9,linetype="dashed")+
  scale_color_manual(values=c("RLCP"="maroon","calLCP"="gold3","baseLCP"="blue"))+
  facet_grid(.~h,labeller=label_bquote(cols=h ==.(h)))+
  labs(x=expression(x[0]))+
  scale_x_continuous(breaks = c(-2,0,2))+
  xlab("Feature X")+
  theme(legend.position = "none",
        panel.spacing = unit(0.5,"cm",data=NULL),
        strip.text = element_text(size = 18),
        axis.text=element_text(size=18),
        axis.title=element_text(size=18))+
  guides(SET = guide_legend(byrow = TRUE))

pdf(file = "../results/figures/local_coverage_setting_1.pdf",width = 12,height = 8)
ggpubr::ggarrange(plot_1.1,plot_2.1,nrow=2,common.legend = TRUE,legend = "right")
dev.off()

pdf(file = "../results/figures/local_coverage_setting_2.pdf",width = 12,height = 8)
ggpubr::ggarrange(plot_1.2,plot_2.2,nrow=2,common.legend = TRUE,legend = "right")
dev.off()

