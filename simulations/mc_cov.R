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
#-----------computing marginal coverage-------------
#---------------------------------------------------
mc_onerep=function(k,d,setting){
  
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
  
  #----------coverage-----------------------
  coverage1=mean(RLCP_res[,1])
  coverage2=mean(calLCP_res[,1])
  coverage3=mean(baseLCP_res[,1])
  
  #---------percentage of infinite cutoff----------
  Inf_RLCP=which(RLCP_res[,2]=="Inf")
  pct_RLCP=length(Inf_RLCP)/ntest
  
  Inf_calLCP=which(calLCP_res[,2]=="Inf")
  pct_calLCP=length(Inf_calLCP)/ntest
  
  Inf_baseLCP=which(baseLCP_res[,2]=="Inf")
  pct_baseLCP=length(Inf_baseLCP)/ntest
  
  return(c(coverage1,coverage2,coverage3,pct_RLCP,pct_calLCP,pct_baseLCP))
}

#--------------------------------------------------------------------
#-----------------------uni-variate experiments----------------------
#--------------------------------------------------------------------

d=1;alpha=0.1

numcores=detectCores()-1
cl=makeCluster(numcores)
registerDoParallel(cl)

#bandwidth choices
hseq=c(0.5,1,1.5)

#setting 1
setting=1;ntrain=2000;ncalib=2000;ntest=2000
#number of repetitions
nrep=100

for(j in 1:3){
  h=hseq[j]
  print(h)
  result_h=foreach(k=1:nrep,.combine=rbind,
                   .packages = c("MASS","parallel","doParallel",
                                 "foreach","mvtnorm")) %dopar% {
                                   mc_onerep(k,d,setting)
                                 }
  print(colMeans(result_h[,1:3]))
  write.csv(result_h,paste0("../results/setting_1_mc_",h,".csv"))
}

#setting 2
setting=2;ntrain=2000;ncalib=2000;ntest=2000
#number of repetitions
nrep=100

for(j in 1:3){
  h=hseq[j]
  result_h=foreach(k=1:nrep,.combine=rbind,
                   .packages = c("MASS","parallel","doParallel",
                                 "foreach","mvtnorm")) %dopar% {
                                   mc_onerep(k,d,setting)
                                 }
  print(colMeans(result_h[,1:3]))
  write.csv(result_h,paste0("../results/setting_2_mc_",h,".csv"))
}

#--------------------------------------------------------------------
#-------------------------Visualization------------------------------
#--------------------------------------------------------------------
plot_result=data.frame()
for(i in 1:3){
  h=hseq[i]
  result=read.csv(paste0("../results/setting_1_mc_",h,".csv"),header=T)
  result=as.data.frame(result)
  result=result[,2:4]
  colnames(result)=c("RLCP","calLCP","baseLCP")
  result=stack(result)
  result=cbind(result,1)
  result=cbind(result,h)
  colnames(result)=c("coverage","method","setting","h")
  plot_result=rbind(plot_result,result)
}
for(i in 1:3){
  h=hseq[i]
  result=read.csv(paste0("../results/setting_2_mc_",h,".csv"),header=T)
  result=as.data.frame(result)
  result=result[,2:4]
  colnames(result)=c("RLCP","calLCP","baseLCP")
  result=stack(result)
  result=cbind(result,2)
  result=cbind(result,h)
  colnames(result)=c("coverage","method","setting","h")
  plot_result=rbind(plot_result,result)
}

colnames(plot_result)=c("coverage","method","setting","h")
level_order=c('baseLCP','calLCP','RLCP')
plot_result$method=factor(plot_result$method,level=level_order)

marginal_plot=ggplot(plot_result, aes(x = method , y = coverage,fill=method)) +
  geom_boxplot() + 
  geom_hline(yintercept=0.9,linetype="dashed")+
  scale_fill_manual(values=c("baseLCP"="blue","calLCP"="gold3","RLCP"="maroon"))+
  facet_grid(setting~h,labeller=label_bquote(rows = paste("Setting ", .(setting)),cols= h ==.(h)))+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        strip.text = element_text(size = 15),
        axis.text=element_text(size=15),
        axis.title=element_text(size=15))+
  labs(x="Method",y="Coverage")


pdf(file = "../results/figures/simul_marginal_coverage.pdf",width = 10,height = 6) 
marginal_plot
dev.off()

