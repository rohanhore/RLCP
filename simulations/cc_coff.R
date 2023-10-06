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
  Vcalib=abs(calib_data$Y-predict(model_lm,newdata=calib_data))
  
  #-------------test data------------------------------
  test_data=conditional_simulation(n=100,r=0,setting)
  Xtest=as.matrix(test_data[,-1])
  Vtest=abs(test_data$Y-predict(model_lm,newdata=test_data))
  
  #-----------evaluating the competing methods-----------
  RLCP_res=RLCP(Xcalib,Vcalib,Xtest,Vtest,"gaussian",h,alpha)
  calLCP_res=calLCP(Xcalib,Vcalib,Xtest,Vtest,"gaussian",h,alpha)
  baseLCP_res=baseLCP(Xcalib,Vcalib,Xtest,Vtest,"gaussian",h,alpha)
  
  #--------------cutoffs for score-------------------
  Vcutoff_RLCP=RLCP_res[,2]
  Vcutoff_RLCP[Vcutoff_RLCP==Inf]=max(Vcalib,Vtest)+0.2
  
  Vcutoff_calLCP=calLCP_res[,2]
  Vcutoff_calLCP[Vcutoff_calLCP==Inf]=max(Vcalib,Vtest)+0.2
  
  Vcutoff_baseLCP=baseLCP_res[,2]
  Vcutoff_baseLCP[Vcutoff_baseLCP==Inf]=max(Vcalib,Vtest)+0.2
  
  Vcutoffs=cbind(Vcutoff_RLCP,Vcutoff_calLCP,Vcutoff_baseLCP)
  
  #------------upper and lower bounds on response-----------
  
  pred_test=cbind(predict(model_lm,newdata=test_data),predict(model_lm,newdata=test_data),
                  predict(model_lm,newdata=test_data))
  upper_cutoffs=pred_test+Vcutoffs
  lower_cutoffs=pred_test-Vcutoffs
  
  return(cbind(upper_cutoffs,lower_cutoffs))
}
#--------------------------------------------------------------------
#-----------------------uni-variate experiments----------------------
#--------------------------------------------------------------------

d=1

numcores=detectCores()-1
cl=makeCluster(numcores)
registerDoParallel(cl)

#bandwidth choices
hseq=c(0.1,0.5,2.5)
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
  result=read.csv(paste0("../results/setting_1_cc_h_",h,".csv"),header=T)[,-1]
  result=as.data.frame(result)
  result_up=cbind(result[,1:3],0.5*xseq+oracle*qnorm(0.95,0,1))
  colnames(result_up)=c("RLCP","calLCP","baseLCP","oracle")
  result_up=stack(result_up)
  result_up=cbind(result_up,1)
  result_up=cbind(result_up,h)
  result_up=cbind(result_up,"U")
  result_up=cbind(result_up,rep(xseq,4))
  result_up=cbind(result_up,"bound")
  colnames(result_up)=c("cutoff","method","setting","h","end","covariate","ID")
  
  result_down=cbind(result[,4:6],0.5*xseq+oracle*qnorm(0.05,0,1))
  colnames(result_down)=c("RLCP","calLCP","baseLCP","oracle")
  result_down=stack(result_down)
  result_down=cbind(result_down,1)
  result_down=cbind(result_down,h)
  result_down=cbind(result_down,"L")
  result_down=cbind(result_down,rep(xseq,4))
  result_down=cbind(result_down,"bound")
  colnames(result_down)=c("cutoff","method","setting","h","end","covariate","ID")
  plot_result=rbind(plot_result,result_up,result_down)
}
for(i in 1:3){
  h=hseq[i]
  oracle=2*dnorm(xseq,0,1.5)
  result=read.csv(paste0("../results/setting_2_cc_h_",h,".csv"),header=T)[,-1]
  result=as.data.frame(result)
  result_up=cbind(result[,1:3],0.5*xseq+oracle*qnorm(0.95,0,1))
  colnames(result_up)=c("RLCP","calLCP","baseLCP","oracle")
  result_up=stack(result_up)
  result_up=cbind(result_up,2)
  result_up=cbind(result_up,h)
  result_up=cbind(result_up,"U")
  result_up=cbind(result_up,rep(xseq,4))
  result_up=cbind(result_up,"bound")
  colnames(result_up)=c("cutoff","method","setting","h","end","covariate","ID")
  
  result_down=cbind(result[,4:6],0.5*xseq+oracle*qnorm(0.05,0,1))
  colnames(result_down)=c("RLCP","calLCP","baseLCP","oracle")
  result_down=stack(result_down)
  result_down=cbind(result_down,2)
  result_down=cbind(result_down,h)
  result_down=cbind(result_down,"L")
  result_down=cbind(result_down,rep(xseq,4))
  result_down=cbind(result_down,"bound")
  colnames(result_down)=c("cutoff","method","setting","h","end","covariate","ID")
  plot_result=rbind(plot_result,result_up,result_down)
}

data=simulation(4000,d=1,1)
for(i in 1:3){
  h=hseq[i]
  result_scatter=data$Y
  result_scatter=cbind(data$Y,"RLCP")
  result_scatter=cbind(result_scatter,1)
  result_scatter=cbind(result_scatter,h)
  result_scatter=cbind(result_scatter,"U")
  result_scatter=cbind(result_scatter,data$X1)
  result_scatter=cbind(result_scatter,"scatter")
  colnames(result_scatter)=c("cutoff","method","setting","h","end","covariate","ID")
  
  plot_result=rbind(plot_result,result_scatter)
}

data=simulation(4000,d=1,2)
for(i in 1:3){
  h=hseq[i]
  result_scatter=data$Y
  result_scatter=cbind(data$Y,"RLCP")
  result_scatter=cbind(result_scatter,2)
  result_scatter=cbind(result_scatter,h)
  result_scatter=cbind(result_scatter,"U")
  result_scatter=cbind(result_scatter,data$X1)
  result_scatter=cbind(result_scatter,"scatter")
  colnames(result_scatter)=c("cutoff","method","setting","h","end","covariate","ID")
  
  plot_result=rbind(plot_result,result_scatter)
}

colnames(plot_result)=c("cutoff","method","setting","h","end","covariate","ID")
level_order=c('oracle','baseLCP','calLCP','RLCP')
plot_result$method=factor(plot_result$method,level=level_order)
plot_result$covariate=as.numeric(plot_result$covariate)
plot_result$cutoff=as.numeric(plot_result$cutoff)


conditional_plot=ggplot(plot_result[plot_result$ID %in% c("bound"),], aes(x = covariate, y = cutoff,group=interaction(method,end),
                                                                          color=method,linetype=method)) +
  geom_point(data=plot_result[plot_result$ID %in% c("scatter"),],aes(x=covariate,y=cutoff),shape=1,col="gray",alpha=0.3)+
  geom_line()+
  scale_color_manual(values=c("oracle"="black","RLCP"="maroon","calLCP"="gold3","baseLCP"="blue"))+
  scale_linetype_manual(values=c("oracle"="dashed",RLCP="solid","calLCP"="solid","baseLCP"="solid"),
                        guide="none")+
  xlim(-3.5,3.5)+
  facet_grid(setting~h ,labeller=label_bquote(rows = paste("Setting ", .(setting)),cols= h ==.(h)))+
  scale_y_continuous(n.breaks=6)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(color="Method",x="Feature X",y=" ")

pdf(file = "../results/simul_marginal_cutoffs.pdf",width = 7,height = 5) 
conditional_plot
dev.off()
