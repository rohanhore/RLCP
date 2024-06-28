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
  
  #----------local coverage-----------------------
  centers=seq(-2.5,2.5,by=0.25);radii=0.4
  local_ball_id=(abs(matrix(rep(Xtest,length(centers)),nrow=ntest)
                     -matrix(rep(centers,each=ntest),nrow=ntest))<=radii)
  ball_size=apply(local_ball_id,2,sum)
  coverage1=as.vector(t(RLCP_res[,1])%*% local_ball_id/ball_size)
  
  return(coverage1)
}
#--------------------------------------------------------------------
#-----------------------uni-variate experiments----------------------
#--------------------------------------------------------------------

d=1;alpha=0.1

numcores=detectCores()-1
cl=makeCluster(numcores)
registerDoParallel(cl)

#bandwidth choices
hseq=c(0.025,0.05)

#setting 1
setting=1;ntrain=2000;ncalib=2000;ntest=2000
nrep=50
xseq=seq(-2.5,2.5,by=0.25)

RLCP_local_cov.1=matrix(0,nrow=length(hseq),ncol=21)
for(j in 1:length(hseq)){
  h=hseq[j]
  print(j)
  result_h=foreach(k=1:nrep,.combine="rbind",
                   .packages = c("MASS","parallel","doParallel",
                                 "foreach","mvtnorm","quantreg")) %dopar% {
                                   local_ball_cov_onerep(k,d,setting)
                                 }
  RLCP_local_cov.1[j,]=colMeans(result_h)
}

#setting 2
setting=2;ntrain=2000;ncalib=2000;ntest=2000
nrep=50

RLCP_local_cov.2=matrix(0,nrow=length(hseq),ncol=21)
for(j in 1:length(hseq)){
  h=hseq[j]
  print(j)
  result_h=foreach(k=1:nrep,.combine="rbind",
                   .packages = c("MASS","parallel","doParallel",
                                 "foreach","mvtnorm","quantreg")) %dopar% {
                                   local_ball_cov_onerep(k,d,setting)
                                 }
  RLCP_local_cov.2[j,]=colMeans(result_h)
}

gibbs_coverage=read.csv("../results/gibbs_et_al_results.csv",header=F,sep=",")

#-----------------------------------------------------
#----------------Visualization------------------------
#-----------------------------------------------------

##plotting local coverage
result=as.data.frame(t(rbind(gibbs_coverage,RLCP_local_cov.1,RLCP_local_cov.2)))
plot_result=matrix(0,nrow=4*21*length(hseq),ncol=5)
plot_result[,1]=stack(result)$values
plot_result[,2]=rep(xseq,2*2*length(hseq))
plot_result[,3]=rep(rep(c("high","low"),each=21),2*2)
plot_result[,4]=rep(c("gibbs_et_al","RLCP"),each=21*length(hseq)*2)
plot_result[,5]=rep(rep(1:2,each=21*2),2)

plot_result=as.data.frame(plot_result)
colnames(plot_result)=c("coverage","center","adaptivity","method","setting")
plot_result$coverage=as.numeric(plot_result$coverage)
plot_result$setting=as.factor(plot_result$setting)
plot_result$center=as.numeric(plot_result$center)
level_order=c('gibbs_et_al','RLCP')



pdf(file = "../results/figures/local_coverage_RLCP_gibbs.pdf",width = 8,height=5)

cols<-c("RLCP"="maroon","gibbs_et_al"="seagreen")
ggplot(plot_result, aes(x = center, y = coverage, col = method)) +
  geom_line() + 
  geom_point()+
  geom_hline(yintercept=0.9,linetype="dashed")+
  scale_color_manual(values=c("RLCP"="maroon","gibbs_et_al"="seagreen"),
                     labels=c("RLCP"="RLCP","gibbs_et_al"="Gibbs et al."))+
  #scale_linetype_manual(values=c("high"="solid","low"="longdash"))+
  facet_grid(adaptivity~setting,labeller=label_bquote(
    rows=paste(.(adaptivity)," adaptivity" ),cols=paste("Setting ", .(setting))))+
  labs(x=expression(x[0]))+
  scale_x_continuous(limits = c(-2.5, 2.5),breaks = c(-2,0,2))+
  #ylim(0.89,0.91)+
  xlab("Feature X")+
  theme(legend.position = "bottom",
        panel.spacing = unit(0.5,"cm",data=NULL),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 15),
        axis.text=element_text(size=14),
        axis.title=element_text(size=18))+
  guides(SET = guide_legend(byrow = TRUE))

dev.off()
