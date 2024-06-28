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
#-----------------------------------------------------------
#------------------multivariate experiments-----------------
#-----------------------------------------------------------
#dimension choices
dseq=c(5*1:6)
#---set prefixed effective sample size----
eff_size=50

#---finding optimum bandwidth choices of RLCP and calLCP---
optimum_RLCP_bandwidth=optimum_calLCP_bandwidth=rep(0,length(dseq))
for(i in 1:length(dseq)){
  if(i==1){h_min_RLCP=h_min_calLCP=0.02}
  else{h_min_RLCP=optimum_RLCP_bandwidth[i-1];h_min_calLCP=optimum_calLCP_bandwidth[i-1]}
  Xtrain=as.matrix(simulation(3000,dseq[i],3)[,-1]) 
  optimum_RLCP_bandwidth[i]=optimum_RLCP_h(Xtrain,"gaussian",h_min_RLCP,eff_size)
  optimum_calLCP_bandwidth[i]=optimum_calLCP_h(Xtrain,"gaussian",h_min_calLCP,eff_size)
}

#-----------------------------------------------------------
#-----------computing marginal coverage on sets-------------
#-----------------------------------------------------------
local_cov_optimized=function(k,setting,j,radii){
  d=dseq[j]
  
  #----------training data and learning score------------
  train_data=simulation(3000,d,setting)
  Xnames=paste("X", 1:d, sep="")
  formula=as.formula(paste("Y ~ ", paste(Xnames, collapse= "+")))
  model_lm=lm(formula,data=train_data)
  
  #---------consider the optimum bandwidth--------------
  Xtrain=as.matrix(train_data[,-1])
  h_opt_RLCP=optimum_RLCP_bandwidth[j]
  h_opt_calLCP=optimum_calLCP_bandwidth[j]
  
  #----------calibration data--------------------------
  calib_data=simulation(3000,d,setting)
  Xcalib=as.matrix(calib_data[,-1])
  scores_calib=abs(calib_data$Y-predict(model_lm,newdata=calib_data))
  
  #-------------test data------------------------------
  test_data=simulation(n,d,setting)
  Xtest=as.matrix(test_data[,-1])
  scores_test=abs(test_data$Y-predict(model_lm,newdata=test_data))
  
  #-----------evaluating the competing methods-----------
  RLCP_res=RLCP(Xcalib,scores_calib,Xtest,scores_test,"gaussian",h_opt_RLCP,alpha)
  calLCP_res=calLCP(Xcalib,scores_calib,Xtest,scores_test,"gaussian",h_opt_calLCP,alpha)
  
  #----------coverage-----------------------
  coverage1=coverage2=rep(0,n_ball)
  for(i in 1:n_ball){
    ID=ind_set(Xtest[,1:3],ball_centers[i,],radii)
    
    coverage1[i]=mean(RLCP_res[ID,1])
    coverage2[i]=mean(calLCP_res[ID,1])
  }
  
  coverages=cbind(coverage1,coverage2)
  return(coverages)
}


#------------------------------------------------------------
#-----------------------experimental results-----------------
#------------------------------------------------------------

numcores=detectCores()-1
cl=makeCluster(numcores)
registerDoParallel(cl)


setting=3
nrep=30;alpha=0.1
dseq=c(5*1:6)
alpha=0.1
n=5000

#setting up the grid
ind_set=function(X,center,radii){which(apply(X,1,FUN=function(x) max(abs(x-center))<=radii)==1)}
grid_centers=seq(-3+0.5,3-0.5,by=1)
grid_length=length(grid_centers)
ball_centers=matrix(0,nrow=grid_length^3,3)
ball_centers[,1]=rep(grid_centers,each=grid_length^2)
ball_centers[,2]=rep(rep(grid_centers,each=grid_length),grid_length)
ball_centers[,3]=rep(grid_centers,grid_length^2)
radii=1;n_ball=grid_length^3

RLCP_local_cov=calLCP_local_cov=matrix(0,nrow=length(dseq),ncol=n_ball)
for(j in 1:length(dseq)){
  d=dseq[j]
  result_j=foreach(k=1:nrep,.combine="+",
                   .packages = c("MASS","parallel","doParallel",
                                 "foreach","mvtnorm")) %dopar% {
                                   local_cov_optimized(k,setting,j,radii)
                                 }
  result_j=result_j/nrep
  RLCP_local_cov[j,]=result_j[,1]
  calLCP_local_cov[j,]=result_j[,2]
}

##plotting
plot_result=matrix(0,nrow=n_ball*2,ncol=2)
histogram_plots=list()
for(i in 1:length(dseq)){
  plot_result[,1]=c(RLCP_local_cov[i,],calLCP_local_cov[i,])
  plot_result[,2]=rep(c("RLCP","calLCP"),each=n_ball)
  plot_result=as.data.frame(plot_result)
  colnames(plot_result)=c("coverage","method")
  plot_result$coverage=as.numeric(plot_result$coverage)
  plot_result$method=as.factor(plot_result$method)
  histogram_plots[[i]]=ggplot(plot_result, aes(x=coverage,fill=method))+
    geom_density(alpha=0.6)+
    scale_x_continuous(n.breaks = 5)+
    scale_fill_manual(values=c("RLCP"="maroon","calLCP"="gold3"))+
    ggtitle(paste("d =",dseq[i]))+
    theme(plot.title = element_text(hjust = 0.5,size=18),
          axis.text=element_text(size=12),
          axis.text.x=element_text(hjust=0.5),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title = element_text(size=16),
          plot.margin = unit(c(0.15, 0.15, 0.15, 0.15), 
                             "inches"),
          strip.text = element_text(size=16),
          legend.text = element_text(size=18))+
    ylab(" ")
}

pdf(file = "../results/figures/local_coverage_highd.pdf",width = 11,height = 8)
ggpubr::ggarrange(histogram_plots[[1]],histogram_plots[[2]],histogram_plots[[3]],
                  histogram_plots[[4]],histogram_plots[[5]],histogram_plots[[6]],
                  nrow=3,ncol=2,common.legend=TRUE,legend="bottom")
dev.off()
