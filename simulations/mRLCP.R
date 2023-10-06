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

#------------------------------------------------------------
#-----------computing marginal coverage of m-rLCP------------
#------------------------------------------------------------
mrand_onerep=function(k,d,setting){
  
  #----------training data and learning score------------
  train_data=simulation(ntrain,d,setting)
  Xnames=paste("X", 1:d, sep="")
  formula=as.formula(paste("Y ~ ", paste(Xnames, collapse= "+")))
  model_lm=lm(formula,data=train_data)
  
  #----------calibration data--------------------------
  calib_data=simulation(ncalib,d,setting)
  Xcalib=as.matrix(calib_data[,-1])
  Vcalib=abs(calib_data$Y-predict(model_lm,newdata=calib_data))
  plot(calib_data$X1,Vcalib)
  
  #-------------test data------------------------------
  test_data=simulation(ntest,d,setting)
  Xtest=as.matrix(test_data[,-1])
  Vtest=abs(test_data$Y-predict(model_lm,newdata=test_data))
  
  #-----------evaluating the competing methods-----------
  m_rLCP_res=m_randomized_LCP(Xcalib,Vcalib,Xtest,Vtest,"gaussian",h,alpha,100)
  rLCP_res=randomized_LCP(Xcalib,Vcalib,Xtest,Vtest,"gaussian",h,alpha)
  
  #----------coverage-----------------------
  coverage1=mean(m_rLCP_res)
  coverage2=mean(rLCP_res[,1])
  
  return(c(coverage1,coverage2))
}


#--------------------------------------------------------------------
#-----------------------multivariate experiments---------------------
#--------------------------------------------------------------------
#choice of dimensions
dseq=c(1,5*1:4)                
numcores=detectCores()-1
cl=makeCluster(numcores)
registerDoParallel(cl)
hseq=c(1,1.5,2,2.5)

#setting up the experiments
setting=1;ntrain=2000;ncalib=2000;ntest=2000
nrep=50
result1=matrix(0,nrow=5*4,ncol=6)
id=1
for(j in 1:4){
  h=hseq[j]
  for(l in 1:5){
    d=dseq[l]
    print(c(h,d))
    result_h=foreach(k=1:nrep,.combine=rbind,
                     .packages = c("MASS","parallel","doParallel",
                                   "foreach","mvtnorm")) %dopar% {
                                     mrand_onerep(k,d,setting)
                                   }
    result1[id,]=c(colMeans(result_h),apply(result_h,2,sd)/sqrt(nrep),d,h)
    print(result1[id,])
    id=id+1
  }
}
stopCluster(cl)
write.csv(result1,"../results/setting_1_mrand_d.csv")

#---------------------------------------------------
#----------------Visualization----------------------
#---------------------------------------------------
result=read.csv(paste0("../results/setting_1_mrand_d",".csv"))[,-1]
plot_result=matrix(0,nrow=length(hseq)*length(dseq)*11,ncol=5)
id=1
for(i in 1:length(hseq)){
  for(j in 1:length(dseq)){
    plot_result[((id-1)*11+1):(id*11),1]=result1[id,1:11]
    plot_result[((id-1)*11+1):(id*11),2]=result1[id,12:22]
    plot_result[((id-1)*11+1):(id*11),3]=result1[id,23]
    plot_result[((id-1)*11+1):(id*11),4]=result1[id,24]
    plot_result[((id-1)*11+1):(id*11),5]=c(10*1:10,1)
    id=id+1
  }
}
plot_result=as.data.frame(plot_result)
colnames(plot_result)=c("coverage","se","dimension","bandwidth","m")
plot_result$m=as.factor(plot_result$m)

pdf(file="../results/mRLCP_coverage_results.pdf",
    width=8,height=6)
ggplot(plot_result[plot_result$m %in% c(1,10,30,60,100),], aes(x=dimension, y = coverage, col=m)) +
  geom_line() + 
  geom_point()+
  geom_errorbar(aes(ymin=coverage-se, ymax=coverage+se), width=.02,
                position=position_dodge(0.01))+
  geom_hline(yintercept=0.9,linetype="dashed")+
  facet_wrap(.~bandwidth,ncol=2)+
  theme(plot.title = element_text(hjust=0.5))
dev.off()
