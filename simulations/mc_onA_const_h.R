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
#-----------computing coverage on sets--------------
#---------------------------------------------------
mccs_coverage_onerep=function(k,setting,d,p_seq){
  
  #----------training data and learning score------------
  train_data=simulation(n,d,setting)
  Xnames=paste("X", 1:d, sep="")
  formula=as.formula(paste("Y ~ ", paste(Xnames, collapse= "+")))
  model_lm=lm(formula,data=train_data)
  
  #----------calibration data--------------------------
  calib_data=simulation(n,d,setting)
  Xcalib=as.matrix(calib_data[,-1])
  Vcalib=abs(calib_data$Y-predict(model_lm,newdata=calib_data))
  
  #-------------test data------------------------------
  test_data=simulation(n,d,setting)
  Xtest=as.matrix(test_data[,-1])
  Vtest=abs(test_data$Y-predict(model_lm,newdata=test_data))
  
  #-----------evaluating the competing methods-----------
  RLCP_res=RLCP(Xcalib,Vcalib,Xtest,Vtest,"gaussian",h,alpha)
  calLCP_res=calLCP(Xcalib,Vcalib,Xtest,Vtest,"gaussian",h,alpha)
  baseLCP_res=baseLCP(Xcalib,Vcalib,Xtest,Vtest,"gaussian",h,alpha)
  
  #----------coverage-----------------------
  for(i in 1:length(p_seq)){
    threshold=qchisq(p_seq[i],d)
    ID=ind_set(Xtest,threshold)
    
    coverage1i=mean(RLCP_res[ID,1])
    coverage2i=mean(calLCP_res[ID,1])
    coverage3i=mean(baseLCP_res[ID,1])
    
    coverage1o=mean(RLCP_res[-ID,1])
    coverage2o=mean(calLCP_res[-ID,1])
    coverage3o=mean(baseLCP_res[-ID,1])
    
    coverage1m=mean(RLCP_res[,1])
    coverage2m=mean(calLCP_res[,1])
    coverage3m=mean(baseLCP_res[,1])
    
    coverages=c(coverage1m,coverage1i,coverage1o,
                coverage2m,coverage2i,coverage2o,
                coverage3m,coverage3i,coverage3o)
  }
  return(coverages)
}

#-----------------------------------------------------------
#------------------multivariate experiments-------------------
#-----------------------------------------------------------

#choices of dimensions
dseq=c(1,5*1:10)

numcores=detectCores()-1
cl=makeCluster(numcores)
registerDoParallel(cl)

#bandwidth and choice of p
hseq=2
p_seq=0.5

setting=1;ntrain=2000;ncalib=2000;ntest=2000
nrep=30 

#the inner set
ind_set=function(X,threshold){which(apply(X,1,FUN=function(x) sum(x^2)<=threshold)==1)}


for(i in 1:length(hseq)){
  resultd_1=matrix(0,nrow=length(dseq),ncol=length(p_seq)*18)
  h=hseq[i]
  for(j in 1:length(dseq)){
    print(dseq[j])
    result_h=foreach(k=1:nrep,.combine="rbind",
                     .packages = c("MASS","parallel","doParallel",
                                   "foreach","mvtnorm")) %dopar% {
                                     mccs_coverage_onerep(k,setting,dseq[j],p_seq)
                                   }
    resultd_1[j,]=c(colMeans(result_h),apply(result_h,2,sd)/sqrt(nrep))
  }
  print(resultd_1)
  write.csv(resultd_1,paste0("../results/setting_1_mc_cov_shift_d_",h,".csv"))
}
stopCluster(cl)


#-----------------------------------------------------------------
#------visualizing multivariate experimental results--------------
#-----------------------------------------------------------------

resultd_1=read.csv(paste0("../results/setting_1_mc_cov_shift_d_",h,".csv"))[,-1]
resultd_1=as.data.frame(resultd_1);dseq=c(1,5*1:10)

plot_result1=matrix(0,nrow=9*length(dseq)*length(p_seq),ncol=5)
plot_result1[,1]=stack(resultd_1[,1:9])$values
plot_result1[,2]=rep(dseq,9*length(p_seq))
plot_result1[,3]=rep(rep(c("whole","inner","outer"),each=length(dseq)),3*length(p_seq))
plot_result1[,4]=rep(rep(c("RLCP","calLCP","baseLCP"),each=3*length(dseq)),length(p_seq))
plot_result1[,5]=stack(resultd_1[,10:18])$values

plot_result1=as.data.frame(plot_result1)
colnames(plot_result1)=c("coverage","d","Set","method","se")
plot_result1$coverage=as.numeric(plot_result1$coverage)
plot_result1$se=as.numeric(plot_result1$se)
plot_result1$d=as.numeric(plot_result1$d)
level_order=c('baseLCP','calLCP','RLCP')

pdf(file = "../results/coverage_trend_bandwidth_2.pdf",width = 8,height = 5) 
ggplot(plot_result1, aes(x = d, y = coverage,linetype = Set,color=method)) +
  geom_line() + 
  geom_point()+
  geom_errorbar(aes(ymin=coverage-se, ymax=coverage+se), width=.05)+
  scale_x_continuous(n.breaks=4)+
  geom_hline(yintercept=0.9,linetype="dashed")+
  scale_color_manual(values=c("RLCP"="maroon","calLCP"="gold3","baseLCP"="blue"),guide="none")+
  scale_linetype_manual(labels=c("whole"= expression(B["in"]*union(B["out"])),
                                 "inner"=expression(B["in"]),outer=expression(B["out"])),
                        values=c("whole"="solid","inner"="dotted","outer"="twodash"))+
  facet_grid(.~factor(method,levels=level_order))+
  theme(legend.text.align = 0,
        panel.spacing = unit(0.5,"cm",data=NULL))+
  labs(linetype="Coverage over set",y="Coverage",x="Dimension d")

dev.off()
