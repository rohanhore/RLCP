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
dseq=c(1,5,10*1:5)
n=2000;setting=4
#---set prefixed effective sample size----
eff_size=50

#---finding optimum bandwidth choices of RLCP and calLCP---
optimum_RLCP_bandwidth=optimum_calLCP_bandwidth=rep(0,length(dseq))
for(i in 1:length(dseq)){
  if(i==1){h_min_RLCP=h_min_calLCP=0.02}
  else{h_min_RLCP=optimum_RLCP_bandwidth[i-1];h_min_calLCP=optimum_calLCP_bandwidth[i-1]}
  Xtrain=as.matrix(simulation(n,dseq[i],setting)[,-1]) 
  optimum_RLCP_bandwidth[i]=optimum_RLCP_h(Xtrain,"gaussian",h_min_RLCP,eff_size)
  optimum_calLCP_bandwidth[i]=optimum_calLCP_h(Xtrain,"gaussian",h_min_calLCP,eff_size)
}

#-----------------------------------------------------------
#-----------computing marginal coverage on sets-------------
#-----------------------------------------------------------
mccs_coverage_optimized=function(k,setting,j,p_seq){
  d=dseq[j]
  
  #----------training data and learning score------------
  train_data=simulation(n,d,setting)
  Xnames=paste("X", 1:d, sep="")
  formula=as.formula(paste("Y ~ ", paste(Xnames, collapse= "+")))
  model_lm=lm(formula,data=train_data)
  
  #---------consider the optimum bandwidth--------------
  Xtrain=as.matrix(train_data[,-1])
  h_opt_RLCP=optimum_RLCP_bandwidth[j]
  h_opt_calLCP=optimum_calLCP_bandwidth[j]
  
  #----------calibration data--------------------------
  calib_data=simulation(n,d,setting)
  Xcalib=as.matrix(calib_data[,-1])
  scores_calib=abs(calib_data$Y-predict(model_lm,newdata=calib_data))
  
  #-------------test data------------------------------
  test_data=simulation(n,d,setting)
  Xtest=as.matrix(test_data[,-1])
  scores_test=abs(test_data$Y-predict(model_lm,newdata=test_data))
  
  #-----------evaluating the competing methods-----------
  RLCP_res=RLCP(Xcalib,scores_calib,Xtest,scores_test,"gaussian",h_opt_RLCP,alpha)
  calLCP_res=calLCP(Xcalib,scores_calib,Xtest,scores_test,"gaussian",h_opt_calLCP,alpha)
  baseLCP_res=baseLCP(Xcalib,scores_calib,Xtest,scores_test,"gaussian",h_opt_calLCP,alpha)
  
  #----------coverage-----------------------
  for(i in 1:length(p_seq)){
    # #for setting 1
    # threshold=qchisq(p_seq[i],d)
    
    #for setting 3/4
    threshold=(3*(1/2)^(1/d))^2
      
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


#------------------------------------------------------------
#-----------------------experimental results-----------------
#------------------------------------------------------------

numcores=detectCores()-1
cl=makeCluster(numcores)
registerDoParallel(cl)

#choice of p
p_seq=0.5

setting=4
nrep=50;alpha=0.1

#for setting 1/4--L_2 balls
ind_set=function(X,threshold){which(apply(X,1,FUN=function(x) sum(x^2)<=threshold)==1)}

# #for setting 3---L_infty balls
# ind_set=function(X,radii){which(apply(X,1,FUN=function(x) max(abs(x))<=radii)==1)}

resultd_1=matrix(0,nrow=length(dseq),ncol=length(p_seq)*18)
for(j in 1:length(dseq)){
  print(dseq[j])
  result_h=foreach(k=1:nrep,.combine="rbind",
                   .packages = c("MASS","parallel","doParallel",
                                 "foreach","mvtnorm")) %dopar% {
                                   mccs_coverage_optimized(k,setting,j,p_seq)
                                 }
  resultd_1[j,]=c(colMeans(result_h[,1:9]),apply(result_h[,1:9],2,sd)/sqrt(nrep))
  print(resultd_1[j,])
}
stopCluster(cl)
print(resultd_1)

#sett1
write.csv(resultd_1,"../results/setting_1_mc_sett1_d_size50_normal.csv")
#------------------------------------------------
#-----------------Visualization------------------
#------------------------------------------------

resultd_1=read.csv("../results/setting_1_mc_sett1_d_size50_normal.csv")[,-1]
resultd_1=as.data.frame(resultd_1)

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

pdf(file = "../results/figures/coverage_trend_effective_size_50_optimized_bandwidth_normal.pdf",width = 8,height=4)

ggplot(plot_result1[plot_result1$method %in% c("calLCP","RLCP"),], aes(x = d, y = coverage,linetype = Set,color=method)) +
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
  labs(linetype="Coverage over set",y="Coverage",x="Dimension d")+
  theme(legend.text.align = 0,
        panel.spacing = unit(0.5,"cm",data=NULL),
        strip.text = element_text(size = 13),
        legend.text=element_text(size=13),
        axis.text=element_text(size=13),
        axis.title=element_text(size=13),
        legend.title=element_text(size=13),
        legend.spacing.y = unit(0.3,'cm'))+
  guides(SET = guide_legend(byrow = TRUE))

dev.off()

