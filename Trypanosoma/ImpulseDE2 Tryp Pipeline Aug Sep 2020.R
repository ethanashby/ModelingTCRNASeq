library(tidyverse)
library(ImpulseDE2)
library(Mfuzz)
library(cluster)
library(factoextra)
library(ggfortify)

#read in data and clean... ONLY first 24 hrs
tccount<-read.table("./Trypanosoma/18-5-23_dedup_Proc_ibet_TC_for_deseq_w3hr.txt", header=TRUE)
tccount<- tccount %>% dplyr::distinct(gene, .keep_all=TRUE)
rownames(tccount)<-tccount$gene
tccount<-tccount[,c(2:31)]

#filter genes that are all 0's
tccount<-tccount[rowSums(tccount)>0,]

#########Read count data for first 24hrs
tccount<-tccount[,c(1:3, 5:6, 7:9, 10:12, 13:15, 16:30)]

#avoid problems with log link by adding 1 to each entry
tccount<-tccount+1

#impulse set up
times<-c(rep(0,3), rep(3,2), rep(6,3), rep(12,3),rep(24,3),rep(48,3), rep(72,3), rep(168,3), rep(240,3), rep(336, 3))

impulse_design<-data.frame(Sample=colnames(tccount), Condition=rep("case",29), Time=times)

#edit count matrix and filter s.t. average count is double digits
tccount<-as.matrix(tccount)

#may take awhile to run ImpulseDE2 (~20 mins), using q-value threshold of 0.01
tryp_impulse<-ImpulseDE2::runImpulseDE2(tccount, impulse_design, boolCaseCtrl=FALSE, boolIdentifyTransients=TRUE, vecConfounders=NULL, scaNProc=8, scaQThres = 0.01)

impulse_DEGs<-tryp_impulse$vecDEGenes


#####work on this tomorrow
impulse_df<-matrix(ncol=6)
sigmoid_df<-matrix(ncol=4)
for (i in 1:length(impulse_DEGs)){
  tmp<-tryp_impulse$dfImpulseDE2Results[grep(impulse_DEGs[1], tryp_impulse$dfImpulseDE2Results$Gene),]
  if(tmp$isTransient==TRUE){
    a<-tryp_impulse@lsModelFits$case[[grep(tmp$Gene, names(tryp_impulse@lsModelFits$case))]]$lsImpulseFit$vecImpulseParam
    impulse_df<-rbind(impulse_df, a)
  }
  if(tmp$isTransient==FALSE){
    a<-tryp_impulse@lsModelFits$case[[grep(tmp$Gene, names(tryp_impulse@lsModelFits$case))]]$lsSigmoidFit$vecSigmoidParam
  }
}
