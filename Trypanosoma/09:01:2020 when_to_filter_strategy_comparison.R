################
#09/01/2020
#Testing which order we should filter conduct clustering in
#Strategy 1: filter genes with significant expression change (padj<0.01), cluster, filter by membership (>0.7), filter by |lfc|>1
#Strategy 2L filter genes with significant and large expression change (padj<0.01 & |lfc|>1), cluster, filter by membership
################

#read in real RNA-Seq data from Schultz
library(DESeq2)
library(tidyverse)
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
tccount<-tccount[,c(1:3, 5:6, 7:9, 10:12, 13:15)]

#avoid problems with log link by adding 1 to each entry
tccount<-tccount+1

#treat time categorically, use LRT to test for any differences between any times (i.e. expression changes)
design<-data.frame(time=as.factor(c(rep("0", 3), rep("3", 2), rep("6", 3), rep("12", 3), rep("24",3))), rep=as.factor(c(1:3, 1:2, rep(1:3, 3))))

#create DESeqDataSet Object, design takes into account time
ddstcCountTable <- DESeqDataSetFromMatrix(
  countData = tccount,
  colData = design,
  design = ~time)

#run DESEQ with LRT and empty reduced model
trypdds<-DESeq(ddstcCountTable, test="LRT", reduced= ~1)

results3<-results(trypdds, name="time_3_vs_0") %>% as.data.frame()
results3<- results3 %>% mutate(gene=rownames(results3))
results6<-results(trypdds, name="time_6_vs_0") %>% as.data.frame()
results6<- results6 %>% mutate(gene=rownames(results6))
results12<-results(trypdds, name="time_12_vs_0") %>% as.data.frame()
results12<- results12 %>% mutate(gene=rownames(results12))
results24<-results(trypdds, name="time_24_vs_0") %>% as.data.frame()
results24<- results24 %>% mutate(gene=rownames(results24))

results<-rbind(results3, results6, results12, results24)

results$time<-rep(c("3 hr", "6 hr", "12 hr", "24 hr"), each=dim(results3)[1])
results$time<-factor(results$time, levels=c("3 hr", "6 hr", "12 hr", "24 hr"))

#volcano plot
results %>% mutate(thresh=as.factor(ifelse(abs(log2FoldChange)>1 & padj<0.01, 2, ifelse(abs(log2FoldChange)<1 & padj>0.01, 0, 1)))) %>% ggplot(aes(x=log2FoldChange, y=-log10(padj), color=thresh))+geom_point()+scale_color_manual(name="Thresholds", values=c("black", "navy", "dodgerblue"), labels=c("padj>0.01 AND log2FoldChange<1", "padj<0.01 OR log2FoldChange>1", "padj<0.01 AND log2FoldChange>1"))+theme_bw()+geom_hline(yintercept=2, color="red", linetype=2)+geom_vline(xintercept=1, color="red", linetype=2)+geom_vline(xintercept=-1, color="red", linetype=2) +facet_wrap(~time)+theme(strip.text = element_text(size=14))



###############################
#Strategy 1: filter by lfc after clustering
###############################


DEG_df<- results %>% filter(padj<0.01)

DEG_list<-unique(DEG_df$gene)

####Mfuzz DEGS

#DEGcounts
DEGcounts<-counts(trypdds, normalized = TRUE)
DEGcounts<-DEGcounts[rownames(DEGcounts) %in% DEG_list,]

###Too many clusters in all TPs. Let's only look at first 24 hours
cluster_data<-matrix(ncol=5, nrow=dim(DEGcounts)[1])
for (i in 1:dim(DEGcounts)[1]){
  cluster_data[i,]<-c(mean(DEGcounts[i, 1:3]), mean(DEGcounts[i, 4:5]), mean(DEGcounts[i, 6:8]), mean(DEGcounts[i, 9:11]), mean(DEGcounts[i, 12:14]))}
rownames(cluster_data)<-rownames(DEGcounts)
colnames(cluster_data)<-c(0,3,6,12,24)
Set<-ExpressionSet(cluster_data)
Set.s<-standardise(Set)

#Dmin to choose # clusters=4
Dmin(Set.s, m=mestimate(Set.s), crange=seq(2,10,1), repeats=10, visu=TRUE)

cl<-mfuzz(Set.s, c=4, m=mestimate(Set.s))

#plot clusters
#pdf("tryp_fuzzy_clusters.pdf", onefile = TRUE)
#mfuzz.plot2(Set.s, cl, mfrow=c(2,2), time.labels = c(0,3,6,12,24), time.points=c(0,3,6,12,24), min.mem=0)
#dev.off()

#return core membership>0.7 genes
cores<-acore(Set.s, cl, min.acore=0.7)

core_members<-lapply(1:4, function(i){cores[[i]]$NAME}) %>% unlist()

#find members in DEG list and filter with high lfc, grep for cluster members
#note "\\b" denote exact matches, as we don't want genes like "Tb927.7.240" matching to "Tb927.7.2400"

lfc_DEGs<-DEG_df[DEG_df$gene %in% core_members & abs(DEG_df$log2FoldChange)>1,]$gene %>% unique()

cluster_numbers<-lapply(1:length(lfc_DEGs), function(i){grep(paste("\\b", lfc_DEGs[i], "\\b", sep=""), cores)})

strat_1_output<-data.frame("DEGs"=lfc_DEGs, "Cluster"=cluster_numbers %>% unlist())
strat_1_output<-strat_1_output %>% arrange(Cluster)

###############################
#Strategy 2: filter by lfc before clustering
###############################

DEG_df<- results %>% filter(padj<0.01 & abs(log2FoldChange)>1)

DEG_list<-unique(DEG_df$gene)

####Mfuzz DEGS

#DEGcounts
DEGcounts<-counts(trypdds, normalized = TRUE)
DEGcounts<-DEGcounts[rownames(DEGcounts) %in% DEG_list,]

###Too many clusters in all TPs. Let's only look at first 24 hours
cluster_data<-matrix(ncol=5, nrow=dim(DEGcounts)[1])
for (i in 1:dim(DEGcounts)[1]){
  cluster_data[i,]<-c(mean(DEGcounts[i, 1:3]), mean(DEGcounts[i, 4:5]), mean(DEGcounts[i, 6:8]), mean(DEGcounts[i, 9:11]), mean(DEGcounts[i, 12:14]))}
rownames(cluster_data)<-rownames(DEGcounts)
colnames(cluster_data)<-c(0,3,6,12,24)
Set<-ExpressionSet(cluster_data)
Set.s<-standardise(Set)

#Dmin to choose # clusters=3, but for consistency we might also want 4
Dmin(Set.s, m=mestimate(Set.s), crange=seq(2,10,1), repeats=10, visu=TRUE)



############
#N_clusters=3
cl<-mfuzz(Set.s, c=3, m=mestimate(Set.s))

#plot clusters
#pdf("tryp_fuzzy_clusters.pdf", onefile = TRUE)
#mfuzz.plot2(Set.s, cl, mfrow=c(2,2), time.labels = c(0,3,6,12,24), time.points=c(0,3,6,12,24), min.mem=0)
#dev.off()

#return core membership>0.7 genes
cores<-acore(Set.s, cl, min.acore=0.7)

#extract core genes 
core_members<-lapply(1:3, function(i){cores[[i]]$NAME}) %>% unlist()

strat_2_output_clust3<-data.frame("DEGs"=core_members, "Cluster"=c(rep(1, dim(cores[[1]])[1]), rep(2, dim(cores[[2]])[1]), rep(3, dim(cores[[3]])[1])))

###########
#N_clusters=4
cl<-mfuzz(Set.s, c=4, m=mestimate(Set.s))

#plot clusters
#pdf("tryp_fuzzy_clusters.pdf", onefile = TRUE)
#mfuzz.plot2(Set.s, cl, mfrow=c(2,2), time.labels = c(0,3,6,12,24), time.points=c(0,3,6,12,24), min.mem=0)
#dev.off()

#return core membership>0.7 genes
cores<-acore(Set.s, cl, min.acore=0.7)

#extract core genes 
core_members<-lapply(1:4, function(i){cores[[i]]$NAME}) %>% unlist()

strat_2_output_clust4<-data.frame("DEGs"=core_members, "Cluster"=c(rep(1, dim(cores[[1]])[1]), rep(2, dim(cores[[2]])[1]), rep(3, dim(cores[[3]])[1]), rep(4, dim(cores[[4]])[1])))

###### write files out as csv
write.csv(strat_1_output, "Strat_1_Output.csv")
write.csv(strat_2_output_clust3, "Strat_2_Output_nclust3.csv")
write.csv(strat_2_output_clust4, "Strat_2_Output_nclust4.csv")


#####################
#How do these methods compare

#overall overlap
area1<-length(strat_1_output$DEGs)
area2<-length(strat_2_output_clust3$DEGs)
area3<-length(strat_2_output_clust4$DEGs)
n12<-length(intersect(strat_1_output$DEGs, strat_2_output_clust3$DEGs))
n23<-length(intersect(strat_2_output_clust3$DEGs, strat_2_output_clust4$DEGs))
n13<-length(intersect(strat_1_output$DEGs, strat_2_output_clust4$DEGs))
n123<-length(intersect(intersect(strat_1_output$DEGs, strat_2_output_clust3$DEGs), strat_2_output_clust4$DEGs))

draw.triple.venn(area1, area2, area3, n12, n23, n13, n123, category = c("Strat_1 C=4", "Strat_2 C=3", "Strat_2 C=4"), fill = c("light blue", "light green", "pink"), scaled=TRUE, cat.pos=2)


