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

DEG_df<- results %>% filter(padj<0.01 & abs(log2FoldChange)>1)

DEG_list<-unique(DEG_df$gene)

DEG_df

#####this is final DEG list containing +/- 2 fold change and padj<0.01 for all genes
DEG_list

#nice venn diagram
set1<-results3 %>% filter(padj<0.01 & abs(log2FoldChange)>1) %>% select(gene) %>% unlist()
set2<- results6 %>% filter(padj<0.01 & abs(log2FoldChange)>1) %>% select(gene) %>% unlist()
set3<- results12 %>% filter(padj<0.01 & abs(log2FoldChange)>1) %>% select(gene) %>% unlist()
set4<- results24 %>% filter(padj<0.01 & abs(log2FoldChange)>1) %>% select(gene) %>% unlist()

area1<-length(set1)
area2<-length(set2)
area3<-length(set3)
area4<-length(set4)
n12 <- length(intersect(set1, set2))
n13<- length(intersect(set1, set3))
n14 <- length(intersect(set1, set4))
n23 <- length(intersect(set2, set3))
n24 <- length(intersect(set2, set4))
n34 <- length(intersect(set3, set4))
n123 <- length(intersect(intersect(set1, set2), set3))
n124 <- length(intersect(intersect(set1, set2), set4))
n134 <- length(intersect(intersect(set1, set3), set4))
n234 <- length(intersect(intersect(set2, set3), set4))
n1234 <- length(intersect(intersect(intersect(set1, set2), set3), set4))

draw.quad.venn(area1, area2, area3, area4, n12, n13, n14, n23, n24, n34, n123, n124, n134, n234, n1234, category = c("3h","6h", "12h","24h"), fill = c("light blue", "light green", "pink", "cyan"))

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

#Dmin to choose # clusters=7
Dmin(Set.s, m=mestimate(Set.s), crange=seq(2,15,1), repeats=10, visu=TRUE)

cl<-mfuzz(Set.s, c=5, m=mestimate(Set.s))

#plot clusters
pdf("tryp_fuzzy_clusters.pdf", onefile = TRUE)
mfuzz.plot2(Set.s, cl, mfrow=c(2,3), time.labels = c(0,3,6,12,24), time.points=c(0,3,6,12,24), min.mem=0)
dev.off()

cores<-acore(Set.s, cl, min.acore=0.8)
