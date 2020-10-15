#####################
#DE DESEQ2 analysis script for E.coli work
#DEGs between WT and RPOS-knockout at every TP except 0 min (30,60,90,120,150)
#by Ethan Ashby (PO Class 2021)
######################

library(DESeq2)
library(tidyverse)
#for Venn diagram
library(systemPipeR) #https://www.bioconductor.org/packages/release/bioc/html/systemPipeR.html

##### Madison's parsing code
allCounts <- read.csv("LB_Time_Course_GCA_000005845.2_ASM584v2_genomic_counts copy.tsv", header = T, sep = "\t")
allCounts$GeneidBackup = allCounts$Geneid
allCounts <- allCounts %>% tidyr::separate(GeneidBackup, c("feature", "rest"), sep="[:]")
allCounts %>% group_by(feature) %>% summarise(number_of_genes = n())
# Now, we must extract the genenames from each Geneid. However, each feature has a slightly different pattern, and we will need to do the features separately. 

# IGR's (this includes AS_IGRSs): 
# IGR stands for intergenic region which means a region between coding sequences or different types of RNA. Therefore, we'll have a start.bnum and end.bnum as well as a start.genename and end.genename. Note that not all genes will have a bnum (only CDS/AS_CDS do)
bnum = "b[0-9]{4}" # what do bnumbers look like?
genename = ",[a-z]{3}[A-Z,]." # what does a genename look like? this is regexp lingo
rna.name = ",rna[0-9].." # what does an RNA name look like?
igr <- allCounts %>% filter(feature %in% c("IGR", "AS_IGR"))
igr$GeneidBackup = igr$Geneid # store the Geneid
igr <- igr %>% separate(GeneidBackup, c("Geneid1", "Geneid2"), sep = "[/]") # separate the first part of the Geneid which talks about the IGR's start gene (Geneid1) and the last part of the IGR Geneid which talks about that IGR's end gene (Geneid2). 
igr$feature1 <- separate(igr, Geneid1, c("feature1", "rest"), sep = "[,]")$feature1
igr$feature1 <- separate(igr, feature1, c("rest", "feature1"), sep = "[()]")$feature1 #start feature
igr$feature2 <- separate(igr, Geneid2, c("feature2", "rest"), sep = "[,]")$feature2
igr$start.gene <- case_when( #start gene name: many possibilities!
  igr$feature1 == "CDS" ~ str_extract(igr$Geneid1, genename), #if the start feature was CDS, then the name is going to be genename style, so we extract a genename-type thing from Geneid1
  TRUE ~ str_extract(igr$Geneid1, rna.name)) #otherwise, it's going to have an RNA-style name, so we extract the rna.name from Geneid1
igr$end.gene <- case_when( #end gene neame: similar to above!
  igr$feature2 == "CDS" ~ str_extract(igr$Geneid2, genename), # if the end feature was CDS, then we're looking for a genename-type bit from Geneid2
  TRUE ~ str_extract(igr$Geneid2, rna.name)) #otherwise, it must be an RNA-style label of some sort. 
igr$start.bnum <- case_when(
  igr$feature1 == "CDS" ~ str_extract(igr$Geneid1, bnum), #bnums only exist for CDS, so we check if the feature is CDS before extracting a bnum from Geneid1
  TRUE ~ "none") # if not CDS, then no bnum exists so we can put "none"
igr$end.bnum <- case_when(
  igr$feature2 == "CDS" ~ str_extract(igr$Geneid2, bnum), #same thing as above but for end bnum
  TRUE ~ "none")
# now get rid of all those pesky commas that got into our start.gene labels. I could have not included the punctuation in my regex pattern, but then str_extract() might have gotten confused with a less specific pattern
igr <- igr %>% separate(start.gene, into = c("comma", "start.gene"), sep = "[,]") %>% dplyr::select(-comma) %>% separate(end.gene, into = c("comma", "end.gene"), sep = "[,]") %>% dplyr::select(-comma)
allCounts <- full_join(igr, allCounts) #add this new information to allCounts!
# CDS
# have bnum and genename columns
# left join to allCounts
genename = ":[a-z]{3}.." #new genename pattern
#bnum pattern stays the same
cds <- allCounts %>% filter(feature %in% c("AS_CDS", "CDS")) 
cds$genename <- str_extract(cds$Geneid, genename) #extract those genenames!
cds$bnum <- str_extract(cds$Geneid, bnum) # extract them bnums!
#get rid of the pesky colon that was part of the pattern
cds <- cds %>% separate(genename, into = c("colon", "genename"), sep = ":") %>%
  dplyr::select(-colon)
allCounts <- full_join(allCounts, cds) #add the new info to allCounts
#ncRNA
#ncRNA doesn't have bnums, but id's which we'll put in the genename column
rna.name = ":rna[0-9].." #new rna.name pattern
rna <- allCounts %>% filter(feature %in% c("ncRNA", "AS_ncRNA"))
rna$genename <- str_extract(rna$Geneid, rna.name) #record those rna.names
rna <- rna %>% separate(genename, into = c("colon", "genename"), sep = ":") %>%
  dplyr::select(-colon) #get rid of colon
allCounts <- full_join(allCounts, rna) #update allCounts with ncRNA names!
#rRNA
rRNA <- allCounts %>% filter(feature %in% c("rRNA", "AS_rRNA"))
rRNA$genename <- str_extract(rRNA$Geneid, rna.name) #same rna.name pattern exists as above, so extract those rna.names! And store the result as genename - I know. This is just for convenience's sake so we have a common column to refer to when we want the short hand name for a gene. 
rRNA <- rRNA %>% separate(genename, into = c("colon", "genename"), sep = ":") %>%
  dplyr::select(-colon) #get rid of that colon
allCounts <- full_join(allCounts, rRNA) #update allCounts
#tRNA -- analogous to rRNA above. 
tRNA <- allCounts %>% filter(feature %in% c("tRNA", "AS_tRNA"))
tRNA$genename <- str_extract(tRNA$Geneid, rna.name)
tRNA <- tRNA %>% separate(genename, into = c("colon", "genename"), sep = ":") %>%
  dplyr::select(-colon)
allCounts <- full_join(tRNA, allCounts)
# remove the NA rows we just created by full_joining while adding the ncRNA, rRNA, tRNA genenames
allCounts <- filter(allCounts, feature %in% c("IGR", "AS_IGR") | genename != "NA")

#filter by CDS, NA, 0 sums, resolve duplicates
#filter
ecounts<-allCounts %>% filter(feature=="CDS")
#remove NA and rows of 0's
ecounts<-ecounts[rowSums(is.na(ecounts)) != ncol(ecounts), ]
ecounts<-subset(ecounts, rowSums(ecounts[,2:37])>0)
#resolve duplicates (takes the first row corresponding to a genename)
ecounts<- ecounts %>% dplyr::distinct(genename, .keep_all=TRUE)
rownames(ecounts)<-ecounts$genename
ecounts<-ecounts[,2:37]




######## Differential expression

#avoid problems with log link by adding 1 to each entry
ecounts<-ecounts+1

#treat time categorically, use LRT to test for any differences between any times (i.e. expression changes)
colidesign<-data.frame(time=as.factor(rep(c("0", "120", "60", "90", "120", "150"), 6)), rep=as.factor(c(rep(1, 12), rep(2,12), rep(3,12))), treatment=as.factor(c(rep("WT",6), rep("treat",6), rep("WT",6), rep("treat",6), rep("WT",6), rep("treat",6)))) %>% mutate(timetreat=paste(treatment, time, sep=""))

#print design matrix
colidesign

#Generate DESEQ dataset with design sourced from here: https://support.bioconductor.org/p/62684/
# According to Michael Love DESEQ2 author: "The design [below] controls for differences between horse at time=0, and then tests to see if there are any differences at any time point."
ddsecoliCountTable <- DESeqDataSetFromMatrix(
  countData = ecounts,
  colData = colidesign,
  design = ~ time + treatment + time:treatment)

#ID genes DE over time (30 min ctrl vs 30 min WT, 60 min ctrl vs 60 min WT, ...)
ecolidds<-DESeq(ddsecoliCountTable, test="LRT", reduced= ~ time + treatment)

#gives us contrasts to look at
resultsNames(ecolidds)

results30<-DESeq2::results(ecolidds, name="time30.treatmentWT") %>% as.data.frame()
results30<- results30 %>% mutate(gene=rownames(results30))
results60<-DESeq2::results(ecolidds, name="time60.treatmentWT") %>% as.data.frame()
results60<- results60 %>% mutate(gene=rownames(results60))
results90<-DESeq2::results(ecolidds, name="time90.treatmentWT") %>% as.data.frame()
results90<- results90 %>% mutate(gene=rownames(results90))
results120<-DESeq2::results(ecolidds, name="time120.treatmentWT") %>% as.data.frame()
results120<- results120 %>% mutate(gene=rownames(results120))
results150<-DESeq2::results(ecolidds, name="time150.treatmentWT") %>% as.data.frame()
results150<- results150 %>% mutate(gene=rownames(results150))

results<-rbind(results30, results60, results90, results120, results150)

results$time<-rep(c("30", "60", "90", "120", "150"), each=dim(results30)[1])
results$time<-factor(results$time, levels=c("30", "60", "90", "120", "150"))
results<-results[complete.cases(results),]

#volcano plot
results %>% mutate(thresh=as.factor(ifelse(abs(log2FoldChange)>1 & padj<0.01, 2, ifelse(abs(log2FoldChange)<1 & padj>0.01, 0, 1)))) %>% ggplot(aes(x=log2FoldChange, y=-log10(padj), color=thresh))+geom_point()+scale_color_manual(name="Thresholds", values=c("grey", "black", "darkorange1"), labels=c("padj>0.01\nAND log2FoldChange<1", "padj<0.01\nOR log2FoldChange>1", "padj<0.01\nAND log2FoldChange>1"))+theme_bw()+geom_hline(yintercept=2, color="blue", linetype=2)+geom_vline(xintercept=1, color="blue", linetype=2)+geom_vline(xintercept=-1, color="blue", linetype=2) +facet_grid(~time)+theme(strip.text = element_text(size=14), legend.text=element_text(size=12), legend.title=element_text(size=12), legend.position="bottom")

#just pull out the significant results, should get 745 unique DEGs
sig.results<-results[results$padj<0.01,]

#turns out they are all diff expressed at every TP
vennset <- overLapper(split(sig.results[,7], sig.results$time)[1:5], type="vennsets")
vennPlot(vennset, mymain="DEG intersections")

#crappy heatmap
sig.results[order(sig.results$padj),] %>% head(100*5) %>% ggplot(aes(x=time, y=gene, fill=log2FoldChange))+geom_tile()+scale_fill_gradient2(low="blue", mid="white", high="orange", midpoint = 5)
