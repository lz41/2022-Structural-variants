library(ggiraphExtra)
library(ggiraph)
#####################################################
#### get genes coordinates #####
gene_annotation<-read.gff("melissa_functional_blastnew.gff", na.strings = c(".", "?"), GFF3 = TRUE)
gene_annotation$scaffold<-paste(paste(gene_annotation$seqid,gene_annotation$start,sep=":"),gene_annotation$end,sep = "-")

genelist<-read.table("melissa_genelistnew.txt")
head(genelist)
length(genelist[,1])
colnames(genelist)<-c("scaffold","start","end","annotation","id")
genelist<-genelist[,-4]
genelist$annotation<-gene_annotation$attributes

genecoord<-read.table("RevisedAnnotationCoords.txt",header=TRUE)
head(genecoord)
length(genecoord[,1])
genelistnew<-merge(genecoord,genelist[,c("id","annotation")],by="id")
head(genelistnew)

for(i in 1:length(genelistnew[,1])){
  genelistnew$scaffold[i]<-strsplit(genelistnew$id[i],":")[[1]][1]
}
colnames(genelistnew)[c(2,3)]<-c("start","end")
##### write.table(genelistnew,"genelistnew.txt")
head(genelistnew)
####################################################
############# don't run #########################
write.table(bg_filter4.sum,"bg_filter4.sum.txt")
bg_filter4.sum<-read.table("bg_filter4.sum.txt",header=TRUE)

##################################################
##### get the number of genes for each SV ##########
scaffold<-unique(bg_filter4.sum$scaffold)
scaffold_annotation<-read.table("genelistnew.txt")
head(scaffold_annotation$id)
bg_filter4.sum$scaffold<-paste("Scaffold",bg_filter4.sum$scaffold,sep="_")
### write.table(scaffold,"scaffold.txt")
overlap.genes<-list()
overlap.alleles<-list()
for(i in 1:length(bg_filter4.sum[,1])){
  annotation<-scaffold_annotation[scaffold_annotation$scaffold==bg_filter4.sum$scaffold[i],]
  overlap.alleles<-annotation[(annotation$start>bg_filter4.sum[i,"start"] & annotation$start<bg_filter4.sum[i,"end"])|(annotation$end>bg_filter4.sum[i,"start"] & annotation$end<bg_filter4.sum[i,"end"]),]
  bg_filter4.sum$n.genes[i]<-length(overlap.alleles[,1])
}
hist(bg_filter4.sum$n.genes)
summary(bg_filter4.sum$n.genes)
############################################################
#### get the list of gene function for deletions ##########
overlap.genes<-list()
overlap.alleles<-list()
overlap.genes<-data.frame(scaffold_annotation[1,])
bg_filter4.DEL<-bg_filter4.sum[bg_filter4.sum$type=="DEL",]
for(i in 1:length(bg_filter4.DEL[,1])){
 annotation<-scaffold_annotation[scaffold_annotation$scaffold==bg_filter4.DEL$scaffold[i],]
 overlap.alleles<-annotation[(annotation$start>bg_filter4.DEL[i,"start"] & annotation$start<bg_filter4.DEL[i,"end"])|
 (annotation$end>bg_filter4.DEL[i,"start"] & annotation$end<bg_filter4.DEL[i,"end"]),]
 overlap.genes<-rbind(overlap.genes,overlap.alleles)
}
overlap.genes<-overlap.genes[-1,]
head(overlap.genes)

for(i in 1:length(overlap.genes[,1])){
  overlap.genes$gene.list[i]<-strsplit(overlap.genes$annotation[i],";Note=")[[1]][2]
  overlap.genes$gene.name[i]<-strsplit(strsplit(overlap.genes$annotation[i],";Alias=")[[1]][2],";Note=")[[1]][1]
}
head(overlap.genes)
overlap.genesnew<-overlap.genes[overlap.genes$gene.list!="Protein of unknown function;",]
length(overlap.genesnew[,1])
### write.csv(overlap.genesnew[,c("gene.name","gene.list")],"deleted.gene.list.csv")
############# check gene functions where SVs exhibit highest degree of Fst #######
overlap.genes<-list()
overlap.alleles<-list()
overlap.genes<-data.frame(scaffold_annotation[1,])
bg_filter4.high<-bg_filter4.sum[order(bg_filter4.sum$Fst,decreasing = TRUE),]
bg_filter4.high<-bg_filter4.high[1:4,]

for(i in 1:length(bg_filter4.high[,1])){
  annotation<-scaffold_annotation[scaffold_annotation$scaffold==bg_filter4.high$scaffold[i],]
  overlap.alleles<-annotation[(annotation$start>bg_filter4.high[i,"start"] & annotation$start<bg_filter4.high[i,"end"])|
                                (annotation$end>bg_filter4.high[i,"start"] & annotation$end<bg_filter4.high[i,"end"]),]
  overlap.genes<-rbind(overlap.genes,overlap.alleles)
}
overlap.genes<-overlap.genes[-1,]
head(overlap.genes)

for(i in 1:length(overlap.genes[,1])){
  overlap.genes$gene.list[i]<-strsplit(overlap.genes$annotation[i],";Note=")[[1]][2]
  overlap.genes$gene.name[i]<-strsplit(strsplit(overlap.genes$annotation[i],";Alias=")[[1]][2],";Note=")[[1]][1]
}

overlap.genes<-list()
overlap.alleles<-list()
overlap.genes<-data.frame(scaffold_annotation[1,])
bg_filter4.highalpha<-bg_filter4.sum[order(bg_filter4.sum$alpha.mean),]
bg_filter4.highalpha<-bg_filter4.highalpha[1:71,]
bg_filter4.highalpha<-bg_filter4.highalpha[bg_filter4.highalpha$type=="DEL",]


for(i in 1:length(bg_filter4.highalpha[,1])){
  annotation<-scaffold_annotation[scaffold_annotation$scaffold==bg_filter4.highalpha$scaffold[i],]
  overlap.alleles<-annotation[(annotation$start>bg_filter4.highalpha[i,"start"] & annotation$start<bg_filter4.highalpha[i,"end"])|
                                (annotation$end>bg_filter4.highalpha[i,"start"] & annotation$end<bg_filter4.highalpha[i,"end"]),]
  overlap.genes<-rbind(overlap.genes,overlap.alleles)
}
overlap.genes<-overlap.genes[-1,]
head(overlap.genes)

for(i in 1:length(overlap.genes[,1])){
  overlap.genes$gene.list[i]<-strsplit(overlap.genes$annotation[i],";Note=")[[1]][2]
  overlap.genes$gene.name[i]<-strsplit(strsplit(overlap.genes$annotation[i],";Alias=")[[1]][2],";Note=")[[1]][1]
}
overlap.genesnew<-overlap.genes[overlap.genes$gene.list!="Protein of unknown function;",]
overlap.genesnew$annotation
######################################################################################
### Figure 6A boxplot of alpha across different SV types ####

par(mfrow=c(2,2),mar=c(7,6,5,5))

boxplot(alpha.mean~type,data=bg_filter4.sum,outline=FALSE,
        xlab="", ylab="", cex.axis=2.2)
mtext(expression(alpha), side=2, line=3, cex=4)
mtext("Types of structural variants", side=1, line=5, cex=2.2)
title("(A)", adj = 0, line =1,cex.main=4)

#### compare alpha of DEL with INV ####
bg_filter4.sum$count<-1
ddply(bg_filter4.sum,"type",summarize,N=sum(count))
sv4.type<-lm(abs(alpha.mean)~type,data=bg_filter4.sum)
emmeans(sv4.type,list(pairwise~type),adjust="tukey")
bg_filter4.del<-bg_filter4.sum[bg_filter4.sum$type=="DEL",]
del_inv.diff<-c()
for (i in 1:10000){
  del.temp<-bg_filter4.del[sample(1:1255,124),]
  del_inv.diff[i]<-mean(del.temp$alpha.mean)-mean(bg_filter4.sum[bg_filter4.sum$type=="INV","alpha.mean"])
}
hist(del_inv.diff)
mean(del_inv.diff)
length(del_inv.diff[del_inv.diff>0])

#######################################################################################
####### correlation between alpha with SV length and number of genes #######
##### standardize length and n. genes ######
bg_filter4.sum$length.sd<-scale(bg_filter4.sum$length)
bg_filter4.sum$n.genes.sd<-scale(bg_filter4.sum$n.genes)
str(bg_filter4.sum)
##### correlation with alpha as response variable ###
inv.m1<-lm(alpha.mean~length.sd*n.genes.sd,data=bg_filter4.sum[bg_filter4.sum$type=="INV",])
summary(inv.m1)
inv.result<-data.frame(summary(lm(data=bg_filter4.sum[bg_filter4.sum$type=="INV",],alpha.mean~length.sd*n.genes.sd))$coefficients)
inv.result

del.result<-data.frame(summary(lm(data=bg_filter4.sum[bg_filter4.sum$type=="DEL",],alpha.mean~length.sd*n.genes.sd))$coefficients)
del.result

ins.result<-data.frame(summary(lm(data=bg_filter4.sum[bg_filter4.sum$type=="INS",],alpha.mean~length.sd*n.genes.sd))$coefficients)
ins.result

dup.result<-data.frame(summary(lm(data=bg_filter4.sum[bg_filter4.sum$type=="DUP",],alpha.mean~length.sd*n.genes.sd))$coefficients)
dup.result

Ngenes.result<-list(del.result[-1,],dup.result[-1,],ins.result[-1,],inv.result[-1,])
names(Ngenes.result)<-c("DEL","DUP","INS","INV")
## write.csv(Ngenes.result,"Ngenes.result.csv")
#### correlation with absolute alpha as response variable ###

inv.result2<-data.frame(summary(lm(data=bg_filter4.sum[bg_filter4.sum$type=="INV",],abs(alpha.mean)~length.sd*n.genes.sd))$coefficients)
inv.result2

del.result2<-data.frame(summary(lm(data=bg_filter4.sum[bg_filter4.sum$type=="DEL",],abs(alpha.mean)~length.sd*n.genes.sd))$coefficients)
del.result2

ins.result2<-data.frame(summary(lm(data=bg_filter4.sum[bg_filter4.sum$type=="INS",],abs(alpha.mean)~length.sd*n.genes.sd))$coefficients)
ins.result2

dup.result2<-data.frame(summary(lm(data=bg_filter4.sum[bg_filter4.sum$type=="DUP",],abs(alpha.mean)~length.sd*n.genes.sd))$coefficients)
dup.result2

Ngenes.result2<-list(del.result2[-1,],dup.result2[-1,],ins.result2[-1,],inv.result2[-1,])
names(Ngenes.result2)<-c("DEL","DUP","INS","INV")
## write.csv(Ngenes.result2,"Ngenes.result2.csv")
########### plot 6B alpha with residuals of n. genes over length ############
m1<-lm(n.genes~length,data=bg_filter4.sum[bg_filter4.sum$type=="INV",])

m2<-lm(bg_filter4.sum$alpha.mean[bg_filter4.sum$type=="INV"]~resid(m1))
summary(m2)
cor.test(bg_filter4.sum$alpha.mean[bg_filter4.sum$type=="INV"],resid(m1))
inv.m1<-lm(alpha.mean~n.genes*length,data=bg_filter4.sum[bg_filter4.sum$type=="INV",])
summary(inv.m1)

plot(resid(m1),bg_filter4.sum$alpha.mean[bg_filter4.sum$type=="INV"],pch=19,
     ylab="",xlab="",ylim=c(-1.5,1.5),cex.axis=2.2)
mtext(expression(alpha), side=2, line=3, cex=4)
mtext("Residuals of the number of genes against length", side=1, line=5, cex=2.2)
abline(m2)
title("(B)", adj = 0, line =1,cex.main=4)
################################################################
##### number of genes in all SVs #####
sv.LKnew2<-read.table("sv.LKnew2.txt",header=TRUE)
SV2<-sv.LKnew2[,c(1:3,121:124)]

for (i in 1:length(SV2[,1])){
  SV2$scaffold[i]<-paste(strsplit(SV2$position[i], "_")[[1]][1],strsplit(SV2$position[i], "_")[[1]][2],sep="_")
  SV2$start[i]<-strsplit(SV2$position[i], "_")[[1]][5]
}
SV2$start<-as.numeric(SV2$start)
SV2$end<-SV2$start+SV2$length

overlap.genes<-list()
overlap.alleles<-list()
for(i in 1:length(SV2[,1])){
  annotation<-scaffold_annotation[scaffold_annotation$scaffold==SV2$scaffold[i],]
  overlap.alleles<-annotation[(annotation$start>SV2[i,"start"] & annotation$start<SV2[i,"end"])|(annotation$end>SV2[i,"start"] & annotation$end<SV2[i,"end"]),]
  SV2$n.genes[i]<-length(overlap.alleles[,1])
}
head(SV2)
cor.test(SV2$length,SV2$n.genes)
### write.table(SV2,"SV2new.txt")
cor.test(bg_filter4.sum$alpha.mean,bg_filter4.sum$Fst)
cor.test(abs(bg_filter4.sum$alpha.mean),bg_filter4.sum$Fst)
cor.m1<-lm(abs(alpha.mean)~Fst,data=bg_filter4.sum)
######## plot 6C correlation between Fst and abs alpha ###
plot(bg_filter4.sum$Fst,abs(bg_filter4.sum$alpha.mean),pch=19,
     ylab="",xlab="",cex.axis=2.2)
mtext(expression(italic(alpha)), side=2, line=3, cex=4)
mtext(expression(italic("Fst")), side=1, line=5, cex=2.2)
text(0.68,1.4,expression(paste(italic("R"),"= 0.171, ",italic("P"),"< 0.001")),cex=3)
abline(cor.m1)
title("(C)", adj = 0, line =1,cex.main=4)

SV2<-read.table("SV2new.txt")
SV.LK2<-read.table("sv.LKfilter2.txt")
head(SV2)
head(SV.LK2)
SV2$Fst<-SV.LK2$Fst
SV2$length.sd<-scale(SV2$length)
SV2$n.genes.sd<-scale(SV2$n.genes)
SV2.new<-SV2[!is.na(SV2$Fst),]
SV2.new$Fst.logit<-log(SV2.new$Fst/(1-SV2.new$Fst)+0.0001)
Fstm.result<-list(summary(lm(Fst.logit~length.sd*n.genes.sd,data=SV2.new[SV2.new$type=="DEL",]))$coefficient,
summary(lm(Fst.logit~length.sd*n.genes.sd,data=SV2.new[SV2.new$type=="DUP",]))$coefficient,
summary(lm(Fst.logit~length.sd*n.genes.sd,data=SV2.new[SV2.new$type=="INS",]))$coefficient,
summary(lm(Fst.logit~length.sd*n.genes.sd,data=SV2.new[SV2.new$type=="INV",]))$coefficient)

names(Fstm.result)<-c("DEL","DUP","INS","INV")
Fstm.result

##### write.csv(Fstm.result,"Fstm.result.csv")

SV2$length.new<-log(SV2$length)

summary(lm(Fst~length.new*n.genes,data=SV2))
#### 