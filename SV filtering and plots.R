###################################################################
library(patchwork)
library(ggplot2)
library(plyr)
library(cowplot)
library(lemon)
#################################################################
sv.DV<-read.table("sv.DV.txt",header=TRUE)
head(sv.DV)
head(sv.DV[,1:3])
sv.DR<-read.table("sv.DR.txt",header=TRUE)
head(sv.DR[,1:3])
#### remove "BND"-"breakend notations" ###
sv.DR.new<-sv.DR[sv.DR$type!="BND",]
sv.DV.new<-sv.DV[sv.DV$type!="BND",]

sv.DRnew2<-sv.DR.new[,-c(1:3)]
sv.DRnew2[sv.DRnew2!=0]<-1

sv.DVnew2<-sv.DV.new[,-c(1:3)]
sv.DVnew2[sv.DVnew2!=0]<-1

sv.allnew2<-sv.DRnew2+sv.DVnew2
sv.allnew2[sv.allnew2!=0]<-1

LI.DRnew<-sv.DRnew2[,27:33]
LM.DRnew<-sv.DRnew2[,34:39]


LI.DVnew<-sv.DVnew2[,27:33]
LM.DVnew<-sv.DVnew2[,34:39]

#### remove loci that are full heterozygotes in one parental species ###
sv.reads<-sv.DR.new[,c(1:3)]
sv.reads$LI.DR<-rowSums(LI.DRnew)
sv.reads$LM.DR<-rowSums(LM.DRnew)


sv.reads$LI.DV<-rowSums(LI.DVnew)
sv.reads$LM.DV<-rowSums(LM.DVnew)


sv.reads$LI<-sv.reads$LI.DR+sv.reads$LI.DV
sv.reads$LM<-sv.reads$LM.DR+sv.reads$LM.DV

sv.reads$all<-rowSums(sv.allnew2)

length(sv.reads[sv.reads$LI==14,1])
length(sv.reads[sv.reads$LM==12,1])
length(sv.reads[,1])

sv.readsnew<-sv.reads[sv.reads$LI!=14,]
length(sv.readsnew[,1])

sv.readsnew<-sv.readsnew[sv.readsnew$LM!=12,]
length(sv.readsnew[,1])
#### remove loci that are missing in 20% of individuals ###
sv.readsnew<-sv.readsnew[!(sv.readsnew$all<31),]
length(sv.readsnew[,1])

sv.DR.new2<-sv.DR.new[sv.DR.new$position %in% sv.readsnew$position,]
sv.DV.new2<-sv.DV.new[sv.DV.new$position %in% sv.readsnew$position,]

sv.DR.new2[1:10,1:3]
sv.DV.new2[1:10,1:3]

length(sv.DR.new2[,1])
length(sv.DV.new2[,1])
########################################################################
###### calculating genotype likelihood with number of reads ##########
sv.LK<-sv.DR.new2[,1:3]
p.ref<-1
p.het<-0.5
p.SV<-0
err<-0.01
################################################3
for(i in 1:39){
  LhomRef<-(p.ref*(1-err) + (1-p.ref)*err )^sv.DR.new2[,i+3] * ((1-p.ref)*(1-err)+p.ref*err)^sv.DV.new2[,i+3]
  Lhet<-(p.het*(1-err) + (1-p.het)*err )^sv.DR.new2[,i+3] * ((1-p.het)*(1-err)+p.het*err)^sv.DV.new2[,i+3]  
  LhomSV<-(p.SV*(1-err) + (1-p.SV)*err )^sv.DR.new2[,i+3] * ((1-p.SV)*(1-err)+p.SV*err)^sv.DV.new2[,i+3] 
  L.sum<-data.frame(LhomRef,Lhet,LhomSV) 
  MaxL<-apply(L.sum, 1, max)
  sv.LK[,i*3+1]<-LhomRef/MaxL
  sv.LK[,i*3+2]<-Lhet/MaxL
  sv.LK[,i*3+3]<-LhomSV/MaxL
}

sv.LKnew<--10*log10(sv.LK[,-c(1,2,3)])
sv.LKnew<-data.frame(sv.LK[,1:3],sv.LKnew)
length(sv.LKnew[,1]) ### 290276 loci
sv.LKnew[1:10,1:10]
###################################
### write.table(sv.LKnew,"sv.LKfilternew.txt")
### sv.LKnew<-read.table("sv.LKfilternew.txt")
############################################################
##################################################################################
########### read allele frequency for each locus ###############
allele.fre<-read.table("sv.allele.txt")
sv.LKnew$allele.fre<-allele.fre[,3]

LI.fre<-read.table("LI.allele.txt")
sv.LKnew$LI.fre<-LI.fre[,3]

LM.fre<-read.table("LM.allele.txt")
sv.LKnew$LM.fre<-LM.fre[,3]

sv.LKnew$allele.diff<-abs(sv.LKnew$LI.fre-sv.LKnew$LM.fre)

#### filter 2: remove loci with minor allele frequency (0.05) ###################
sv.LKnew2<-sv.LKnew[!(sv.LKnew$allele.fre<0.05),]
sv.LKnew2<-sv.LKnew2[!(sv.LKnew2$allele.fre>0.95),]
length(sv.LKnew2[,1]) #### 127574 loci ####
### write.table(sv.LKnew2,"sv.LKnew2.txt")
######## get scaffold position for each locus #############
head(sv.LKnew2)
sv2.position<-sv.LKnew2[,1:3]
sv2.position$scaffold<-NA
sv2.position$start<-NA
for (i in 1:length(sv2.position[,1])){
  sv2.position$scaffold[i]<-paste(strsplit(sv2.position$position[i], "_")[[1]][2])
  sv2.position$start[i]<-strsplit(sv2.position$position[i], "_")[[1]][5]
}
head(sv2.position)
sv2.position$start<-as.numeric(sv2.position$start)
sv2.position$end<-sv2.position$start+sv2.position$length
####### match chromosome number to scaffolds #########
lgs<-read.table("lgs.txt",header=TRUE)
head(lgs)

for (i in 1:length(lgs[,1])){
  sv2.position[sv2.position$scaffold==lgs$NewScaf[i],"linkage"]<-lgs$LG[i]
}
head(sv2.position)

sv2.position$end<-as.numeric(sv2.position$end)
sv2.position[,8:11]<-sv.LKnew2[,121:124] ###align allele frequencies to each locus ####
#### filter 3: remove loci size<1000 bp (0.05) ###################
sv.LKnew3<-sv.LKnew2[!(sv.LKnew2$length<1000),]
length(sv.LKnew3[,1]) ### 15319 loci ###
sv3.position<-sv.LKnew3[,1:3]
for (i in 1:length(sv3.position[,1])){
  sv3.position$scaffold[i]<-strsplit(sv3.position$position[i], "_")[[1]][2]
  sv3.position$start[i]<-strsplit(sv3.position$position[i], "_")[[1]][5]
}
head(sv3.position)
sv3.position$start<-as.numeric(sv3.position$start)
sv3.position$end<-sv3.position$start+sv3.position$length
sv3.position$linkage<-NA

for (i in 1:length(lgs[,1])){
  sv3.position[sv3.position$scaffold==lgs$NewScaf[i],"linkage"]<-lgs$LG[i]
}
head(sv3.position)

sv3.position$start<-as.numeric(sv3.position$start)
sv3.position$end<-as.numeric(sv3.position$end)
sv3.position[,8:11]<-sv.LKnew3[,121:124]
#### filter 4: only include loci with allele frequency differences>0.3 ###################
sv.LKnew4<-sv.LKnew3[!(sv.LKnew3$allele.diff<0.3),]
length(sv.LKnew4[,1]) ### 1419 loci ###
#### write.table(sv.LKnew4,"sv.LKnew4.txt") 
###############################################################################
################# error rate = 0.05 #################
sv.LK0.05<-sv.DR.new2[,1:3]
p.ref<-1
p.het<-0.5
p.SV<-0
err<-0.05
################################################3
for(i in 1:39){
  LhomRef<-(p.ref*(1-err) + (1-p.ref)*err )^sv.DR.new2[,i+3] * ((1-p.ref)*(1-err)+p.ref*err)^sv.DV.new2[,i+3]
  Lhet<-(p.het*(1-err) + (1-p.het)*err )^sv.DR.new2[,i+3] * ((1-p.het)*(1-err)+p.het*err)^sv.DV.new2[,i+3]  
  LhomSV<-(p.SV*(1-err) + (1-p.SV)*err )^sv.DR.new2[,i+3] * ((1-p.SV)*(1-err)+p.SV*err)^sv.DV.new2[,i+3] 
  L.sum<-data.frame(LhomRef,Lhet,LhomSV) 
  MaxL<-apply(L.sum, 1, max)
  sv.LK0.05[,i*3+1]<-LhomRef/MaxL
  sv.LK0.05[,i*3+2]<-Lhet/MaxL
  sv.LK0.05[,i*3+3]<-LhomSV/MaxL
}

sv.LK0.05new<--10*log10(sv.LK0.05[,-c(1,2,3)])
sv.LK0.05new<-data.frame(sv.LK0.05[,1:3],sv.LK0.05new)
length(sv.LK0.05new[,1]) ### 290276 loci
sv.LK0.05new[1:10,1:10]

write.table(sv.LK0.05new,"sv.LK0.05new.txt")
sv.LK0.05new<-read.table("sv.LK0.05new.txt")
#### filter 2.1: remove loci with minor allele frequency (0.05) ###################
allele.fre_0.05<-read.table("sv0.05_allele.txt")
LI.fre_0.05<-read.table("LI0.05_allele.txt")
LM.fre_0.05<-read.table("LM0.05_allele.txt")
DBS.fre_0.05<-read.table("DBS0.05_allele.txt")
sv.LK0.05new$allele.fre<-allele.fre_0.05[,3]
sv.LK0.05new$allele.diff<-abs(LI.fre_0.05[,3]-LM.fre_0.05[,3])
sv.LK0.05new$LI.fre<-LI.fre_0.05[,3]
sv.LK0.05new$LM.fre<-LM.fre_0.05[,3]
sv.LK0.05new$DBS.fre<-DBS.fre_0.05[,3]
sv.LKnew2.1<-sv.LK0.05new[!(sv.LK0.05new$allele.fre<0.05),]
sv.LKnew2.1<-sv.LKnew2.1[!(sv.LKnew2.1$allele.fre>0.95),]
length(sv.LKnew2.1[,1]) #### 86788 loci ####
#### write.table(sv.LKnew2.1,"sv.LK0.05new2.txt")
head(sv.LKnew2.1)
sv2.1.position<-sv.LKnew2.1[,1:3]
sv2.1.position$scaffold<-NA
sv2.1.position$start<-NA
for (i in 1:length(sv2.1.position[,1])){
  sv2.1.position$scaffold[i]<-strsplit(sv2.1.position$position[i], "_")[[1]][2]
  sv2.1.position$start[i]<-strsplit(sv2.1.position$position[i], "_")[[1]][5]
}
head(sv2.1.position)
sv2.1.position$start<-as.numeric(sv2.1.position$start)
sv2.1.position$end<-sv2.1.position$start+sv2.1.position$length

for (i in 1:length(lgs[,1])){
  sv2.1.position[sv2.1.position$scaffold==lgs$NewScaf[i],"linkage"]<-lgs$LG[i]
}
head(sv2.1.position)

sv2.1.position$end<-as.numeric(sv2.1.position$end)
sv2.1.position[,8:12]<-sv.LKnew2.1[,121:125]
######### write.csv(sv2.1.position,"sv2.1.position.csv") #########
ddply(sv2.1.position,"type",summarize,N=sum(count))
#### filter 3: remove loci size<1000 bp (0.05) ###################
sv.LKnew3.1<-sv.LKnew2.1[!(sv.LKnew2.1$length<1000),]
length(sv.LKnew3.1[,1]) ### 11199 loci ###
sv3.1.position<-sv.LKnew3.1[,1:3]
for (i in 1:length(sv3.1.position[,1])){
  sv3.1.position$scaffold[i]<-strsplit(sv3.1.position$position[i], "_")[[1]][2]
  sv3.1.position$start[i]<-strsplit(sv3.1.position$position[i], "_")[[1]][5]
}
head(sv3.1.position)
sv3.1.position$start<-as.numeric(sv3.1.position$start)
sv3.1.position$end<-sv3.1.position$start+sv3.1.position$length
sv3.1.position$linkage<-NA

for (i in 1:length(lgs[,1])){
  sv3.1.position[sv3.1.position$scaffold==lgs$NewScaf[i],"linkage"]<-lgs$LG[i]
}
head(sv3.1.position)

sv3.1.position$start<-as.numeric(sv3.1.position$start)
sv3.1.position$end<-as.numeric(sv3.1.position$end)
sv3.1.position[,8:11]<-sv.LKnew3.1[,121:124]
head(sv3.1.position)
ddply(sv3.1.position,"type",summarize,N=sum(count))
#### filter 4: only include loci with allele frequency differences>0.3 ###################
sv.LKnew4.1<-sv.LKnew3.1[!(sv.LKnew3.1$allele.diff<0.3),]
length(sv.LKnew4.1[,1]) ### 1403 loci ###
### write.table(sv.LKnew4.1,"sv.LKnew4.1.txt") 
################################################################################
##### figure 1. sv coverage #############
### read the scaffold length ###
scaffold.length<-read.table("scaffold_lengths.txt")
head(scaffold.length)
colnames(scaffold.length)<-c("scaffold","length")
head(lgs)
for(i in 1:length(lgs[,1])){
  lgs$length[i]<-scaffold.length[scaffold.length$scaffold==lgs$NewScaf[i],"length"]
}
LG.length<-ddply(lgs,c("LG"),summarize,length=sum(length))
LG.length$LG<-as.factor(LG.length$LG)
LG.length$LG<-factor(LG.length$LG,levels=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,
                                           15,16,17,18,19,20,21,22,"Z"))
GM.size<-sum(LG.length$length) ###
#################################################################################
########### filter 2, all SVs ########
###### plot figure 1A ########
sv2.position$count<-1
type.sum.filter2<-ddply(sv2.position,c("linkage","type"),summarize, N=sum(count),length=sum(length))
head(type.sum.filter2)
type.sum.filter2$linkage<-as.factor(type.sum.filter2$linkage)
type.sum.filter2$type<-factor(type.sum.filter2$type,levels=c("INV","DUP","INS","DEL"))
type.sum.filter2<-type.sum.filter2[!is.na(type.sum.filter2$linkage),]
type.sum.filter2$linkage<-factor(type.sum.filter2$linkage,levels=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,
                                                                   15,16,17,18,19,20,21,22,"Z"))

sv2.sum<-ddply(sv2.position,c("type"),summarize, N=sum(count),mean.length=sum(length)/N,
      lowerCI=quantile(length,0.025),higherCI=quantile(length,0.975))

ddply(sv2.position,c("type"),summarize, N=sum(count),length=sum(length))[,3]/
  ddply(sv2.position,c("type"),summarize, N=sum(count),length=sum(length))[,2]

coverage.numbersfilter2<-ggplot(data=type.sum.filter2,aes(x=linkage,y=N,fill=type))+
  geom_bar(colour="black",stat="identity")+
  ylab("Number of SVs")+xlab("")+
  scale_fill_manual(values=c("black","white","grey90","grey60"))+
  scale_y_continuous(expand = c(0,0),limits=c(0,10000),breaks=seq(0,10000,2500))+
  theme_bw()+
  theme(
    plot.margin = unit(c(0.8,0.8,0.5,0.8),"cm"),
    axis.title.x = element_text(size = 36),
    axis.title.y = element_text(size = 36),
    axis.text.y = element_text(colour = 'black',size=32),
    axis.text.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(colour = "black",size=2),
    legend.text = element_text(size = 36),legend.spacing.y = unit(0.6, 'cm'),
    legend.title  = element_blank(),legend.position=c(0.8,0.85))
coverage.numbersfilter2
################ plot figure 1B chromosome coverage of SVs ###############################
#### calculating the overall length covered by SVs #####
coverage.length.filter2<-c()
sv2.position$linkage<-as.factor(sv2.position$linkage)
sv2.position$linkage<-factor(sv2.position$linkage,levels=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,
                                                           15,16,17,18,19,20,21,22,"Z"))
linkage<-levels(as.factor(sv2.position$linkage))
sv2.positionnew<-sv2.position[!is.na(sv2.position$linkage),]
coverage.length.filter2<-c()
for (i in 1:23){
  chromo1<-sv2.positionnew[sv2.positionnew$linkage==linkage[i],]
  chromo1<-chromo1[order(chromo1$start,chromo1$end),]
  chromo1$overlap[1]<-TRUE
  chromo1$segment[1]<-1
  for(j in 2:length(chromo1[,1])){
    chromo1$overlap[j]<-(chromo1$start[j]<chromo1$end[j-1])
    chromo1$segment[j]<-chromo1$segment[j-1]+length(which(!((chromo1$start[j]<chromo1$end[j-1]))))
  }
  contig<-data.frame(matrix(ncol=3,nrow=max(chromo1$segment)))
  contig[,1]<-1:max(chromo1$segment)
  colnames(contig)<-c("segment","start","end")
  for(p in 1:max(chromo1$segment)){
    segment<-chromo1[chromo1$segment==p,]
    contig[p,2]<-min(segment$start)
    contig[p,3]<-max(segment$end)
  }
  
  contig$length<-contig$end-contig$start+1
  coverage.length.filter2[i]<-sum(contig$length)
}

sum(coverage.length.filter2)/GM.size ### totoal SV coverage ###
### calculate the length by each SV type #####
coverage.SV.length.filter2<-data.frame(matrix(ncol=3,nrow=23*4))
coverage.SV.length.filter2[,1]<-rep(linkage, each=4)
coverage.SV.length.filter2[,2]<-rep(unique(sv2.positionnew$type),times=23)
colnames(coverage.SV.length.filter2)<-c("linkage","type","length")
for (i in 1:23){
  for (j in 1:4){
    chromo1<-sv2.positionnew[sv2.positionnew$linkage==linkage[i],]
    chromo1<-chromo1[chromo1$type==unique(sv2.positionnew$type)[j],]
    chromo1<-chromo1[order(chromo1$start,chromo1$end),]
    chromo1$overlap[1]<-TRUE
    chromo1$segment[1]<-1
    for(p in 2:length(chromo1[,1])){
      chromo1$overlap[p]<-(chromo1$start[p]<chromo1$end[p-1])
      chromo1$segment[p]<-chromo1$segment[p-1]+length(which(!((chromo1$start[p]<chromo1$end[p-1]))))
    }
    contig<-data.frame(matrix(ncol=3,nrow=max(chromo1$segment)))
    contig[,1]<-1:max(chromo1$segment)
    colnames(contig)<-c("segment","start","end")
    for(q in 1:max(chromo1$segment)){
      segment<-chromo1[chromo1$segment==q,]
      contig[q,2]<-min(segment$start)
      contig[q,3]<-max(segment$end)
    }
    contig$length<-contig$end-contig$start+1
    coverage.SV.length.filter2[coverage.SV.length.filter2$linkage==linkage[i] & coverage.SV.length.filter2$type==unique(sv2.position$type)[j],"length"]<-sum(contig$length)
  }
}

ggplot(data=coverage.SV.length.filter2,aes(x=linkage,y=length,fill=type))+
  geom_bar(stat="identity")
#################### coverage of SVs per type ########################
coverage.SV.length.filter2$count<-1
 ddply(coverage.SV.length.filter2,"type",summarize,length=sum(length))[,2]/
  GM.size
#############################################################
###### figure 1B ##########
for (i in 1:23){
  coverage.SV.length.filter2$LG.length[coverage.SV.length.filter2$linkage==LG.length$LG[i]]<-LG.length$length[i]
}
coverage.SV.length.filter2$proportion<-coverage.SV.length.filter2$length/coverage.SV.length.filter2$LG.length

coverage.SV.length.filter2$linkage<-as.factor(coverage.SV.length.filter2$linkage)
coverage.SV.length.filter2$linkage<-factor(coverage.SV.length.filter2$linkage,levels=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,
                                                                                       15,16,17,18,19,20,21,22,"Z"))
for(i in 1:23){
  temp<-coverage.SV.length.filter2[coverage.SV.length.filter2$linkage==LG.length$LG[i],]
  temp$cumu.percentage[1]<-temp$proportion[1]
  for (j in 2:4){
    temp$cumu.percentage[j]<-temp$proportion[j]+temp$cumu.percentage[j-1]
  }
  coverage.SV.length.filter2$cumu.percentage[coverage.SV.length.filter2$linkage==LG.length$LG[i]]<-temp$cumu.percentage
}
coverage.SV.length.filter2$percentage<-paste(round(coverage.SV.length.filter2$proportion*100,digits = 1),"%",sep="")
coverage.SV.length.filter2$type<-factor(coverage.SV.length.filter2$type,levels=c("INV","DUP","INS","DEL"))
coverage.lengthfilter2<-ggplot(data=coverage.SV.length.filter2,aes(x=linkage,y=proportion,fill=type,label=percentage))+
  geom_bar(color="black",stat="identity")+
  ylab("Proportion of chromosome covered by SVs")+xlab("")+
  scale_fill_manual(values=c("black","white","grey90","grey60"))+
  scale_y_continuous(expand=c(0,0),limits=c(0,0.6),breaks=seq(0,0.6,by=0.1))+
  theme_bw()+
  theme(
    plot.margin = unit(c(0.8,0.8,0.5,0.8),"cm"),
    axis.title.x = element_text(size = 36),
    axis.title.y = element_text(size = 36),
    axis.text.y = element_text(colour = 'black',size=32),
    axis.text.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(colour = "black",size=2),
    legend.position="none")
coverage.lengthfilter2
###### filter 3, SVs>1000 bp ####
sv3.position$count<-1
type.sum.filter3<-ddply(sv3.position,c("linkage","type"),summarize, N=sum(count),length=sum(length))
head(type.sum.filter3)
type.sum.filter3$linkage<-as.factor(type.sum.filter3$linkage)
type.sum.filter3$type<-factor(type.sum.filter3$type,levels=c("INV","DUP","INS","DEL"))
type.sum.filter3<-type.sum.filter3[!is.na(type.sum.filter3$linkage),]
type.sum.filter3$linkage<-factor(type.sum.filter3$linkage,levels=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,
                                                                   15,16,17,18,19,20,21,22,"Z"))
##### summary of SV number and range of length #######
sv3.sum<-ddply(sv3.position,c("type"),summarize, N=sum(count),mean.length=sum(length)/N,
               lowerCI=quantile(length,0.025),higherCI=quantile(length,0.975))

sv.sum<-list(sv2.sum,sv3.sum)
### write.csv(sv.sum,"sv.sum.csv")
####### plot 1C ############
coverage.numbersfilter3<-ggplot(data=type.sum.filter3,aes(x=linkage,y=N,fill=type))+
  geom_bar(colour="black",stat="identity")+
  ylab("Number of large SVs")+xlab("Chromosome number")+
  scale_fill_manual(values=c("black","white","grey90","grey60"))+
  scale_y_continuous(expand = c(0,0),limits=c(0,1250),breaks=seq(0,1250,250))+
  theme_bw()+
  theme(
    plot.margin = unit(c(1.2,0.8,0.5,0.8),"cm"),
    axis.title.x = element_text(size = 36),
    axis.title.y = element_text(size = 36),
    axis.text.y = element_text(colour = 'black',size=32),
    axis.text.x = element_text(colour = 'black',size=28),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(colour = "black",size=2),
    legend.position = "none")
coverage.numbersfilter3
#### calculating the overall length covered by each SV type #####
coverage.length.filter3<-c()
sv3.position$linkage<-as.factor(sv3.position$linkage)
sv3.position$linkage<-factor(sv3.position$linkage,levels=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,
                                                           15,16,17,18,19,20,21,22,"Z"))
linkage<-levels(as.factor(sv3.position$linkage))
sv3.positionnew<-sv3.position[!is.na(sv3.position$linkage),]
coverage.length.filter3<-c()
for (i in 1:23){
  chromo1<-sv3.positionnew[sv3.positionnew$linkage==linkage[i],]
  chromo1<-chromo1[order(chromo1$start,chromo1$end),]
  chromo1$overlap[1]<-TRUE
  chromo1$segment[1]<-1
  for(j in 2:length(chromo1[,1])){
    chromo1$overlap[j]<-(chromo1$start[j]<chromo1$end[j-1])
    chromo1$segment[j]<-chromo1$segment[j-1]+length(which(!((chromo1$start[j]<chromo1$end[j-1]))))
  }
  contig<-data.frame(matrix(ncol=3,nrow=max(chromo1$segment)))
  contig[,1]<-1:max(chromo1$segment)
  colnames(contig)<-c("segment","start","end")
  for(p in 1:max(chromo1$segment)){
    segment<-chromo1[chromo1$segment==p,]
    contig[p,2]<-min(segment$start)
    contig[p,3]<-max(segment$end)
  }
  
  contig$length<-contig$end-contig$start+1
  coverage.length.filter3[i]<-sum(contig$length)
}

sum(coverage.length.filter3)/GM.size #### SV coverage ####
### calculate the length by each SV type #####
coverage.SV.length.filter3<-data.frame(matrix(ncol=3,nrow=23*4))
coverage.SV.length.filter3[,1]<-rep(linkage, each=4)
coverage.SV.length.filter3[,2]<-rep(unique(sv3.positionnew$type),times=23)
colnames(coverage.SV.length.filter3)<-c("linkage","type","length")
for (i in 1:23){
  for (j in 1:4){
    chromo1<-sv3.positionnew[sv3.positionnew$linkage==linkage[i],]
    chromo1<-chromo1[chromo1$type==unique(sv3.positionnew$type)[j],]
    chromo1<-chromo1[order(chromo1$start,chromo1$end),]
    chromo1$overlap[1]<-TRUE
    chromo1$segment[1]<-1
    if (length(chromo1[,1])>1){
    for(p in 2:length(chromo1[,1])){
      chromo1$overlap[p]<-(chromo1$start[p]<chromo1$end[p-1])
      chromo1$segment[p]<-chromo1$segment[p-1]+length(which(!((chromo1$start[p]<chromo1$end[p-1]))))
    }
    }
    contig<-data.frame(matrix(ncol=3,nrow=max(chromo1$segment)))
    contig[,1]<-1:max(chromo1$segment)
    colnames(contig)<-c("segment","start","end")
    for(q in 1:max(chromo1$segment)){
      segment<-chromo1[chromo1$segment==q,]
      contig[q,2]<-min(segment$start)
      contig[q,3]<-max(segment$end)
    }
    contig$length<-contig$end-contig$start+1
    coverage.SV.length.filter3[coverage.SV.length.filter3$linkage==linkage[i] & coverage.SV.length.filter3$type==unique(sv3.position$type)[j],"length"]<-sum(contig$length)
  }
}

ggplot(data=coverage.SV.length.filter3,aes(x=linkage,y=length,fill=type))+
  geom_bar(stat="identity")
################## plot 1D #####################
### read the scaffold length ###
for (i in 1:23){
  coverage.SV.length.filter3$LG.length[coverage.SV.length.filter3$linkage==LG.length$LG[i]]<-LG.length$length[i]
}
coverage.SV.length.filter3$proportion<-coverage.SV.length.filter3$length/coverage.SV.length.filter3$LG.length
coverage.SV.length.filter3$linkage<-as.factor(coverage.SV.length.filter3$linkage)
coverage.SV.length.filter3$linkage<-factor(coverage.SV.length.filter3$linkage,levels=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,
                                                                                       15,16,17,18,19,20,21,22,"Z"))
for(i in 1:23){
  temp<-coverage.SV.length.filter3[coverage.SV.length.filter3$linkage==LG.length$LG[i],]
  temp$cumu.percentage[1]<-temp$proportion[1]
  for (j in 2:4){
    temp$cumu.percentage[j]<-temp$proportion[j]+temp$cumu.percentage[j-1]
  }
  coverage.SV.length.filter3$cumu.percentage[coverage.SV.length.filter3$linkage==LG.length$LG[i]]<-temp$cumu.percentage
}
coverage.SV.length.filter3$percentage<-paste(round(coverage.SV.length.filter3$proportion*100,digits = 1),"%",sep="")
coverage.SV.length.filter3$type<-factor(coverage.SV.length.filter3$type,levels=c("INV","DUP","INS","DEL"))
coverage.lengthfilter3<-ggplot(data=coverage.SV.length.filter3,aes(x=linkage,y=proportion,fill=type,label=percentage))+
  geom_bar(color="black",stat="identity")+
  ylab("Proportion of chromosome covered by large SVs")+xlab("Chromosome number")+
  scale_fill_manual(values=c("black","white","grey90","grey60"))+
  scale_y_continuous(expand=c(0,0),limits=c(0,0.6),breaks=seq(0,0.6,by=0.1))+
  theme_bw()+
  theme(
    plot.margin = unit(c(1.2,0.8,0.5,0.8),"cm"),
    axis.title.x = element_text(size = 36),
    axis.title.y = element_text(size = 36),
    axis.text.y = element_text(colour = 'black',size=32),
    axis.text.x = element_text(colour = 'black',size=28),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(colour = "black",size=2),
    legend.position="none")
coverage.lengthfilter3
### combine four plots ####
plot_grid(coverage.numbersfilter2,coverage.lengthfilter2,
          coverage.numbersfilter3,coverage.lengthfilter3,
          nrow=2,labels=c("A","B","C","D"),align="v",label_size=44)
#### save as 30*26 ####
############################################################################
sv2.position$minor.fre<-sv2.position$allele.fre
sv2.position$minor.fre[sv2.position$allele.fre>0.5]<-1-sv2.position$allele.fre[sv2.position$allele.fre>0.5]

sv2.position$LI.minor<-sv2.position$LI.fre
sv2.position$LI.minor[sv2.position$LI.fre>0.5]<-1-sv2.position$LI.fre[sv2.position$LI.fre>0.5]

sv2.position$LM.minor<-sv2.position$LM.fre
sv2.position$LM.minor[sv2.position$LM.fre>0.5]<-1-sv2.position$LM.fre[sv2.position$LM.fre>0.5]
##### figure S1. Histograms of distribution of allele frequency differences for different types of SVs ####
allele.diff.type<-ggplot(data=sv2.position,aes(allele.diff))+
  facet_rep_grid(~type)+
  geom_histogram(aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]),binwidth = 0.02)+
  scale_y_continuous(expand=c(0,0),limits = c(0,0.15),breaks=seq(0,0.15,by=0.05))+
  scale_x_continuous(expand=c(0,0),limits = c(-0.05,1.08),breaks=seq(0,1,by=0.2))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),  
        panel.background=element_blank(),
        axis.line=element_line(colour="black"),
        text=element_text(colour="black",size=32),
        axis.text=element_text(colour="black",size=28),
        axis.title.y = element_text(colour="black",size=36),
        strip.text = element_text(size=32))+
  xlab("Allele frequency difference ")+ylab("Frequency")
allele.diff.type
##### figure 2. correlation between length of SVs and minor allele frequency  #####
SV.length.type<-ggplot(data=sv2.position,aes(log(length)))+
  facet_rep_grid(~type)+
  geom_histogram(aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]),binwidth =1)+
  scale_y_continuous(expand=c(0,0),limits = c(0,0.41),breaks=seq(0,0.4,by=0.1))+
  scale_x_continuous(expand=c(0,0),limits = c(3,16),breaks=seq(3,15,by=3))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),  
        panel.background=element_blank(),
        axis.line=element_line(colour="black"),
        text=element_text(colour="black",size=32),
        axis.text=element_text(colour="black",size=30),
        axis.title= element_text(colour="black",size=36),
        strip.text = element_text(size=34))+
  xlab("Length of SV (log)")+ylab("Frequency")
SV.length.type

minor.type<-ggplot(data=sv2.position,aes(minor.fre))+
  facet_rep_grid(~type)+
  geom_histogram(aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]))+
  scale_y_continuous(expand=c(0,0),limits = c(0,0.21),breaks=seq(0,0.2,by=0.1))+
  scale_x_continuous(expand=c(0,0),limits = c(0.05,0.53),breaks=seq(0.1,0.5,by=0.1))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),  
        panel.background=element_blank(),
        axis.line=element_line(colour="black"),
        text=element_text(colour="black",size=32),
        axis.text=element_text(colour="black",size=30),
        axis.title = element_text(colour="black",size=36),
        strip.text = element_blank(),
        strip.background =element_blank())+
  xlab("Minor allele frequency")+ylab("Frequency")
minor.type

plot(sort(sv2.position$minor.fre))

plot_grid(SV.length.type,minor.type,
          nrow=2,align="v",labels=c("A","B"),label_size=38)

### write.table(sv3.position,"sv3.position.txt") ###
##########################################################################################
########################### Fst calculation ###############################
N<-39
L<-127574 ## no. SNPs #

DBSf2.af<-read.table("DBS_f2.allele.txt")
head(DBSf2.af)
DBS.af<-DBSf2.af[,c(1,3)]

LIf2.af<-read.table("LI_f2.allele.txt")
head(LIf2.af)
LI.af<-LIf2.af[,c(1,3)]

LMf2.af<-read.table("LM_f2.allele.txt")
head(LMf2.af)
LM.af<-LMf2.af[,c(1,3)]
####### calculation of Fst between different populations ##########
## population allele frequencies
## allele frequency through species
Npop <- 3 ## number of pops ##
Pmat.filter2 <- matrix(NA, nrow = Npop, ncol = L)

Pmat.filter2[1, ] <- DBS.af[,2]
Pmat.filter2[2, ] <- LI.af[,2]
Pmat.filter2[3, ] <- LM.af[,2]

Pmat.filter2[1:3,1:6]
## Fst
pimat.filter2 <- 2 * Pmat.filter2 * (1 - Pmat.filter2)
pimat.filter2[1:3,1:6]
## pairwise Fst
calcPair <- function(pia = NA, pib = NA, ppbar = NA) {
  piT <- 2 * ppbar * (1 - ppbar)
  piS <- (pia + pib)/2
  Fst <- mean(piT - piS)/mean(piT)
  return(Fst)
}

piT<-2 * ((Pmat.filter2[2,] + Pmat.filter2[3, ])/2) * (1 - (Pmat.filter2[2,] + Pmat.filter2[3, ])/2)

pairFst.filter2 <- matrix(0, nrow = Npop, ncol = Npop)
for (i in 1:(Npop - 1)) {
  for (j in (i + 1):Npop) {
    ## only loop over upper-triange, its symmetric pass the function pi for
    ## the two populations, and their mean allele freq.
    pairFst.filter2[i, j] <- calcPair(pia = pimat.filter2[i, ], pib = pimat.filter2[j, ], ppbar = (Pmat.filter2[i,] + Pmat.filter2[j, ])/2)
    pairFst.filter2[j, i] <- pairFst.filter2[i, j]
  }
}

colnames(pairFst.filter2) <- c("DBS","LI","LM")
rownames(pairFst.filter2) <- c("DBS","LI","LM")
pairFst.filter2
#### write.table(pairFst.filter2,"pairFst.filter2.txt")
############# calculation of Fst for each locus ##################
calcPair1 <- function(pia = NA, pib = NA, ppbar = NA) {
  piT <- 2 * ppbar * (1 - ppbar)
  piS <- (pia + pib)/2
  Fst <- (piT - piS)/piT
  return(Fst)
}

Fst.filter2new<-calcPair1(pia = pimat.filter2[2, ], pib = pimat.filter2[3, ], ppbar = (Pmat.filter2[2,] + Pmat.filter2[3, ])/2)
sv.LKnew2$Fst<-Fst.filter2new
### write.table(sv.LKnew2,"sv.LKfilter2.txt")
#############################################################
########## calculation of Fst for each SV types #######
Pmat.filter2[1:3,1:10]
Pmat.filter2t<-t(Pmat.filter2)
length(Pmat.filter2t[,1])
Pmat.filter2t<-data.frame(Pmat.filter2t)
Pmat.filter2t$type<-sv.LKnew2$type
head(Pmat.filter2t)
###########################################################
Pmat.filter2.INV<-Pmat.filter2t[Pmat.filter2t$type=="INV",]
Pmat.filter2.INV<-as.matrix(Pmat.filter2.INV[,-4])
Pmat.filter2.INV<-t(Pmat.filter2.INV)
Pmat.filter2.INV[1:3,1:10]

## Fst
pimat.filter2.INV <- 2 * Pmat.filter2.INV * (1 - Pmat.filter2.INV)
## pairwise Fst
pairFst.filter2.INV <- matrix(0, nrow = Npop, ncol = Npop)
for (i in 1:(Npop - 1)) {
  for (j in (i + 1):Npop) {
    ## only loop over upper-triange, its symmetric pass the function pi for
    ## the two populations, and their mean allele freq.
    pairFst.filter2.INV[i, j] <- calcPair(pia = pimat.filter2.INV[i, ], pib = pimat.filter2.INV[j, ],
                                          ppbar = (Pmat.filter2.INV[i,] + Pmat.filter2.INV[j, ])/2)
    pairFst.filter2.INV[j, i] <- pairFst.filter2.INV[i, j]
  }
}
colnames(pairFst.filter2.INV) <- c("DBS","LI","LM")
rownames(pairFst.filter2.INV) <- c("DBS","LI","LM")
pairFst.filter2.INV

######################################################################
Pmat.filter2.DUP<-Pmat.filter2t[Pmat.filter2t$type=="DUP",]
Pmat.filter2.DUP<-as.matrix(Pmat.filter2.DUP[,-4])
Pmat.filter2.DUP<-t(Pmat.filter2.DUP)
Pmat.filter2.DUP[1:3,1:10]

## Fst
pimat.filter2.DUP <- 2 * Pmat.filter2.DUP * (1 - Pmat.filter2.DUP)
## pairwise Fst
pairFst.filter2.DUP <- matrix(0, nrow = Npop, ncol = Npop)
for (i in 1:(Npop - 1)) {
  for (j in (i + 1):Npop) {
    ## only loop over upper-triange, its symmetric pass the function pi for
    ## the two populations, and their mean allele freq.
    pairFst.filter2.DUP[i, j] <- calcPair(pia = pimat.filter2.DUP[i, ], pib = pimat.filter2.DUP[j, ], 
                                          ppbar = (Pmat.filter2.DUP[i,] + Pmat.filter2.DUP[j, ])/2)
    pairFst.filter2.DUP[j, i] <- pairFst.filter2.DUP[i, j]
  }
}
colnames(pairFst.filter2.DUP) <- c("DBS","LI","LM")
rownames(pairFst.filter2.DUP) <- c("DBS","LI","LM")
pairFst.filter2.DUP
###################################################################################
Pmat.filter2.INS<-Pmat.filter2t[Pmat.filter2t$type=="INS",]
Pmat.filter2.INS<-as.matrix(Pmat.filter2.INS[,-4])
Pmat.filter2.INS<-t(Pmat.filter2.INS)
Pmat.filter2.INS[1:3,1:10]

## Fst
pimat.filter2.INS <- 2 * Pmat.filter2.INS * (1 - Pmat.filter2.INS)
## pairwise Fst
pairFst.filter2.INS <- matrix(0, nrow = Npop, ncol = Npop)
for (i in 1:(Npop - 1)) {
  for (j in (i + 1):Npop) {
    ## only loop over upper-triange, its symmetric pass the function pi for
    ## the two populations, and their mean allele freq.
    pairFst.filter2.INS[i, j] <- calcPair(pia = pimat.filter2.INS[i, ], pib = pimat.filter2.INS[j, ], 
                                          ppbar = (Pmat.filter2.INS[i,] + Pmat.filter2.INS[j, ])/2)
    pairFst.filter2.INS[j, i] <- pairFst.filter2.INS[i, j]
  }
}
colnames(pairFst.filter2.INS) <- c("DBS","LI","LM")
rownames(pairFst.filter2.INS) <- c("DBS","LI","LM")
pairFst.filter2.INS
######################################################################################
Pmat.filter2.DEL<-Pmat.filter2t[Pmat.filter2t$type=="DEL",]
Pmat.filter2.DEL<-as.matrix(Pmat.filter2.DEL[,-4])
Pmat.filter2.DEL<-t(Pmat.filter2.DEL)
Pmat.filter2.DEL[1:3,1:10]

## Fst
pimat.filter2.DEL <- 2 * Pmat.filter2.DEL * (1 - Pmat.filter2.DEL)
## pairwise Fst
pairFst.filter2.DEL <- matrix(0, nrow = Npop, ncol = Npop)
for (i in 1:(Npop - 1)) {
  for (j in (i + 1):Npop) {
    ## only loop over upper-triange, its symmetric pass the function pi for
    ## the two populations, and their mean allele freq.
    pairFst.filter2.DEL[i, j] <- calcPair(pia = pimat.filter2.DEL[i, ], pib = pimat.filter2.DEL[j, ], 
                                          ppbar = (Pmat.filter2.DEL[i,] + Pmat.filter2.DEL[j, ])/2)
    pairFst.filter2.DEL[j, i] <- pairFst.filter2.DEL[i, j]
  }
}
colnames(pairFst.filter2.DEL) <- c("DBS","LI","LM")
rownames(pairFst.filter2.DEL) <- c("DBS","LI","LM")
pairFst.filter2.DEL

pairFst.filter2all<-list(pairFst.filter2,pairFst.filter2.INV,pairFst.filter2.DUP,pairFst.filter2.INS,pairFst.filter2.DEL)
names(pairFst.filter2all)<-c("all","INV","DUP","INS","DEL")
### write.csv(pairFst.filter2all,"pairFst.filter2all.csv")
##############################################################################################
