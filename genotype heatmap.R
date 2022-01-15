########## loading packages ########
library(ggplot2)
library(lemon)
####### genotypes for ancestral SVs ######
genotype4<-read.table("SVfilter4_gprobk2.txt",header=TRUE,sep=",")
genotype4[1:10,1:10]
length(genotype4[,1])
genotype4new<-t(genotype4[,-1])
genotype4new[1:10,1:10]

###### information of ancestral SVs ########
SV.LKnew4<-read.table("sv.LKnew4.txt",header=TRUE)
length(SV.LKnew4[,1])
SV.LKnew4[1:10,1:10]
genotype4new<-data.frame(genotype4new)
genotype4new<-round(genotype4new)
genotype4new$type<-SV.LKnew4$type
genotype4new$length<-SV.LKnew4$length
genotype4new$allele.diff<-SV.LKnew4$allele.diff
genotype4new$position<-SV.LKnew4$position

genotypenew<-genotype4new[genotype4new$allele.diff>=0.5,]
length(genotypenew[,1])
genotypenew$length<-as.numeric(genotypenew$length)
for (i in 1:231){
  genotypenew$scaffold[i]<-strsplit(genotypenew$position[i], "_")[[1]][2]
  genotypenew$start[i]<-strsplit(genotypenew$position[i], "_")[[1]][5]
}

mappos<-read.table("lgs.txt", header=T)
colnames(mappos)[1]<-"scaffold"
genotypenew$locus<-rownames(genotypenew)
genotypenew<-merge(mappos, genotypenew,by="scaffold")
genotypenew$LG<-factor(genotypenew$LG,levels=c("1","2","3","4","5","6","7","8",
                                               "9","10","11","12","13","14","15",
                                               "16","17","18","19","20","21","22","Z"))
genotypenew$start<-as.numeric(genotypenew$start)
genotypenew<-genotypenew[order(genotypenew$LG,genotypenew$start,genotypenew$type,genotypenew$length),]

genotypenew2<-as.matrix(genotypenew[,-c(1,2,42,43,44,45,46,47)])
ind<-c("DBS-13-04","DBS-13-08","DBS-13-09","DBS-14-02","DBS-14-06","DBS-14-08", 
       "DBS-14-15","DBS-14-17","DBS-16-07","DBS-16-08","DBS-16-10","DBS-16-14", 
       "DBS-16-15","DBS-16-16","DBS-16-18","DBS-16-19","DBS-16-20","DBS-16-21",
       "DBS-16-22","DBS-16-23","DBS-16-24","DBS-16-26","DBS-16-27","DBS-16-28",
       "DBS-18-14","DBS-18-15","frc-16-01","frc-16-12","frc-16-15","frc-16-16",
       "mrf-20-01","mrf-20-03","mrf-20-04","SIN-20-01","SIN-20-02","SIN-20-04",
       "SIN-20-05","vic-17-01","vic-18-12")

ind2<-c("frc-16-01","frc-16-12","frc-16-15","frc-16-16",
        "mrf-20-01","mrf-20-03","mrf-20-04","DBS-13-04","DBS-13-08","DBS-13-09","DBS-14-02","DBS-14-06","DBS-14-08", 
        "DBS-14-15","DBS-14-17","DBS-16-07","DBS-16-08","DBS-16-10","DBS-16-14", 
        "DBS-16-15","DBS-16-16","DBS-16-18","DBS-16-19","DBS-16-20","DBS-16-21",
        "DBS-16-22","DBS-16-23","DBS-16-24","DBS-16-26","DBS-16-27","DBS-16-28",
        "DBS-18-14","DBS-18-15","SIN-20-01","SIN-20-02","SIN-20-04",
        "SIN-20-05","vic-17-01","vic-18-12")

hyb.ind<-c("DBS-13-04","DBS-13-08","DBS-13-09","DBS-14-02","DBS-14-06","DBS-14-08", 
           "DBS-14-15","DBS-14-17","DBS-16-07","DBS-16-08","DBS-16-10","DBS-16-14", 
           "DBS-16-15","DBS-16-16","DBS-16-18","DBS-16-19","DBS-16-20","DBS-16-21",
           "DBS-16-22","DBS-16-23","DBS-16-24","DBS-16-26","DBS-16-27","DBS-16-28",
           "DBS-18-14","DBS-18-15")

SV.ids<-c("DBS","DBS","DBS","DBS","DBS","DBS", 
          "DBS","DBS","DBS","DBS","DBS","DBS", 
          "DBS","DBS","DBS","DBS","DBS","DBS", 
          "DBS","DBS","DBS","DBS","DBS","DBS", "DBS","DBS","JH","JH","JH","JH","JH","JH","JH","LM","LM","LM",
          "LM","LM","LM")

SV.ids2<-c("JH","JH","JH","JH","JH","JH","JH","DBS","DBS","DBS",
           "DBS","DBS","DBS","DBS","DBS","DBS","DBS","DBS","DBS", 
           "DBS","DBS","DBS","DBS","DBS","DBS", "DBS","DBS","DBS",
           "DBS","DBS","DBS", "DBS","DBS","LM","LM","LM","LM","LM","LM")

genotypenew3<-data.frame(t(genotypenew2))
genotypenew3$species<-SV.ids
length(genotypenew3[1,])
LM.genotype<-colMeans(genotypenew3[genotypenew3$species=="LM",-233])
LI.genotype<-colMeans(genotypenew3[genotypenew3$species=="JH",-233])
geno<-data.frame(t(rbind(LM.genotype,LI.genotype)))
geno$change<-"N"
geno$change[geno$LI.genotype<geno$LM.genotype]<-"Y"

genotypenew2<-data.frame(genotypenew2)
genotypenew2$change<-geno$change

genotypenew2[genotypenew2$change=="Y",-40]<-2-genotypenew2[genotypenew2$change=="Y",-40]
genotypenew2<-as.matrix(genotypenew2[,-40])

colnames(genotypenew2)<-ind
length(genotypenew[,1])

hybrid.ind<-data.frame(ind)
for (i in 1:39){
  hybrid.ind$index[i]<-mean(genotypenew2[,i])/2
}
hybrid.ind$species<-SV.ids
hybrid.ind$species<-factor(hybrid.ind$species,levels=c("LM","DBS","JH"))
hybrid.indnew<-hybrid.ind[order(hybrid.ind$species,hybrid.ind$index,decreasing=TRUE),]

genotypenew$chrome<-"Autosomes"
genotypenew$chrome[genotypenew$LG=="Z"]<-"Z chromosome"

genotypedata<-data.frame(matrix(nrow=9048,ncol=5))
genotypedata[,1]<-rep(genotypenew$locus,39)
genotypedata[,2]<-rep(ind,each=232)
genotypedata[,3]<-as.vector(genotypenew2)
genotypedata[,4]<-rep(genotypenew$type,39)
genotypedata[,5]<-rep(genotypenew$chrome,39)
genotypedata$species<-rep(SV.ids,each=232)
colnames(genotypedata)<-c("Locus","Individual","genotype","SV.type","chromo","species")
genotypedata$genotype<-as.factor(genotypedata$genotype)
genotypedata$Individual<-factor(genotypedata$Individual,levels=hybrid.indnew$ind)
genotypedata$Locus<-factor(genotypedata$Locus,levels=genotypenew$locus)
levels(genotypedata$Individual)

heatmap<-ggplot(genotypedata) +
  facet_rep_grid(~chromo, scales = "free_x", space = "free_x")+
  geom_tile(aes(x =Locus, y = Individual, fill = genotype))+
  scale_fill_manual(values=c("lightgreen","forestgreen","darkgreen"),guide = guide_legend(),
                    labels=c("Homozygous reference","Heterozygous","Homozygous variant"))+
  theme_bw()+theme(axis.title = element_text(size=42),
                   axis.title.x = element_blank(),
                   axis.text.y = element_text(size=32),
                   axis.text.x= element_blank(),
                   axis.ticks= element_blank(),
                   plot.title = element_text(size=40,hjust=0.5),
                   legend.title = element_blank(),
                   legend.text = element_text(size=36),
                   legend.position="bottom",
                   legend.spacing.x = unit(0.8, 'cm'),
                   panel.spacing.x = unit(0, "lines"),
                   strip.text.x = element_text(size = 40),
                   strip.background =element_rect(fill="white",colour = NA))+
  scale_y_discrete(labels=SV.ids2)

heatmap

