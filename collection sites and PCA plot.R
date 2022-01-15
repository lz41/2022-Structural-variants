################ plot Figure 4 #################
library(maps)
library(ggplot2)
library(cowplot)
library(rcartocolor)
###### collection site map ##########
NCstate=map_data("state")
WYstate<-NCstate[NCstate$region=="wyoming",]
location<-read.csv("site coordinates.csv")
head(location)
US.map<-ggplot(NCstate, aes(long, lat))+geom_polygon()+coord_equal()
US.map<-ggplot()+coord_fixed(1.5)
US.map<-US.map+geom_polygon(data=NCstate, aes(x=long, y=lat, group=group),
                            fill="white", col="black", lwd=0.5)
US.map<-US.map+theme_bw()+
  geom_rect(aes(xmin =-113, xmax =-103, ymin =40, ymax=46),
            fill = "transparent", color = "red", size = 1.5)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background=element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.line=element_blank(),axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    plot.margin=unit(c(0.3,0.3,0.3,0.3),"lines"),
    legend.justification = c(1, 1),
    legend.background = element_rect(),
    legend.key.size=unit(2,"cm"),
    legend.position = c(0.95,0.95),
    legend.title=element_blank(),
    legend.text=element_text(face="italic",size=20),
    legend.key=element_blank())
US.map

p<-ggplot(WYstate, aes(long, lat))+geom_polygon()+coord_equal()
p<-ggplot()+coord_fixed(1.5)
p<-p+geom_polygon(data=WYstate, aes(x=long, y=lat, group=group),
                  fill="white", col="black", lwd=0.5)
p<-p+theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background=element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.line=element_blank(),axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    plot.margin=unit(c(0.3,0.3,0.3,0.3),"lines"),
    legend.justification = c(1, 1),
    legend.background = element_rect(),
    legend.key.size=unit(2,"cm"),
    legend.position = c(0.95,0.95),
    legend.title=element_blank(),
    legend.text=element_text(face="italic",size=32),
    legend.key=element_blank())
p

location$Species<-factor(location$Species,levels=c("Jacksonhole","L. melissa","Hybrids"))

############################################################################
pt_colors=c("maroon2", "slateblue2", "olivedrab2", "black")
library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
family<-as.factor(pops$V1)
usa <- map_data("usa")
state_dat<-map_data("state")
map_dat<-rbind(state_dat,usa)

p0 <- ggplot() +
  geom_polygon(data=map_dat,aes(x=long,y=lat,group=group, fill=region),fill="white",color="black", show.legend=FALSE) +
  coord_map("gilbert",xlim=c(-112,-105),ylim=c(40,45)) +
  labs(x=expression("Longitude"*~degree*W), y=expression("Latitude"*~degree*N)) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.margin = unit(c(0.25,0.5,0.5,0.25),'inches'),
        legend.position='none',axis.title.x=element_text(hjust=0.9,size=28),
        axis.title.y=element_text(size=36),
        axis.text=element_text(size=28)) +
  theme(rect = element_blank())
p0 
site.Map<-p0+geom_point(data=location, 
            aes(x=Longitude, y=Latitude,
                shape=factor(Species),fill=factor(Species)),size=7.5)+
  scale_shape_manual(values=c(22,23,21))+
  scale_fill_manual(values=c("blue4","skyblue1","grey"))

gg_inset_map1<- ggdraw() +
  draw_plot(site.Map) +
  draw_plot(US.map, x =0.00, y = -0.07, width = 0.55, height = 0.55)

gg_inset_map1
############################################################################
######### PCA with different filtered dataset ############
################ PCA SNPs ##################
N<-37
SNP.n<-1164 ## no. SNPs #

SNP.gen <- read.table("SNP_gprobk2.txt",header=TRUE,sep=",")

SNP.gen[1:10, 1:10]
ind<-c("dbs-13-04","dbs-13-08","dbs-13-09","dbs-14-02","dbs-14-06","dbs-14-08", 
       "dbs-14-15","dbs-14-17","dbs-16-07","dbs-16-08","dbs-16-10","dbs-16-14", 
       "dbs-16-15","dbs-16-16","dbs-16-18","dbs-16-19","dbs-16-20","dbs-16-21",
       "dbs-16-22","dbs-16-23","dbs-16-24","dbs-16-26","dbs-16-27","dbs-16-28",
       "frc-16-01","frc-16-12","frc-16-15","frc-16-16","mrf-20-01","mrf-20-03",
       "mrf-20-04","SIN-20-01","SIN-20-02","SIN-20-04","SIN-20-05","vic-17-01",
       "vic-18-12")
SNP.ids<-c("dbs","dbs","dbs","dbs","dbs","dbs","dbs","dbs","dbs","dbs","dbs","dbs", 
           "dbs","dbs","dbs","dbs","dbs","dbs","dbs","dbs","dbs","dbs","dbs","dbs",
           "frc","frc","frc","frc","mrf","mrf","mrf","SIN","SIN","SIN","SIN","vic","vic")

SNP.gen<-data.frame(ind=ind,SNP.gen)
SNP.gen1<-SNP.gen[,-c(1,2)]
SNP.gen1[1:5,1:5]
################################################################################
################### Fst for SNPs ################
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
#####
SNP.gmns <- apply(SNP.gen1, 2, mean)
SNP.p <- (apply(SNP.gen1, 2, sum) + 1)/(2 + 2 * N)
SNP.genS <- SNP.gen1
for (i in 1:SNP.n) {
  ## sqrt(p*(1-p)) is the relevant standard deviation
  SNP.genS[, i] <- (SNP.gen1[, i] - SNP.gmns[i])/sqrt(SNP.p[i] * (1 - SNP.p[i]))
}

SNP.pc <- prcomp(SNP.genS, center = FALSE, scale = FALSE)
summary(SNP.pc)
par(pty = "s")
SNP.pc1mn<- tapply(X = SNP.pc$x[, 1], INDEX = ind, mean)
SNP.pc2mn<- tapply(X = SNP.pc$x[, 2], INDEX = ind, mean)
SNP.pc3mn<- tapply(X = SNP.pc$x[, 3], INDEX = ind, mean)
summary(SNP.pc1mn)
summary(SNP.pc2mn)
summary(SNP.pc3mn)
pop.order<-dimnames(SNP.pc1mn)[[1]]
par(mfrow = c(1, 1))
par(pty = "s")
#### PCA plot ####
SNP.ids<-data.frame(SNP.ids)

SNP.ids$fill[SNP.ids$SNP.ids=="dbs"]<-"grey"
SNP.ids$fill[SNP.ids$SNP.ids %in% c("mrf","frc")]<-"white"
SNP.ids$fill[SNP.ids$SNP.ids %in% c("SIN","vic")]<-"black"

SNP.ids$fill[SNP.ids$SNP.ids=="dbs"]<-"grey"
SNP.ids$fill[SNP.ids$SNP.ids %in% c("mrf","frc")]<-"white"
SNP.ids$fill[SNP.ids$SNP.ids %in% c("SIN","vic")]<-"black"

SNP.ids$pch[SNP.ids$SNP.ids=="dbs"]<-21
SNP.ids$pch[SNP.ids$SNP.ids %in% c("mrf","frc")]<-22
SNP.ids$pch[SNP.ids$SNP.ids %in% c("SIN","vic")]<-23
summary(SNP.pc1mn)
summary(SNP.pc2mn)
plot(SNP.pc1mn, SNP.pc2mn, type = "n", xlab = "PC1 (21%)", ylab = "PC2 (10.2%)",xlim=c(-50,50),ylim=c(-30,30))
points(SNP.pc1mn,SNP.pc2mn,bg=SNP.ids$fill,pch=SNP.ids$pch,cex = 1)
###########################################################
#### PCA filter 4 #########
N<-39
L4<-1419 ## no. SNPs #
SV4.gen<-read.table("SVfilter4_gprobk2.txt",sep =",",header = TRUE)
SV4.gen[1:10, 1:10]
SV4.gen<-data.frame(SV4.gen)
ind<-c("dbs-13-04","dbs-13-08","dbs-13-09","dbs-14-02","dbs-14-06","dbs-14-08", 
       "dbs-14-15","dbs-14-17","dbs-16-07","dbs-16-08","dbs-16-10","dbs-16-14", 
       "dbs-16-15","dbs-16-16","dbs-16-18","dbs-16-19","dbs-16-20","dbs-16-21",
       "dbs-16-22","dbs-16-23","dbs-16-24","dbs-16-26","dbs-16-27","dbs-16-28",
       "dbs-18-14","dbs-18-15","frc-16-01","frc-16-12","frc-16-15","frc-16-16",
       "mrf-20-01","mrf-20-03","mrf-20-04","SIN-20-01","SIN-20-02","SIN-20-04",
       "SIN-20-05","vic-17-01","vic-18-12")
SV.ids<-c("dbs","dbs","dbs","dbs","dbs","dbs", 
          "dbs","dbs","dbs","dbs","dbs","dbs", 
          "dbs","dbs","dbs","dbs","dbs","dbs",
          "dbs","dbs","dbs","dbs","dbs","dbs","dbs","dbs","frc","frc","frc","frc","mrf","mrf","mrf","SIN","SIN","SIN",
          "SIN","vic","vic")
ind.data<-data.frame(ind,SV.ids)
ind.data$count<-1
ind.sum<-ddply(ind.data,"SV.ids",summarize,N=sum(count))
write.csv(ind.sum,"ind.sum.csv")
SV4.gen<-data.frame(ind=ind,SV4.gen)
SV4.gen1<-SV4.gen[,-c(1,2)]
SV4.gen1[1:10,1:10]
#####
SV4.gmns <- apply(SV4.gen1, 2, mean)
SV4.p <- (apply(SV4.gen1, 2, sum) + 1)/(2 + 2 * N)
SV4.genS <- SV4.gen1
for (i in 1:L4) {
  ## sqrt(p*(1-p)) is the relevant standard deviation
  SV4.genS[, i] <- (SV4.gen1[, i] - SV4.gmns[i])/sqrt(SV4.p[i] * (1 - SV4.p[i]))
}

SV4.pc <- prcomp(SV4.genS, center = FALSE, scale = FALSE)
summary(SV4.pc)
par(pty = "s")
SV4.pc1mn<- tapply(X = SV4.pc$x[, 1], INDEX = ind, mean)
SV4.pc2mn<- tapply(X = SV4.pc$x[, 2], INDEX = ind, mean)
SV4.pc3mn<- tapply(X = SV4.pc$x[, 3], INDEX = ind, mean)
summary(SV4.pc1mn)
summary(SV4.pc2mn)
summary(SV4.pc3mn)
pop.order<-dimnames(SV4.pc1mn)[[1]]
par(mfrow = c(1, 1))
par(pty = "s")
#### PCA plot ####
SV4.ids<-data.frame(SV.ids)
SV4.ids$color[SV4.ids$SV.ids=="dbs"]<-"green"
SV4.ids$color[SV4.ids$SV.ids %in% c("mrf","frc")]<-"blue"
SV4.ids$color[SV4.ids$SV.ids %in% c("SIN","vic")]<-"red"
plot(SV4.pc1mn, SV4.pc2mn, type = "n", xlab = "PC1 (24.3%)", ylab = "PC2 (8.1%)",xlim=c(-50,40),ylim=c(-20,25))
points(SV4.pc1mn,SV4.pc2mn,col=SV4.ids$color,pch=19,cex = 1)
######## PCA with SNP and SV data ############
SV4.pc1mnnew<-SV4.pc1mn[-c(25,26)]
cor.test(SNP.pc2mn,SV4.pc1mnnew)
PCA.data<-data.frame(SNP.ids,SNP.PC=SNP.pc2mn,SV.PC=SV4.pc1mnnew)
PCA.data$species<-"hybrid"
PCA.data$species[PCA.data$SNP.ids %in% c("frc","mrf")]<-"L.idas"
PCA.data$species[PCA.data$SNP.ids %in% c("SIN","vic")]<-"L.melissa"
PCA.data$species<-factor(PCA.data$species,levels=c("L.idas","L.melissa","hybrid"))


PCAplot<-ggplot(data=PCA.data,aes(x=SNP.pc2mn,y=-SV4.pc1mnnew,fill=species,shape=species))+
  geom_point(size=8)+
  scale_shape_manual(values=c(22,23,21))+
  scale_fill_manual(values=c("blue4","skyblue1","grey"))+
  ylab("PC of SVs (24.3%)")+ xlab("PC of SNPs (10.2%)")+
  scale_y_continuous(limits=c(-40,45),breaks=seq(-40,40,20))+
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background=element_blank(),
        text=element_text(colour="black",size=32),
        axis.text=element_text(colour="black",size=28),
        strip.text = element_text(face = "italic",size=32),legend.position = "none")+
  annotate(geom="text",x=-10,y=40,label=expression(paste(italic(R),"= 0.777, ",italic(P),"< 0.0001")),size=10)

PCAplot
#########################################################################################
############################ entropy plot ##############################
filter4.q2<-read.table("SV_filter4_q2.txt",header=T,sep=",")
ind<-c("dbs-13-04","dbs-13-08","dbs-13-09","dbs-14-02","dbs-14-06","dbs-14-08", 
       "dbs-14-15","dbs-14-17","dbs-16-07","dbs-16-08","dbs-16-10","dbs-16-14", 
       "dbs-16-15","dbs-16-16","dbs-16-18","dbs-16-19","dbs-16-20","dbs-16-21",
       "dbs-16-22","dbs-16-23","dbs-16-24","dbs-16-26","dbs-16-27","dbs-16-28",
       "dbs-18-14","dbs-18-15","frc-16-01","frc-16-12","frc-16-15","frc-16-16",
       "mrf-20-01","mrf-20-03","mrf-20-04","SIN-20-01","SIN-20-02","SIN-20-04",
       "SIN-20-05","vic-17-01","vic-18-12")

sp<-c("DBS","DBS","DBS","DBS","DBS","DBS","DBS","DBS","DBS","DBS","DBS","DBS", 
      "DBS","DBS","DBS","DBS","DBS","DBS","DBS","DBS","DBS","DBS","DBS","DBS",
      "DBS","DBS","JH","JH","JH","JH","JH","JH","JH","LM","LM","LM","LM","LM","LM")

filter4.q2$ID[1:39]<-ind
filter4.q2$ID[40:78]<-ind
filter4.q2$sp[1:39]<-sp
filter4.q2$sp[40:78]<-sp
filter4.q2$q[1:39]<-1
filter4.q2$q[40:78]<-2

filter4.q2$sp<-factor(filter4.q2$sp,levels=c("JH","DBS","LM"))
filter4.q2<-filter4.q2[order(filter4.q2$q,filter4.q2$sp),]
filter4.q1<-filter4.q2[1:39,]
filter4.q1$sp<-factor(filter4.q1$sp,levels=c("LM","DBS","JH"))
filter4.q1<-filter4.q1[order(filter4.q1$sp,filter4.q1$mean,decreasing = TRUE),]

filter4.q12<-filter4.q2[40:78,]
filter4.q12$sp<-factor(filter4.q12$sp,levels=c("JH","DBS","LM"))
filter4.q12<-filter4.q12[order(filter4.q12$sp,filter4.q12$mean),]

filter4.q1$ID 
filter4.q12$ID

filter4.q2<-rbind(filter4.q1,filter4.q12)

filter4.q2new<-matrix(ncol=2,nrow=39)
filter4.q2new[,1]<-filter4.q2$mean[1:39]
filter4.q2new[,2]<-filter4.q2$mean[40:78]

filter4.q2newt<-t(filter4.q2new)

colors2 <-  c("red","darkblue")
par(pty = "s",mai=c(0.5,0,0.5,0))
barplot(filter4.q2newt, col = colors2, beside = FALSE,names.arg=filter4.q2$sp[1:39], las = 2, ylab = "k = 2", yaxt = "n",cex.names =0.5)
axis(2, at = c(0, 0.5, 1.0), las = 1, lwd = 2)

filter4.q2$ID<-factor(filter4.q2$ID,levels=c(filter4.q2$ID[1:39]))

filter4_ENT<-ggplot(data=filter4.q2,aes(x=ID,y=mean,fill=q))+
  geom_bar(stat="identity")+
  scale_y_continuous(expand = c(0,0),limits=c(0,1),breaks=seq(0,1,0.25))+
  theme_bw()+
  theme(
    plot.margin = unit(c(2.5,1,0.5,2.5),"cm"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 32),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text = element_text(colour = 'black',size=28),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(colour = "black",size=2),
    legend.position = "none")+
  ylab("Admixture proportion")

filter4_ENT1<-ggdraw(add_sub(filter4_ENT,c("Jackson Hole","Dubois hybrid",expression(italic("L. melissa"))),
                           x=c(0.09,0.52,0.93),y=0.4,size=16,color="white"))
filter4_ENT1

boxes <- data.frame( x = c(0.116,0.27,0.846),y = c(0,0,0))
filter4_ENT2<-ggdraw(filter4_ENT1) + 
  geom_rect(data = boxes, aes(xmin = x, xmax = x+c(0.154,0.576,0.134), alpha=0.05,
                              ymin = y, ymax = y+0.09),colour = "black", 
            fill = c("deepskyblue3","grey","skyblue1"))
filter4_ENT2

filter4_ENT3<-ggdraw(add_sub(filter4_ENT2,c("Jacksonhole","Dubois hybrid",expression(italic("L. melissa"))),
                             x=c(0.2,0.56,0.91),y=1.4,size=28))
filter4_ENT3

mergeplot1<-plot_grid(gg_inset_map1,PCAplot,nrow=1,labels = c('(A)', '(B)'),
          label_size = 42) 
mergeplot1
mergeplot2<-plot_grid(mergeplot1,filter4_ENT3,nrow=2,labels = c('', '(C)'),
                      label_size = 42,rel_heights= c(1,0.8))
mergeplot2 ### save as 20*18 ###
