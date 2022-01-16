################################
library(plyr)
library(cowplot)
####################################################################
#### make data file document alpha and beta "bg_filter4.sum" #######
####################################################################
#### combine different chains of alpha and beta files ########
bg_alpha4_1<-read.table("bgcout_pct3_filter4_1_alpha.txt",sep=",")
bg_alpha4_2<-read.table("bgcout_pct3_filter4_2_alpha.txt",sep=",")
bg_alpha4_3<-read.table("bgcout_pct3_filter4_3_alpha.txt",sep=",")
bg_alpha4_4<-read.table("bgcout_pct3_filter4_4_alpha.txt",sep=",")

bg_beta4_1<-read.table("bgcout_pct3_filter4_1_beta.txt",sep=",")
bg_beta4_2<-read.table("bgcout_pct3_filter4_2_beta.txt",sep=",")
bg_beta4_3<-read.table("bgcout_pct3_filter4_3_beta.txt",sep=",")
bg_beta4_4<-read.table("bgcout_pct3_filter4_4_beta.txt",sep=",")

bg_alpha4<-cbind(bg_alpha4_1[,-1],bg_alpha4_2[,-1],bg_alpha4_3[,-1],bg_alpha4_4[,-1])
bg_beta4<-cbind(bg_beta4_1[,-1],bg_beta4_2[,-1],bg_beta4_3[,-1],bg_beta4_4[,-1])

bg_filter4.sum<-data.frame(alpha.mean=rowMeans(bg_alpha4))
bg_filter4.sum$alpha.LB<-apply(bg_alpha4,1,quantile,probs=0.025)
bg_filter4.sum$alpha.UB<-apply(bg_alpha4,1,quantile,probs=0.975)
bg_filter4.sum$beta.mean<-rowMeans(bg_beta4)
bg_filter4.sum$beta.LB<-apply(bg_beta4,1,quantile,probs=0.025)
bg_filter4.sum$beta.UB<-apply(bg_beta4,1,quantile,probs=0.975)
bg_filter4.sum$alpha.median<-apply(bg_alpha4,1,median)
bg_filter4.sum$beta.median<-apply(bg_beta4,1,median)

##### filter 4 position ######
sv.LKnew4<-read.table("sv.LKnew4.txt") 
SV.LKnew4<-SV.LKnew4[,-1]
length(SV.LKnew4[,1])
bg_filter4.sum$type<-SV.LKnew4$type

bg_filter4.sum$position<-NA
bg_filter4.sum$position<-SV.LKnew4$position

bg_filter4.sum$scaffold<-NA
bg_filter4.sum$start<-NA

for (i in 1:1419){
  bg_filter4.sum$scaffold[i]<-strsplit(bg_filter4.sum$position[i], "_")[[1]][2]
  bg_filter4.sum$start[i]<-strsplit(bg_filter4.sum$position[i], "_")[[1]][5]
}

head(bg_filter4.sum)

bg_filter4.sum$length<-SV.LKnew4$length
bg_filter4.sum$start<-as.numeric(bg_filter4.sum$start)
bg_filter4.sum$length<-as.numeric(bg_filter4.sum$length)
bg_filter4.sum$end<-bg_filter4.sum$start+bg_filter4.sum$length

head(bg_filter4.sum)
mean(bg_filter4.sum$alpha.mean)
sd(bg_filter4.sum$alpha.mean)/sqrt(1419)
###   write.table(bg_filter4.sum,"bg_filter4.sum.txt")
################# multipanel plot ###############################################################

par(mfrow=c(1,2))
### Figure 5A hybrid index plot ####
h<-read.table("cline7b_hi_mean.txt", header=T)
#par(mar=c(7,7,3,0.3))
par(mar=c(5.5,6,5,1))

par(mai = c(1, 1.1, 1, 0.1))

hist(h[,2],xlim=c(0,1),col="lavender",xlab="Hybrid index",ylab="Frequency",main="",cex.lab=2.5, cex.axis=2)
title(main="(A)",cex.main = 4, adj=0)
### Figure 5B cline plot ####
## function for cline plots
##getting hi and low alpha and beta
hiA<-which(bg_filter4.sum$alpha.LB > 0) #number of loci have significant positive alpha value: 152
lowA<-which(bg_filter4.sum$alpha.UB < 0) #number of loci have significant negative alpha value: 410
hiB<-which(bg_filter4.sum$beta.LB > 0) #number of loci have significant positive beta value: 0
lowB<-which(bg_filter4.sum$beta.UB < 0) #number of loci have significant negative beta value: 0
length(bg_filter4.sum[,1])
length(bg_filter4.sum[bg_filter4.sum$alpha.mean<0,1])

#neutral
n<-length(bg_filter4.sum[,1])
xx<-1:n ###
notneu<-which(xx %in% c(lowA,lowB,hiA,hiB))
neu<-xx[-notneu]
subNeu<-sample(neu,500,replace=FALSE)

comlowA<-which(lowA %in% c(lowB,hiB))
comlowB<-which(lowB %in% c(lowB,hiB))
comhiA<-which(hiA %in% c(lowB,hiB))
##only these loci
ohiA<-which(bg_filter4.sum$alpha.LB> 0 & bg_filter4.sum$beta.LB <= 0 & bg_filter4.sum$beta.UB >= 0) #152
olowA<-which(bg_filter4.sum$alpha.UB < 0 & bg_filter4.sum$beta.LB <= 0 & bg_filter4.sum$beta.UB >= 0) #410
ohiB<-which(bg_filter4.sum$alpha.LB> 0 & bg_filter4.sum$alpha.LB  <= 0 & bg_filter4.sum$alpha.UB >= 0) #0
#####################
calcCline<-function(h=NULL, a=NULL, b=NULL){
  z <- 2 * h * (1-h) * (a + (2 * h - 1) * b)
  p <- h + z
  p[p<0]<-0
  p[p>1]<-1
  return(p)
}
#genomic clines plot
#par(mar=c(7,7,3,0.3))
h<-seq(0,1,0.01)
phi<-calcCline(h,bg_filter4.sum$alpha.mean[1],bg_filter4.sum$beta.mean[1])
plot(h,phi,type='l',xlab="Hybrid index",ylab="Probs. of LM anc.",col="darkgray",lwd=.7,cex.lab=2.5, cex.axis=2)
for(i in subNeu){
  phi<-calcCline(h,bg_filter4.sum$alpha.mean[i],bg_filter4.sum$beta.mean[i])
  lines(h,phi,col="darkgray",lwd=2)
}
for(i in ohiA){
  phi<-calcCline(h,bg_filter4.sum$alpha.mean[i],bg_filter4.sum$beta.mean[i])
  lb<-which(phi == 1)[1]
  ub<-length(phi)
  phi[c(lb:ub)]<-1
  lines(h,phi,col="darkgreen",lwd=2)
}
for(i in olowA){
  phi<-calcCline(h,bg_filter4.sum$alpha.mean[i],bg_filter4.sum$beta.mean[i])
  lines(h,phi,col="darkolivegreen3",lwd=2)
}
for(i in ohiB){
  phi<-calcCline(h,bg_filter4.sum$alpha.mean[i],bg_filter4.sum$beta.mean[i])
  lines(h,phi,col="darkorchid",lwd=3)
}
abline(a=0,b=1,lty=2,lwd=2)
title(main="(B)",cex.main = 4, adj=0)
p1 <- recordPlot()  
p1
##################################################################################
#### assign scaffold to different chromosome #####
lgs<-read.table("lgs.txt",header=TRUE)
colnames(lgs)[1]<-"scaffold"
### filter 4  ####
bg_filter4.sumnew2<-merge(bg_filter4.sum,lgs,by="scaffold",all.x=TRUE)
length(bg_filter4.sumnew2[!is.na(bg_filter4.sumnew2$LG),1])
bg_filter4.sumnew2<-bg_filter4.sumnew2[!is.na(bg_filter4.sumnew2$LG),]
summary(as.factor(bg_filter4.sumnew2$LG))
str(bg_filter4.sumnew2)
length(bg_filter4.sumnew2[bg_filter4.sumnew2$high.CI<0,1])

bg_filter4.sumnew2$LG<-as.factor(bg_filter4.sumnew2$LG)
summary(as.factor(bg_filter4.sumnew2$LG))
bg_filter4.sumnew2$LG<-factor(bg_filter4.sumnew2$LG,levels=c("1","2","3","4","5","6","7","8","9","10","11","12",
                                                        "13","14","15","16","17","18","19","20","21","22","Z"))


LGavg<-ddply(bg_filter4.sumnew2,"LG",summarize,meanalpha=mean(alpha.mean),meanbeta=mean(beta.mean))
bg_filter4.sumnew2$significant<-"N"
bg_filter4.sumnew2$significant[bg_filter4.sumnew2$alpha.LB>0 | bg_filter4.sumnew2$alpha.UB<0]<-"Y"
bg_filter4.sumnew2$significant<-factor(bg_filter4.sumnew2$significant,levels=c("Y","N"))
alpha4.plot<-ggplot()+
  geom_point(data=bg_filter4.sumnew2,aes(x=LG,y=alpha.mean,fill=significant,shape=significant,color=significant),size=4) +
  scale_shape_manual(values=c(21,21))+
  scale_fill_manual(values=c("mediumpurple3","white"))+
  scale_color_manual(values=c("mediumpurple3","mediumpurple3"))+
  scale_y_continuous(limits=c(-2,2),breaks=seq(-2,2,1))+
  geom_segment(data=LGavg,aes(x=seq(0.75,23.25,1),xend=seq(1.25,23.25,1),y=meanalpha,yend=meanalpha),color="purple3",size=4)+
  geom_hline(yintercept=0, linetype="dashed",color = "black", size=1)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background=element_blank(),
        plot.margin=unit(c(1,1,1,1),"cm"),
        text=element_text(colour="black",size=32),
        axis.text=element_text(colour="black",size=34),
        strip.text = element_text(face = "italic",size=32),
        legend.position = "none",
        axis.title.y=element_text(face="italic",size=48),
        axis.title.x=element_blank())+
  ylab(expression(alpha))
alpha4.plot

beta4.plot<-ggplot()+
  geom_point(data=bg_filter4.sumnew2,aes(x=LG,y=beta.mean,fill=significant,shape=significant,color=significant),size=4) +
  scale_shape_manual(values=c(21,21))+
  scale_fill_manual(values=c("green4","white"))+
  scale_color_manual(values=c("green4","green4"))+
  scale_y_continuous(limits=c(-2,2),breaks=seq(-2,2,1))+
  geom_segment(data=LGavg,aes(x=seq(0.75,23.25,1),xend=seq(1.25,23.25,1),y=meanbeta,yend=meanbeta),color="darkseagreen4",size=4)+
  geom_hline(yintercept=0, linetype="dashed",color = "black", size=1)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background=element_blank(),
        plot.margin=unit(c(1,1,1,1),"cm"),
        text=element_text(colour="black",size=32),
        axis.text=element_text(colour="black",size=34),
        strip.text = element_text(face = "italic",size=32),
        legend.position = "none",
        axis.title.y=element_text(face="italic",size=48),
        axis.title.x=element_blank())+
  ylab(expression(alpha))
beta4.plot
#### combine plots for figure 5 #####
plot_grid(p1,plot_grid(alpha4.plot, beta4.plot,labels = c('(C)', '(D)'),
                       align = "v", ncol = 1, rel_heights = c(0.5, 0.5),
                       label_size = 42),nrow=2,rel_heights=c(0.35,0.65)) 

head(bg_filter4.sumnew2)
####################################################################################
############## cline parameter analysis ######################
### comparing alpha of autosome with alpha of Z chromosome ###
bg_filter4.sumnew3<-bg_filter4.sumnew2[bg_filter4.sumnew2$LG!="Z",]
length(bg_filter4.sumnew3[,1])
length(bg_filter4.sumnew2[bg_filter4.sumnew2$LG=="Z",1])
alpha.auto<-c()
for(i in 1:1000){
  temp<-bg_filter4.sumnew3[sample(1:1131,236),"alpha.mean"]
  alpha.auto[i]<-mean(temp)
}
alpha.Z<-mean(bg_filter4.sumnew2[bg_filter4.sumnew2$LG=="Z","alpha.mean"])
mean(alpha.auto-alpha.Z)

length(alpha.auto[alpha.auto<alpha.Z])## number of alpha.auto smaller than alpha.Z
############################################################################

