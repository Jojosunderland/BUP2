rm(list=ls(all.names = T))
gc()
graphics.off()

quartz()

library(tidyverse)
library(vegan)
library(betapart)


dat.all <- read.csv("PaciFlora/Species_list_full_2905.csv",sep=";")
dat.all <- dat.all[-which(is.na(dat.all$species)),]
dat.all <- dat.all[-which(is.na(dat.all$island)),]

dat.soc <- dat.all[which(dat.all$islandgroup=="Society"),]
dat.haw <- dat.all[which(dat.all$islandgroup=="Hawaiian"),]
dat.sam <- dat.all[which(dat.all$islandgroup=="Samoa"),]
dat.mar <- dat.all[which(dat.all$islandgroup=="Marquesas"),]
dat.fij <- dat.all[which(dat.all$islandgroup=="Fiji"),]
dat.comb <- rbind(dat.soc,dat.haw,dat.sam,dat.mar,dat.fij)

length(unique(dat.soc$island))
length(unique(dat.haw$island))
length(unique(dat.sam$island))
length(unique(dat.mar$island))
length(unique(dat.fij$island))

dat.soc.red <- dat.soc[,c("species","island")]
dat.haw.red <- dat.haw[,c("species","island")]
dat.sam.red <- dat.sam[,c("species","island")]
dat.mar.red <- dat.mar[,c("species","island")]
dat.fij.red <- dat.fij[,c("species","island")]
dat.comb.red <- dat.comb[,c("species","island")]
dat.soc.red$presence <- 1
dat.haw.red$presence <- 1
dat.sam.red$presence <- 1
dat.mar.red$presence <- 1
dat.fij.red$presence <- 1
dat.comb.red$presence <- 1

##reshape - pivot matrix
dat.soc.pa <- dat.soc.red %>% 
  pivot_wider(names_from=species,values_from=c(presence))
list0 <- as.list(rep(0,ncol(dat.soc.pa)))
names(list0) <- names(dat.soc.pa)
dat.soc.pa <- as.data.frame(dat.soc.pa %>% replace_na(list0))
row.names(dat.soc.pa) <- dat.soc.pa$island
dat.soc.pa <- dat.soc.pa[,-1]

dat.haw.pa <- dat.haw.red %>% 
  pivot_wider(names_from=species,values_from=c(presence))
list0 <- as.list(rep(0,ncol(dat.haw.pa)))
names(list0) <- names(dat.haw.pa)
dat.haw.pa <- as.data.frame(dat.haw.pa %>% replace_na(list0))
row.names(dat.haw.pa) <- dat.haw.pa$island
dat.haw.pa <- dat.haw.pa[,-1]

dat.sam.pa <- dat.sam.red %>% 
  pivot_wider(names_from=species,values_from=c(presence))
list0 <- as.list(rep(0,ncol(dat.sam.pa)))
names(list0) <- names(dat.sam.pa)
dat.sam.pa <- as.data.frame(dat.sam.pa %>% replace_na(list0))
row.names(dat.sam.pa) <- dat.sam.pa$island
dat.sam.pa <- dat.sam.pa[,-1]

dat.mar.pa <- dat.mar.red %>% 
  pivot_wider(names_from=species,values_from=c(presence))
list0 <- as.list(rep(0,ncol(dat.mar.pa)))
names(list0) <- names(dat.mar.pa)
dat.mar.pa <- as.data.frame(dat.mar.pa %>% replace_na(list0))
row.names(dat.mar.pa) <- dat.mar.pa$island
dat.mar.pa <- dat.mar.pa[,-1]

dat.fij.pa <- dat.fij.red %>% 
  pivot_wider(names_from=species,values_from=c(presence))
list0 <- as.list(rep(0,ncol(dat.fij.pa)))
names(list0) <- names(dat.fij.pa)
dat.fij.pa <- as.data.frame(dat.fij.pa %>% replace_na(list0))
row.names(dat.fij.pa) <- dat.fij.pa$island
dat.fij.pa <- dat.fij.pa[,-1]

dat.comb.pa <- dat.comb.red %>% 
  pivot_wider(names_from=species,values_from=c(presence))
list0 <- as.list(rep(0,ncol(dat.comb.pa)))
names(list0) <- names(dat.comb.pa)
dat.comb.pa <- as.data.frame(dat.comb.pa %>% replace_na(list0))
row.names(dat.comb.pa) <- dat.comb.pa$island
dat.comb.pa <- dat.comb.pa[,-1]

dim(dat.soc.pa)
dim(dat.mar.pa)
dim(dat.haw.pa)
dim(dat.fij.pa)
dim(dat.sam.pa)
dim(dat.comb.pa)

###Gamma and alpha diversity
##Gamma
ncol(dat.soc.pa)
ncol(dat.haw.pa)
ncol(dat.sam.pa)
ncol(dat.mar.pa)
ncol(dat.fij.pa)
ncol(dat.comb.pa)
##Alpha
mean(rowSums(dat.soc.pa))
mean(rowSums(dat.haw.pa))
mean(rowSums(dat.sam.pa))
mean(rowSums(dat.mar.pa))
mean(rowSums(dat.fij.pa))
mean(rowSums(dat.comb.pa))


##species accumulation curves & Chao estimator
SAC.soc <- specaccum(dat.soc.pa)
SAC.haw <- specaccum(dat.haw.pa)
SAC.sam <- specaccum(dat.sam.pa)
SAC.mar <- specaccum(dat.mar.pa)
SAC.fij <- specaccum(dat.fij.pa)
SAC.comb <- specaccum(dat.comb.pa)

Estim.soc <- poolaccum(dat.soc.pa)
Estim.haw <- poolaccum(dat.haw.pa)
Estim.sam <- poolaccum(dat.sam.pa)
Estim.mar <- poolaccum(dat.mar.pa)
Estim.fij <- poolaccum(dat.fij.pa)
Estim.comb <- poolaccum(dat.comb.pa)


par(mfrow=c(2,3))
plot(SAC.soc$richness,pch=1,lty=1,lwd=2,type="b",col="blue",ylim=c(0,max(rowMeans(Estim.soc$chao))),ylab="Richness",main="Society")
points(3:nrow(dat.soc.pa),rowMeans(Estim.soc$chao),pch=2,lty=2,lwd=2,type="b",col="red")
plot(SAC.haw$richness,pch=1,lty=1,lwd=2,type="b",col="blue",ylim=c(0,max(rowMeans(Estim.haw$chao))),ylab="Richness",main="Hawai'i")
points(3:nrow(dat.haw.pa),rowMeans(Estim.haw$chao),pch=2,lty=2,lwd=2,type="b",col="red")
plot(SAC.sam$richness,pch=1,lty=1,lwd=2,type="b",col="blue",ylim=c(0,max(rowMeans(Estim.sam$chao))),ylab="Richness",main="Samoa")
points(3:nrow(dat.sam.pa),rowMeans(Estim.sam$chao),pch=2,lty=2,lwd=2,type="b",col="red")
plot(SAC.mar$richness,pch=1,lty=1,lwd=2,type="b",col="blue",ylim=c(0,max(rowMeans(Estim.mar$chao))),ylab="Richness",main="Marquesas")
points(3:nrow(dat.mar.pa),rowMeans(Estim.mar$chao),pch=2,lty=2,lwd=2,type="b",col="red")
plot(SAC.fij$richness,pch=1,lty=1,lwd=2,type="b",col="blue",ylim=c(0,max(rowMeans(Estim.fij$chao))),ylab="Richness",main="Fiji")
points(3:nrow(dat.fij.pa),rowMeans(Estim.fij$chao),pch=2,lty=2,lwd=2,type="b",col="red")
plot(SAC.comb$richness,pch=1,lty=1,lwd=2,type="b",col="blue",ylim=c(0,max(rowMeans(Estim.comb$chao))),ylab="Richness",main="All data")
points(3:nrow(dat.comb.pa),rowMeans(Estim.comb$chao),pch=2,lty=2,lwd=2,type="b",col="red")

last(rowMeans(Estim.soc$chao))/last(SAC.soc$richness)
last(rowMeans(Estim.haw$chao))/last(SAC.haw$richness)
last(rowMeans(Estim.sam$chao))/last(SAC.sam$richness)
last(rowMeans(Estim.mar$chao))/last(SAC.mar$richness)
last(rowMeans(Estim.fij$chao))/last(SAC.fij$richness)
last(rowMeans(Estim.comb$chao))/last(SAC.comb$richness)


##Beta diversity
beta.soc <- beta.pair(dat.soc.pa)
beta.haw <- beta.pair(dat.haw.pa)
beta.sam <- beta.pair(dat.sam.pa)
beta.mar <- beta.pair(dat.mar.pa)
beta.fij <- beta.pair(dat.fij.pa)
beta.comb <- beta.pair(dat.comb.pa)

mean(beta.soc$beta.sim)
mean(beta.haw$beta.sim)
mean(beta.sam$beta.sim)
mean(beta.mar$beta.sim)
mean(beta.fij$beta.sim)
mean(beta.comb$beta.sim)

mean(beta.soc$beta.sor)
mean(beta.haw$beta.sor)
mean(beta.sam$beta.sor)
mean(beta.mar$beta.sor)
mean(beta.fij$beta.sor)
mean(beta.comb$beta.sor)

coord.comb.sim <- data.frame(cmdscale(beta.comb$beta.sim))
coord.comb.sim$col <- c(rep("blue",nrow(dat.soc.pa)),rep("red",nrow(dat.haw.pa)),rep("darkgreen",nrow(dat.sam.pa)),rep("orange",nrow(dat.mar.pa)),rep("purple",nrow(dat.fij.pa)))
coord.comb.sim$pch <- c(rep(0,nrow(dat.soc.pa)),rep(1,nrow(dat.haw.pa)),rep(2,nrow(dat.sam.pa)),rep(3,nrow(dat.mar.pa)),rep(4,nrow(dat.fij.pa)))

coord.comb.sor <- data.frame(cmdscale(beta.comb$beta.sor))
coord.comb.sor$col <- c(rep("blue",nrow(dat.soc.pa)),rep("red",nrow(dat.haw.pa)),rep("darkgreen",nrow(dat.sam.pa)),rep("orange",nrow(dat.mar.pa)),rep("purple",nrow(dat.fij.pa)))
coord.comb.sor$pch <- c(rep(0,nrow(dat.soc.pa)),rep(1,nrow(dat.haw.pa)),rep(2,nrow(dat.sam.pa)),rep(3,nrow(dat.mar.pa)),rep(4,nrow(dat.fij.pa)))

par(mfrow=c(1,3))
plot(coord.comb.sim$X1,coord.comb.sim$X2,col=coord.comb.sim$col,pch=coord.comb.sim$pch,lwd=2,main="Simpson")
plot(coord.comb.sor$X1,coord.comb.sor$X2,col=coord.comb.sor$col,pch=coord.comb.sor$pch,lwd=2,main="Sorensen")
plot.new()
legend(x="topleft",legend=c("Society","Hawai'i","Samoa","Marquesas","Fiji"),col = c("blue","red","darkgreen","orange","purple"),pch=0:4,lwd=2,bty="n",lty=0,cex=2)



##OFD
freq.soc <- colSums(dat.soc.pa)/nrow(dat.soc.pa)
freq.haw <- colSums(dat.haw.pa)/nrow(dat.haw.pa)
freq.sam <- colSums(dat.sam.pa)/nrow(dat.sam.pa)
freq.mar <- colSums(dat.mar.pa)/nrow(dat.mar.pa)
freq.fij <- colSums(dat.fij.pa)/nrow(dat.fij.pa)
freq.comb <- colSums(dat.comb.pa)/nrow(dat.comb.pa)

par(mfrow=c(2,3))
hist(freq.soc,main="Society",xlab="Occupancy",breaks = seq(0,1,0.1))
hist(freq.haw,main="Hawa'i",xlab="Occupancy",breaks = seq(0,1,0.1))
hist(freq.sam,main="Samoa",xlab="Occupancy",breaks = seq(0,1,0.1))
hist(freq.mar,main="Marquesas",xlab="Occupancy",breaks = seq(0,1,0.1))
hist(freq.fij,main="Fiji",xlab="Occupancy",breaks = seq(0,1,0.1))
hist(freq.comb,main="Combined",xlab="Occupancy",breaks = seq(0,1,0.1))
























