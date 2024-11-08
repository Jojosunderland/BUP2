## Today, we will compute some of the community patterns seen in class for data 
## on alien plant species in Pacific islands, to explore how alien plant communities are structured.

# install.packages
install.packages("tidyverse")
install.packages("vegan")
install.packages("betapart")
library(tidyverse)
library(vegan)
library(betapart)

## we will use the paciflora database
#load the databases

dat.all <- read.csv("Species_list_full_2905.csv",sep=";")
View(dat.all)
## this data frame is in a “tidy” format, i.e. each row represents a unique combination of site and species, with additional information about these two elements. 

# remove NAs for species or islands

dat.all <- dat.all[-which(is.na(dat.all$species)),] #which species have NAs
dat.all <- dat.all[-which(is.na(dat.all$island)),] #which islands have NAs

# create 6 data frames
dat.soc <- dat.all[which(dat.all$islandgroup=="Society"),] # data for islands belong to the Society archipelago
dat.haw <- dat.all[which(dat.all$islandgroup=="Hawaiian"),] # Hawaiian archipelago
dat.sam <- dat.all[which(dat.all$islandgroup=="Samoa"),] # Samoan archipelago
dat.mar <- dat.all[which(dat.all$islandgroup=="Marquesas"),] # Marquesas archipelago
dat.fij <- dat.all[which(dat.all$islandgroup=="Fiji"),] # Fiji archepelago

# combine all the dataframes above
dat.comb <- rbind(dat.soc,dat.haw,dat.sam,dat.mar,dat.fij)

# this tells you how many unique/different islands there are in each area
length(unique(dat.soc$island))
length(unique(dat.haw$island))
length(unique(dat.sam$island))
length(unique(dat.mar$island))
length(unique(dat.fij$island))

## We need site by species data frames to compute community patterns

## create a site by species presence absence matrix for each data frame

# Society archepelago
dat.soc.red <- dat.soc[,c("species","island")] 
dat.soc.red$presence <- 1  # first two lines create a data frame with 3 columns, islands, species and presences (set to 1s)
dat.soc.pa <- dat.soc.red %>% 
  pivot_wider(names_from=species,values_from=c(presence)) # pivot_wider creates a wide-format data frame, each island becomes a row, each species becomes a column
# the values inside the matrix indicate presence (1) or absence (0)
list0 <- as.list(rep(0,ncol(dat.soc.pa))) # creates a list where each column name is assigned a value of 0
names(list0) <- names(dat.soc.pa)
dat.soc.pa <- as.data.frame(dat.soc.pa %>% replace_na(list0)) # replace NA values with 0, representing species absence
row.names(dat.soc.pa) <- dat.soc.pa$island #island columns are assigned as row names, columns are species
dat.soc.pa <- dat.soc.pa[,-1] # delete the first column, island column as it is now redundant after being set as row names

View(dat.soc.pa)
## dat.soc.pa is the presence absence matrix ##

# Marquesas archipelago
dat.mar.red <- dat.mar[,c("species","island")] 
dat.mar.red$presence <- 1 
dat.mar.pa <- dat.mar.red %>% 
  pivot_wider(names_from=species,values_from=c(presence))
list0 <- as.list(rep(0,ncol(dat.mar.pa)))
names(list0) <- names(dat.mar.pa)
dat.mar.pa <- as.data.frame(dat.mar.pa %>% replace_na(list0))
row.names(dat.mar.pa) <- dat.mar.pa$island
dat.mar.pa <- dat.mar.pa[,-1]

View(dat.mar.pa)

# Hawaiian archipelago
dat.haw.red <- dat.haw[,c("species","island")] 
dat.haw.red$presence <- 1 
dat.haw.pa <- dat.haw.red %>% 
  pivot_wider(names_from=species,values_from=c(presence))
list0 <- as.list(rep(0,ncol(dat.haw.pa)))
names(list0) <- names(dat.haw.pa)
dat.haw.pa <- as.data.frame(dat.haw.pa %>% replace_na(list0))
row.names(dat.haw.pa) <- dat.haw.pa$island
dat.haw.pa <- dat.haw.pa[,-1]

View(dat.haw.pa)

# Samoan archipelago
dat.sam.red <- dat.sam[,c("species","island")] 
dat.sam.red$presence <- 1 
dat.sam.pa <- dat.sam.red %>% 
  pivot_wider(names_from=species,values_from=c(presence))
list0 <- as.list(rep(0,ncol(dat.sam.pa)))
names(list0) <- names(dat.sam.pa)
dat.sam.pa <- as.data.frame(dat.sam.pa %>% replace_na(list0))
row.names(dat.sam.pa) <- dat.sam.pa$island
dat.sam.pa <- dat.sam.pa[,-1]

View(dat.sam.pa)

# Fiji archipelago
dat.fij.red <- dat.fij[,c("species","island")] 
dat.fij.red$presence <- 1 
dat.fij.pa <- dat.fij.red %>% 
  pivot_wider(names_from=species,values_from=c(presence))
list0 <- as.list(rep(0,ncol(dat.fij.pa)))
names(list0) <- names(dat.fij.pa)
dat.fij.pa <- as.data.frame(dat.fij.pa %>% replace_na(list0))
row.names(dat.fij.pa) <- dat.fij.pa$island
dat.fij.pa <- dat.fij.pa[,-1]

View(dat.fij.pa)

# Combined archipelago
dat.comb.red <- dat.comb[,c("species","island")] 
dat.comb.red$presence <- 1 
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

## Calculate gamma and alpha diversity
##Gamma - count number of columns across all islands
ncol(dat.soc.pa)
ncol(dat.haw.pa)
ncol(dat.sam.pa)
ncol(dat.mar.pa)
ncol(dat.fij.pa)
ncol(dat.comb.pa)

##Alpha - mean number of species at each island added together
mean(rowSums(dat.soc.pa))
mean(rowSums(dat.haw.pa))
mean(rowSums(dat.sam.pa))
mean(rowSums(dat.mar.pa))
mean(rowSums(dat.fij.pa))
mean(rowSums(dat.comb.pa))

## Species accumulation curves
?specaccum() 

SAC.soc <- specaccum(dat.soc.pa)
SAC.haw <- specaccum(dat.haw.pa)
SAC.sam <- specaccum(dat.sam.pa)
SAC.mar <- specaccum(dat.mar.pa)
SAC.fij <- specaccum(dat.fij.pa)
SAC.comb <- specaccum(dat.comb.pa)

plot(SAC.soc)
plot(SAC.haw)
plot(SAC.sam)
plot(SAC.mar)
plot(SAC.fij)
plot(SAC.comb)

## Chao2 estimator
?poolaccum() 
Estim.soc <- poolaccum(dat.soc.pa)
Estim.haw <- poolaccum(dat.haw.pa)
Estim.sam <- poolaccum(dat.sam.pa)
Estim.mar <- poolaccum(dat.mar.pa)
Estim.fij <- poolaccum(dat.fij.pa)
Estim.comb <- poolaccum(dat.comb.pa)

# plot the SAC with the chao2 estimator
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

## Beta diversity - using the betapart package.
?beta.pair() 

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

## 2D plot based on the dissimilarity between them, using the combined data.

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

quartz()

## Occupancy frequency distribution (OFD)
# compute the occupancy of each species relative to the total number of islands. 
# For each archipelago, you should obtain one vector with one value in [0,1] for each island.

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



