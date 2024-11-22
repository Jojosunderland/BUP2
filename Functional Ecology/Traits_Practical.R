library(tidyverse)
library(betapart)
library(BAT)
library(ggplot2)
library(ggrepel)

dat.all <- read.csv("PaciFlora/Species_list_full_2905.csv",sep=";")
dat.all <- dat.all[-which(is.na(dat.all$species)),]
dat.all <- dat.all[-which(is.na(dat.all$island)),]

dat.soc <- dat.all[which(dat.all$islandgroup=="Society"),]
length(unique(dat.soc$island))

dat.soc.red <- dat.soc[,c("species","island")]
dat.soc.red$presence <- 1

##reshape - pivot matrix
dat.soc.pa <- dat.soc.red %>% 
  pivot_wider(names_from=species,values_from=c(presence))
list0 <- as.list(rep(0,ncol(dat.soc.pa)))
names(list0) <- names(dat.soc.pa)
dat.soc.pa <- as.data.frame(dat.soc.pa %>% replace_na(list0))
row.names(dat.soc.pa) <- dat.soc.pa$island
dat.soc.pa <- dat.soc.pa[,-1] 
dat.soc.pa <- dat.soc.pa[order(row.names(dat.soc.pa)),] ##ordering rows alphabetically

##Gamma
gamma.soc <- ncol(dat.soc.pa)
##Alpha
alpha.soc <- rowSums(dat.soc.pa)
alpha.soc.mean <- mean(alpha.soc)

##The ratio of alpha diversity compared to gamma diversity
alpha.soc.ratio <- alpha.soc/gamma.soc
alpha.soc.mean.ratio <- alpha.soc.mean/gamma.soc

# Load the betapart and the BAT packages, and use the functions to compute the beta diversity indices for the two partitions 

# Beta
beta.soc1 <- beta.pair(dat.soc.pa)
beta.soc2 <- beta(dat.soc.pa,func = "sorensen")

plot(beta.soc1$beta.sor,beta.soc2$Btotal,xlim=c(0,1),ylim=c(0,1))
lines(c(0,1),c(0,1),col="red")
plot(beta.soc1$beta.sim,beta.soc2$Brepl,xlim=c(0,1),ylim=c(0,1))
lines(c(0,1),c(0,1),col="red")

## We will use the function beta() from the BAT R package, as it allows to compute species, trait and phylogenetic diversity

## Trait Diversity ##
# To compute trait diversity, we first need to get trait data for the species in the community
# we will use seed mass and plant height, which are amongst the most commonly available traits


## Here we will use the TR8 R package, to extract this information from multiple plant databases. 
## You can see all databases and the traits they contain included in the package by typing:

install.packages('TR8')
library(TR8)
available_traits()

names.spp.soc <- unique(dat.soc$species)
Soc.sp.tr8 <- tr8(species_list = names.spp.soc, download_list=c("seed_mas_cal","seed_mass","seed_wght","SeedMass","Height","h_max","max_height_cal"),allow_persistent=FALSE) 
View(Soc.sp.tr8@results)

library(stringi)
traits <- Soc.sp.tr8@results
traits$h_max <- as.numeric(stri_extract_first_regex(traits$h_max, "[0-9]+"))
traits$seed_wght <- as.numeric(stri_extract_first_regex(traits$seed_wght, "[0-9]+"))
traits$Height <- as.numeric(stri_extract_first_regex(traits$Height, "[0-9]+"))
traits$SeedMass <- as.numeric(stri_extract_first_regex(traits$SeedMass, "[0-9]+"))
View(traits)

traits.mean <- data.frame(height=numeric(nrow(traits)),seed.mass=numeric(nrow(traits)))
row.names(traits.mean) <- row.names(traits)
traits.mean$height <- apply(traits[,c(1,4)],1,mean,na.rm=TRUE)
traits.mean$seed.mass <- apply(traits[,c(2,3,5)],1,mean,na.rm=TRUE)
View(traits.mean)

traits.mean <- traits.mean[-which(is.na(traits.mean$height) | is.na(traits.mean$seed.mass)),]

