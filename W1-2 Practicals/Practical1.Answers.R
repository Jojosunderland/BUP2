rm(list=ls(all.names = TRUE))
gc()

library(dplyr)

##Import csv files
landuse <- read.csv("landuse/assessments.csv")
overexploitation <- read.csv("overexploitation/assessments.csv")
invasion <- read.csv("invasion/assessments.csv")
pollution <- read.csv("pollution/assessments.csv")
climate <- read.csv("climate/assessments.csv")

unlist(lapply(landuse, class))
unlist(lapply(overexploitation, class))
unlist(lapply(invasion, class))
unlist(lapply(pollution, class))
unlist(lapply(climate, class))

str(landuse)
str(overexploitation)
str(invasion)
str(pollution)
str(climate)

landuse$landuse <- 1
landuse$overexploitation <- 0
landuse$invasion <- 0
landuse$pollution <- 0
landuse$climate <- 0

overexploitation$landuse <- 0
overexploitation$overexploitation <- 1
overexploitation$invasion <- 0
overexploitation$pollution <- 0
overexploitation$climate <- 0

invasion$landuse <- 0
invasion$overexploitation <- 0
invasion$invasion <- 1
invasion$pollution <- 0
invasion$climate <- 0

pollution$landuse <- 0
pollution$overexploitation <- 0
pollution$invasion <- 0
pollution$pollution <- 1
pollution$climate <- 0

climate$landuse <- 0
climate$overexploitation <- 0
climate$invasion <- 0
climate$pollution <- 0
climate$climate <- 1

threats.dat <- rbind(landuse, overexploitation, invasion, pollution, climate)
dim(threats.dat)
length(unique(threats.dat$scientificName))

threats.dup <- which(duplicated(threats.dat$scientificName))

threats.redundant <- distinct(threats.dat[threats.dup,1:23])

for(i in 1:nrow(threats.redundant)){
  row.number <- which(threats.dat$scientificName==threats.redundant$scientificName[i])
  row.selected <- threats.dat[row.number,]
  new.row <- data.frame(c(row.selected[1,1:23],colSums(row.selected[,24:28])))
  threats.dat <- threats.dat[-row.number,]
  threats.dat <- rbind(threats.dat,new.row)
}
dim(threats.dat)


##General patterns
#if (!require(devtools)) install.packages("devtools")
#if (!require(ggVennDiagram)) devtools::install_github("gaospecial/ggVennDiagram")
library("ggVennDiagram")

colSums(threats.dat[,24:28])
barplot(colSums(threats.dat[,24:28]))

#Venn diagram
x <- list(
  landuse=landuse$scientificName,
  overexploitation=overexploitation$scientificName,
  invasion=invasion$scientificName,
  pollution=pollution$scientificName,
  climate=climate$scientificName
)
ggVennDiagram(x, label_alpha = 0, set_color = c("blue","black","red","purple","green"))

##another way to do it
y <- list()
for(i in 1:5){
  y[[i]] <- threats.dat$scientificName[which(threats.dat[23+i]==1)]
}
names(y) <- names(threats.dat[24:28])
ggVennDiagram(y, label_alpha = 0, set_color = c("blue","black","red","purple","green"))


##UpSet plot
if (!require(UpSetR)) install.packages("UpSetR")
library(UpSetR)

upset(data=threats.dat[c(3,24:28)],nsets=5)



#######################
##Per taxonomic group##
landuse.tax <- read.csv("landuse/taxonomy.csv")
overexploitation.tax <- read.csv("overexploitation/taxonomy.csv")
invasion.tax <- read.csv("invasion/taxonomy.csv")
pollution.tax <- read.csv("pollution/taxonomy.csv")
climate.tax <- read.csv("climate/taxonomy.csv")

all.tax <- rbind(landuse.tax,overexploitation.tax,invasion.tax,pollution.tax,climate.tax)
all.tax <- distinct(all.tax)
dim(all.tax)

all.tax <- all.tax[order(all.tax$scientificName),]
threats.dat <- threats.dat[order(threats.dat$scientificName),]

table(all.tax$className)

threats.dat.mammals <- threats.dat[which(all.tax$className=="MAMMALIA"),]
threats.dat.birds <- threats.dat[which(all.tax$className=="AVES"),]
threats.dat.fish <- threats.dat[which(all.tax$className=="ACTINOPTERYGII"),]
threats.dat.amphibian <- threats.dat[which(all.tax$className=="AMPHIBIA"),]
threats.dat.gastropode <- threats.dat[which(all.tax$className=="GASTROPODA"),]
threats.dat.flower <- threats.dat[which(all.tax$className=="MAGNOLIOPSIDA"),]
threats.dat.reptile <- threats.dat[which(all.tax$className=="REPTILIA"),]


upset(data=threats.dat[c(3,24:28)],nsets=5)
upset(data=threats.dat.mammals[c(3,24:28)],nsets=5)
upset(data=threats.dat.birds[c(3,24:28)],nsets=5)
upset(data=threats.dat.fish[c(3,24:28)],nsets=5)
upset(data=threats.dat.amphibian[c(3,24:28)],nsets=5)
upset(data=threats.dat.gastropode[c(3,24:28)],nsets=5)
upset(data=threats.dat.flower[c(3,24:28)],nsets=5)
upset(data=threats.dat.reptile[c(3,24:28)],nsets=5)


###################
##temporal trends##
range(threats.dat$yearPublished)

threats.dat.temp <- threats.dat[which(threats.dat$yearPublished <= 1996),]
threats.dat.temp.mammal <- threats.dat.mammals[which(threats.dat.mammals$yearPublished <= 1996),]
threats.dat.temp.birds <- threats.dat.birds[which(threats.dat.birds$yearPublished <= 1996),]
threats.dat.temp.fish <- threats.dat.fish[which(threats.dat.fish$yearPublished <= 1996),]
threats.dat.temp.amphibian <- threats.dat.amphibian[which(threats.dat.amphibian$yearPublished <= 1996),]
threats.dat.temp.gastropode <- threats.dat.gastropode[which(threats.dat.gastropode$yearPublished <= 1996),]
threats.dat.temp.flower <- threats.dat.flower[which(threats.dat.flower$yearPublished <= 1996),]
threats.dat.temp.reptile <- threats.dat.reptile[which(threats.dat.reptile$yearPublished <= 1996),]
threats.time <- colSums(threats.dat.temp[,24:28])
threats.time.mammal <- colSums(threats.dat.temp.mammal[,24:28])
threats.time.birds <- colSums(threats.dat.temp.birds[,24:28])
threats.time.fish <- colSums(threats.dat.temp.fish[,24:28])
threats.time.amphibian <- colSums(threats.dat.temp.amphibian[,24:28])
threats.time.gastropode <- colSums(threats.dat.temp.gastropode[,24:28])
threats.time.flower <- colSums(threats.dat.temp.flower[,24:28])
threats.time.reptile <- colSums(threats.dat.temp.reptile[,24:28])
for(y in 1997:2024){
  threats.dat.temp <- threats.dat[which(threats.dat$yearPublished <= y),]
  threats.dat.temp.mammal <- threats.dat.mammals[which(threats.dat.mammals$yearPublished <= y),]
  threats.dat.temp.birds <- threats.dat.birds[which(threats.dat.birds$yearPublished <= y),]
  threats.dat.temp.fish <- threats.dat.fish[which(threats.dat.fish$yearPublished <= y),]
  threats.dat.temp.amphibian <- threats.dat.amphibian[which(threats.dat.amphibian$yearPublished <= y),]
  threats.dat.temp.gastropode <- threats.dat.gastropode[which(threats.dat.gastropode$yearPublished <= y),]
  threats.dat.temp.flower <- threats.dat.flower[which(threats.dat.flower$yearPublished <= y),]
  threats.dat.temp.reptile <- threats.dat.reptile[which(threats.dat.reptile$yearPublished <= y),]
  threats.time <- rbind(threats.time,colSums(threats.dat.temp[,24:28]))
  threats.time.mammal <- rbind(threats.time.mammal,colSums(threats.dat.temp.mammal[,24:28]))
  threats.time.birds <- rbind(threats.time.birds,colSums(threats.dat.temp.birds[,24:28]))
  threats.time.fish <- rbind(threats.time.fish,colSums(threats.dat.temp.fish[,24:28]))
  threats.time.amphibian <- rbind(threats.time.amphibian,colSums(threats.dat.temp.amphibian[,24:28]))
  threats.time.gastropode <- rbind(threats.time.gastropode,colSums(threats.dat.temp.gastropode[,24:28]))
  threats.time.flower <- rbind(threats.time.flower,colSums(threats.dat.temp.flower[,24:28]))
  threats.time.reptile <- rbind(threats.time.reptile,colSums(threats.dat.temp.reptile[,24:28]))
}

quartz()
par(mfrow=c(2,4))
plot(x=1996:2024,y=threats.time[,1],type = "b",pch=15,ylim=c(1,max(threats.time)),xlab="Year",ylab="Number of species",main="All species",log="y")
lines(x=1996:2024,y=threats.time[,2],type = "b",pch=16,col="gray10")
lines(x=1996:2024,y=threats.time[,3],type = "b",pch=17,col="gray20")
lines(x=1996:2024,y=threats.time[,4],type = "b",pch=18,col="gray30")
lines(x=1996:2024,y=threats.time[,5],type = "b",pch=19,col="gray40")
legend(x="topleft",legend=c("landuse","overexploitation","invasion","pollution","climate"),bty="n",lty=1,pch=15:19,col = c("black","gray10","gray20","gray30","gray40"))

plot(x=1996:2024,y=threats.time.mammal[,1],type = "b",pch=15,ylim=c(1,max(threats.time)),xlab="Year",ylab="Number of species",main="Mammals",log="y")
lines(x=1996:2024,y=threats.time.mammal[,2],type = "b",pch=16,col="gray10")
lines(x=1996:2024,y=threats.time.mammal[,3],type = "b",pch=17,col="gray20")
lines(x=1996:2024,y=threats.time.mammal[,4],type = "b",pch=18,col="gray30")
lines(x=1996:2024,y=threats.time.mammal[,5],type = "b",pch=19,col="gray40")

plot(x=1996:2024,y=threats.time.birds[,1],type = "b",pch=15,ylim=c(1,max(threats.time)),xlab="Year",ylab="Number of species",main="Birds",log="y")
lines(x=1996:2024,y=threats.time.birds[,2],type = "b",pch=16,col="gray10")
lines(x=1996:2024,y=threats.time.birds[,3],type = "b",pch=17,col="gray20")
lines(x=1996:2024,y=threats.time.birds[,4],type = "b",pch=18,col="gray30")
lines(x=1996:2024,y=threats.time.birds[,5],type = "b",pch=19,col="gray40")

plot(x=1996:2024,y=threats.time.fish[,1],type = "b",pch=15,ylim=c(1,max(threats.time)),xlab="Year",ylab="Number of species",main="Fish",log="y")
lines(x=1996:2024,y=threats.time.fish[,2],type = "b",pch=16,col="gray10")
lines(x=1996:2024,y=threats.time.fish[,3],type = "b",pch=17,col="gray20")
lines(x=1996:2024,y=threats.time.fish[,4],type = "b",pch=18,col="gray30")
lines(x=1996:2024,y=threats.time.fish[,5],type = "b",pch=19,col="gray40")

plot(x=1996:2024,y=threats.time.amphibian[,1],type = "b",pch=15,ylim=c(1,max(threats.time)),xlab="Year",ylab="Number of species",main="Amphibians",log="y")
lines(x=1996:2024,y=threats.time.amphibian[,2],type = "b",pch=16,col="gray10")
lines(x=1996:2024,y=threats.time.amphibian[,3],type = "b",pch=17,col="gray20")
lines(x=1996:2024,y=threats.time.amphibian[,4],type = "b",pch=18,col="gray30")
lines(x=1996:2024,y=threats.time.amphibian[,5],type = "b",pch=19,col="gray40")

plot(x=1996:2024,y=threats.time.gastropode[,1],type = "b",pch=15,ylim=c(1,max(threats.time)),xlab="Year",ylab="Number of species",main="Gastropodes",log="y")
lines(x=1996:2024,y=threats.time.gastropode[,2],type = "b",pch=16,col="gray10")
lines(x=1996:2024,y=threats.time.gastropode[,3],type = "b",pch=17,col="gray20")
lines(x=1996:2024,y=threats.time.gastropode[,4],type = "b",pch=18,col="gray30")
lines(x=1996:2024,y=threats.time.gastropode[,5],type = "b",pch=19,col="gray40")

plot(x=1996:2024,y=threats.time.flower[,1],type = "b",pch=15,ylim=c(1,max(threats.time)),xlab="Year",ylab="Number of species",main="Flowers",log="y")
lines(x=1996:2024,y=threats.time.flower[,2],type = "b",pch=16,col="gray10")
lines(x=1996:2024,y=threats.time.flower[,3],type = "b",pch=17,col="gray20")
lines(x=1996:2024,y=threats.time.flower[,4],type = "b",pch=18,col="gray30")
lines(x=1996:2024,y=threats.time.flower[,5],type = "b",pch=19,col="gray40")

plot(x=1996:2024,y=threats.time.reptile[,1],type = "b",pch=15,ylim=c(1,max(threats.time)),xlab="Year",ylab="Number of species",main="Reptiles",log="y")
lines(x=1996:2024,y=threats.time.reptile[,2],type = "b",pch=16,col="gray10")
lines(x=1996:2024,y=threats.time.reptile[,3],type = "b",pch=17,col="gray20")
lines(x=1996:2024,y=threats.time.reptile[,4],type = "b",pch=18,col="gray30")
lines(x=1996:2024,y=threats.time.reptile[,5],type = "b",pch=19,col="gray40")






threats.time.all <- threats.time.mammal+threats.time.birds+threats.time.amphibian+threats.time.fish+threats.time.flower+threats.time.gastropode+threats.time.reptile









