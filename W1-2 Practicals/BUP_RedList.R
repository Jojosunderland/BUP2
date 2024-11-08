# set working directory
setwd("~/Documents/WorkingD")

# this code clears environment
rm(list=ls())

?read.csv

# read in the IUCN documents and name them

climate <- read.csv("IUCN data/climate/assessments.csv")
invasion <-  read.csv("IUCN data/invasion/assessments.csv")
landuse <- read.csv("IUCN data/landuse/assessments.csv")
overexploitation <- read.csv("IUCN data/overexploitation/assessments.csv")
pollution <- read.csv("IUCN data/pollution/assessments.csv")

# You can check the class of each variable with either lines of code:
unlist(lapply(pollution, class))
str(pollution)

# selecting a column in the dataframe uses $
climate$scientificName

# we need to add threat category to the data frames, here are two way we can do it
# we can add a single column to each data frame indicating "land use" or "overexploitation" etc, 
# however some species are affected by multiple threats 
# therefore, we will add 5 columns one for each threat, and fill it with 0's (if this is not the threat 
#corresponding to the data frame) or 1's (if this is the threat)

# to do this you need to assign the desired values to the new columns

climate$landuse <- 0
climate$overexploitation <- 0
climate$invasion <- 0
climate$pollution <- 0
climate$climate <- 1

# because a value is assigned to it a column is automatically created for each at the end of the data frame
# do this for each dataframe

invasion$landuse <- 0
invasion$overexploitation <- 0
invasion$invasion <- 1
invasion$pollution <- 0
invasion$climate <- 0

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

pollution$landuse <- 0
pollution$overexploitation <- 0
pollution$invasion <- 0
pollution$pollution <- 1
pollution$climate <- 0


# now we have data frames with all the same variables, so we can combine the data frames together by stacking them
# we will do this using rbind()

threats.dat <- rbind(landuse, overexploitation, invasion, pollution, climate)

# we have one problem however:
# some species have multiple threats so appear multiple times, we can check the number of unique species name

length(unique(threats.dat$scientificName)) # 558 unique species

# which will be less than the total number of rows, can get the dimensions (rows and columns) using dim()
dim(threats.dat) # 919 total rows = duplicates

# we need to merge duplicated rows and delete redundancies
# find rows which are duplicates

threats.dup <- which(duplicated(threats.dat$scientificName))

# this extracts the subset of rows that are duplicates in the first 23 columns
threats.dat[threats.dup,1:23] # but there can still be duplicates in this

# we need unique rows from this so we use distinct() from package dplyr()

library(dplyr)
# install.packages("dplyr")
# library(dplyr)
threats.redundant <- distinct(threats.dat[threats.dup,1:23])

# now we need to merge all rows that are duplicates. 
# we will use a loop to do so, a for loop is used to iterate over a sequence of values
# here we want to iterate over all unique rows (in threats.redundant)
# i will iterate from 1 to the number of rows of threats.redundant with an increment of 1
# whatever is within {} will be iterated in the loop

for(i in 1:nrow(threats.redundant)) {
  
  # we need to first select the rows that correspond to the row i in threats.redundant, using which()
  # saying the redundant duplicate rows are the same as in the rows in the main data frame 
  
  row.number <- which(threats.dat$scientificName == threats.redundant$scientificName[i])
  
  row.selected <- threats.dat[row.number,] # extracts the rows from threats.dat corresponding to the indices stored in row.number
  #these rows are the duplicates that need to be merged
  new.row <-
    data.frame(c(row.selected[1,1:23], colSums(row.selected[,24:28]))) # this line creates a new merged row
  # it does this by taking the first 23 columns from the first row of duplicates and summing the values of columns 24-28 across all duplicated rows
  
  # removing the redundant duplicated rows and then combining the new row with the big data frame
  threats.dat <- threats.dat[-row.number,] # -row.number, syntax means "remove the rows with these indices"
  threats.dat <- rbind(threats.dat,new.row) # this adds the merged rows back into the data frame
}

dim(threats.dat) # number of rows now match the number of unique species

# use barplot to summarise which threat is most importance, do this by comparing number of rows each threat has
# the more rows, the more species it affects
barplot(colSums (threats.dat[,24:28]))

## VENN DIAGRAMS are useful at showing all possible combinations of threats that affects species:

# this installs and loads the package
if(!require(devtools)) install.packages("devtools") 
devtools::install_github("gaospecial/ggVennDiagram")
library("ggVennDiagram")

# the function ggVennDiagram() needs a list as an input, so we need to create one
# create a list containing vectors of species names affected by each threat

x <- list(
  landuse=landuse$scientificName,
  overexploitation=overexploitation$scientificName,
  invasion=invasion$scientificName,
  pollution=pollution$scientificName,
  climate=climate$scientificName)

# then plot the Venn diagram

ggVennDiagram(x, label_alpha = 0, set_color = 
                c("blue", "black", "red", "purple", "green"))

# Venn diagrams can be hard to read, so we can use an UpSet plot using funciton upset() from UpSetR

if (!require(UpSetR)) install.packages("UpSetR")
library(UpSetR)
upset(data=threats.dat[c(3,24:28)], nsets=5)


## PATTERNS FOR DIFFERENT TAXONOMIC GROUPS

# download taxonomic data to the data frames

climate_tax <- read.csv("IUCN data/climate/taxonomy.csv")
invasion_tax <-  read.csv("IUCN data/invasion/taxonomy.csv")
landuse_tax <- read.csv("IUCN data/landuse/taxonomy.csv")
overexploitation_tax <- read.csv("IUCN data/overexploitation/taxonomy.csv")
pollution_tax <- read.csv("IUCN data/pollution/taxonomy.csv")

# combine the dataframes using rbind() and distinct()

all.tax <- rbind(landuse_tax, overexploitation_tax, invasion_tax, 
                 pollution_tax, climate_tax)

all.tax <- distinct(all.tax)

# need to relate it to threats.dat
# make sure the species described by each row of all.tax is also described by the same row number in threats.dat
#sort data frames by species name
all.tax <- all.tax[order(all.tax$scientificName),]
threats.dat <- threats.dat[order(threats.dat$scientificName),]

# check how many classes are included in this dataset and how many species belong to each class

table(all.tax$className)

# analyses for the best represented classes

# we can extract the data from threats.dat for each taxa by selecting the rows in threats.dat that correspond to the desired class in all.tax
threats.dat.mammals <- threats.dat[which(all.tax$className=="MAMMALIA"),]
threats.dat.birds <- threats.dat[which(all.tax$className=="AVES"),]
threats.dat.fish <- threats.dat[which(all.tax$className=="ACTINOPTERYGII"),]
threats.dat.amphibian <- threats.dat[which(all.tax$className=="AMPHIBIA"),]
threats.dat.gastropode <- threats.dat[which(all.tax$className=="GASTROPODA"),]
threats.dat.flower <- threats.dat[which(all.tax$className=="MAGNOLIOPSIDA"),]
threats.dat.reptile <- threats.dat[which(all.tax$className=="REPTILIA"),]

#plot the upset plots for each taxonomic groups and compare visually to each other
# are the different taxonomic groups affected by the same threats? 

upset(data=threats.dat[c(3,24:28)],nsets=5)
upset(data=threats.dat.mammals[c(3,24:28)],nsets=5)
upset(data=threats.dat.birds[c(3,24:28)],nsets=5)
upset(data=threats.dat.fish[c(3,24:28)],nsets=5)
upset(data=threats.dat.amphibian[c(3,24:28)],nsets=5)
upset(data=threats.dat.gastropode[c(3,24:28)],nsets=5)
upset(data=threats.dat.flower[c(3,24:28)],nsets=5)
upset(data=threats.dat.reptile[c(3,24:28)],nsets=5)


## TEMPORAL PATTERNS
# look at temporal trends for all data and taxa group, to see if there are changes in threat prevalence over another 
# will use the column yearPublished

# first look at the period covered
range(threats.dat$yearPublished)

# we have data from 1996-2024, so we will use a for loop to select the data published before each year

# first select only the data before 1996
threats.dat.temp <- threats.dat[which(threats.dat$yearPublished <= 1996),] 
# the line above selects rows from threats.dat where the publication is less than or equalto 1996
# it is then repeated for all taxa groups
threats.dat.temp.mammal <- threats.dat.mammals[which(threats.dat.mammals$yearPublished <= 1996),]
threats.dat.temp.birds <- threats.dat.birds[which(threats.dat.birds$yearPublished <= 1996),]
threats.dat.temp.fish <- threats.dat.fish[which(threats.dat.fish$yearPublished <= 1996),]
threats.dat.temp.amphibian <- threats.dat.amphibian[which(threats.dat.amphibian$yearPublished <= 1996),]
threats.dat.temp.gastropode <- threats.dat.gastropode[which(threats.dat.gastropode$yearPublished <= 1996),]
threats.dat.temp.flower <- threats.dat.flower[which(threats.dat.flower$yearPublished <= 1996),]
threats.dat.temp.reptile <- threats.dat.reptile[which(threats.dat.reptile$yearPublished <= 1996),]

# then compute the sum of the columns 24-28 of the df and store in a new variable
threats.time <- colSums(threats.dat.temp[,24:28])
# this line calculates the sum of the values in columns 24-48 for all records up to 1996
# same sum operation is done for each taxa dataset
threats.time.mammal <- colSums(threats.dat.temp.mammal[,24:28])
threats.time.birds <- colSums(threats.dat.temp.birds[,24:28])
threats.time.fish <- colSums(threats.dat.temp.fish[,24:28])
threats.time.amphibian <- colSums(threats.dat.temp.amphibian[,24:28])
threats.time.gastropode <- colSums(threats.dat.temp.gastropode[,24:28])
threats.time.flower <- colSums(threats.dat.temp.flower[,24:28])
threats.time.reptile <- colSums(threats.dat.temp.reptile[,24:28])

#then start a loop for years 1997 to 2024
#in each iteratino the code selects the data where yearPublished is less than or equal to the current year (y)
# meaning it incudes data up to that year
for(y in 1997:2024){
  # select only the data before the year in question (less than or equal to year = y)
  threats.dat.temp <- threats.dat[which(threats.dat$yearPublished <= y),]
  threats.dat.temp.mammal <- threats.dat.mammals[which(threats.dat.mammals$yearPublished <= y),]
  threats.dat.temp.birds <- threats.dat.birds[which(threats.dat.birds$yearPublished <= y),]
  threats.dat.temp.fish <- threats.dat.fish[which(threats.dat.fish$yearPublished <= y),]
  threats.dat.temp.amphibian <- threats.dat.amphibian[which(threats.dat.amphibian$yearPublished <= y),]
  threats.dat.temp.gastropode <- threats.dat.gastropode[which(threats.dat.gastropode$yearPublished <= y),]
  threats.dat.temp.flower <- threats.dat.flower[which(threats.dat.flower$yearPublished <= y),]
  threats.dat.temp.reptile <- threats.dat.reptile[which(threats.dat.reptile$yearPublished <= y),]
  
  #compute the sum of columns 24-28 of the df and append it to threats.time using rbind()
  threats.time <- rbind(threats.time,colSums(threats.dat.temp[,24:28]))
  threats.time.mammal <- rbind(threats.time.mammal,colSums(threats.dat.temp.mammal[,24:28]))
  threats.time.birds <- rbind(threats.time.birds,colSums(threats.dat.temp.birds[,24:28]))
  threats.time.fish <- rbind(threats.time.fish,colSums(threats.dat.temp.fish[,24:28]))
  threats.time.amphibian <- rbind(threats.time.amphibian,colSums(threats.dat.temp.amphibian[,24:28]))
  threats.time.gastropode <- rbind(threats.time.gastropode,colSums(threats.dat.temp.gastropode[,24:28]))
  threats.time.flower <- rbind(threats.time.flower,colSums(threats.dat.temp.flower[,24:28]))
  threats.time.reptile <- rbind(threats.time.reptile,colSums(threats.dat.temp.reptile[,24:28]))
}
# the loop has progessively calculated the cumulative sum of threats up to each year and appended the rsult to a data frame
# this process was repeated for each taxa group, allowing for a comparison of threat trends over time

# then plot changes in the number of species impacted by each threat over the years

quartz()
par(mfrow=c(2,4)) # if you run these with all the code below it creates a combined document with all plots together

# for all species
plot(x=1996:2024,y=threats.time[,1],type = "b",pch=15,ylim=c(1,max(threats.time)),xlab="Year",ylab="Number of species",main="All species",log="y")
lines(x=1996:2024,y=threats.time[,2],type = "b",pch=16,col="gray10")
lines(x=1996:2024,y=threats.time[,3],type = "b",pch=17,col="gray20")
lines(x=1996:2024,y=threats.time[,4],type = "b",pch=18,col="gray30")
lines(x=1996:2024,y=threats.time[,5],type = "b",pch=19,col="gray40")
legend(x="topleft",legend=c("landuse","overexploitation","invasion","pollution","climate"),bty="n",lty=1,pch=15:19,col = c("black","gray10","gray20","gray30","gray40"))

# for mammals
plot(x=1996:2024,y=threats.time.mammal[,1],type = "b",pch=15,ylim=c(1,max(threats.time)),xlab="Year",ylab="Number of species",main="Mammals",log="y")
lines(x=1996:2024,y=threats.time.mammal[,2],type = "b",pch=16,col="gray10")
lines(x=1996:2024,y=threats.time.mammal[,3],type = "b",pch=17,col="gray20")
lines(x=1996:2024,y=threats.time.mammal[,4],type = "b",pch=18,col="gray30")
lines(x=1996:2024,y=threats.time.mammal[,5],type = "b",pch=19,col="gray40")

# birds
plot(x=1996:2024,y=threats.time.birds[,1],type = "b",pch=15,ylim=c(1,max(threats.time)),xlab="Year",ylab="Number of species",main="Birds",log="y")
lines(x=1996:2024,y=threats.time.birds[,2],type = "b",pch=16,col="gray10")
lines(x=1996:2024,y=threats.time.birds[,3],type = "b",pch=17,col="gray20")
lines(x=1996:2024,y=threats.time.birds[,4],type = "b",pch=18,col="gray30")
lines(x=1996:2024,y=threats.time.birds[,5],type = "b",pch=19,col="gray40")

# fish
plot(x=1996:2024,y=threats.time.fish[,1],type = "b",pch=15,ylim=c(1,max(threats.time)),xlab="Year",ylab="Number of species",main="Fish",log="y")
lines(x=1996:2024,y=threats.time.fish[,2],type = "b",pch=16,col="gray10")
lines(x=1996:2024,y=threats.time.fish[,3],type = "b",pch=17,col="gray20")
lines(x=1996:2024,y=threats.time.fish[,4],type = "b",pch=18,col="gray30")
lines(x=1996:2024,y=threats.time.fish[,5],type = "b",pch=19,col="gray40")

# amphibians
plot(x=1996:2024,y=threats.time.amphibian[,1],type = "b",pch=15,ylim=c(1,max(threats.time)),xlab="Year",ylab="Number of species",main="Amphibians",log="y")
lines(x=1996:2024,y=threats.time.amphibian[,2],type = "b",pch=16,col="gray10")
lines(x=1996:2024,y=threats.time.amphibian[,3],type = "b",pch=17,col="gray20")
lines(x=1996:2024,y=threats.time.amphibian[,4],type = "b",pch=18,col="gray30")
lines(x=1996:2024,y=threats.time.amphibian[,5],type = "b",pch=19,col="gray40")

#gastropods
plot(x=1996:2024,y=threats.time.gastropode[,1],type = "b",pch=15,ylim=c(1,max(threats.time)),xlab="Year",ylab="Number of species",main="Gastropodes",log="y")
lines(x=1996:2024,y=threats.time.gastropode[,2],type = "b",pch=16,col="gray10")
lines(x=1996:2024,y=threats.time.gastropode[,3],type = "b",pch=17,col="gray20")
lines(x=1996:2024,y=threats.time.gastropode[,4],type = "b",pch=18,col="gray30")
lines(x=1996:2024,y=threats.time.gastropode[,5],type = "b",pch=19,col="gray40")

#flowers
plot(x=1996:2024,y=threats.time.flower[,1],type = "b",pch=15,ylim=c(1,max(threats.time)),xlab="Year",ylab="Number of species",main="Flowers",log="y")
lines(x=1996:2024,y=threats.time.flower[,2],type = "b",pch=16,col="gray10")
lines(x=1996:2024,y=threats.time.flower[,3],type = "b",pch=17,col="gray20")
lines(x=1996:2024,y=threats.time.flower[,4],type = "b",pch=18,col="gray30")
lines(x=1996:2024,y=threats.time.flower[,5],type = "b",pch=19,col="gray40")

#reptiles
plot(x=1996:2024,y=threats.time.reptile[,1],type = "b",pch=15,ylim=c(1,max(threats.time)),xlab="Year",ylab="Number of species",main="Reptiles",log="y")
lines(x=1996:2024,y=threats.time.reptile[,2],type = "b",pch=16,col="gray10")
lines(x=1996:2024,y=threats.time.reptile[,3],type = "b",pch=17,col="gray20")
lines(x=1996:2024,y=threats.time.reptile[,4],type = "b",pch=18,col="gray30")
lines(x=1996:2024,y=threats.time.reptile[,5],type = "b",pch=19,col="gray40")
