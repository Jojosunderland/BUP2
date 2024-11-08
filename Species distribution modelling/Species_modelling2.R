## Practical 2

#install packages
install.packages("geodata",dependencies=TRUE,repos="https://cloud.r-project.org")
install.packages("predicts",dependencies=TRUE,repos="https://cloud.r-project.org")
install.packages("terra",dependencies=TRUE,repos="https://cloud.r-project.org")

#load packages
library(geodata)
library(predicts)
library(terra)

# species chosen:3 toed sloth
# download species data
occdata <- geodata::sp_occurrence("Bradypus", "Linnaeus*", 
                                  geo=FALSE,removeZeros=TRUE,start=1,end=10000)

dim(occdata)

occdata[1:10,]

#plot the global distribution, make sure it fits with expectation
wrld <- world(path=".")
#this function gives us an outline of the world's political boundaries. Reminder, if ever you want to know more about an R function, you can write ?function.name, e.g., ?world
plot(wrld, xlim=c(-180,180), ylim=c(-80,80), col="light yellow", border="light gray")
# add the points
points(occdata$lon, occdata$lat, col='blue', pch=20)

## Cleaning up occurrence data
occdata<-subset(occdata, lat<20) # removed point from N america

dups <- duplicated(occdata[, c('lon', 'lat')])
#This identifies observations that have already appeared above
sum(dups)

#remove duplicates
occ <- occdata[!dups, ]

# download worldclim data
output_dir<-"~/Documents/WorkingD/BUP/Species distribution modelling"

bio_glob<-worldclim_global(var="bio", res=10,path=output_dir, version="2.1")

dim(bio_glob)

#we will also clip the spatraster so it only covers the spatial extent of our study species. First its longitudes then latitudes
summary(occ$lon)
summary(occ$lat)

e <- ext(-90, -30, -55,20) #min/max lon/lat values


predictors <- crop(bio_glob, e)

names(predictors)<-substring(names(predictors),11,16)
#here we're just shortening the names of predictors by taking the 11th to 16th characters.

#plot first 9 climate variables
plot(predictors,1:9)

# plot climate data for first variable
plot(predictors,1)
points(occ$lon,occ$lat, col='firebrick1',pch=16,cex=0.2)

#here I'm setting the spatial extent to be broadly consistent with that of my study species (you need to make sure it is sampling from the same extent). Remember to find out how a function works you can do ?function
bg<-spatSample(predictors,5000,"random", na.rm=TRUE, as.points=TRUE,ext=e)

#Here we'll plot our background points on a map of climwin variable 1 (you could change this to any of the worldclim variables)
plot(predictors, 1)
points(bg, cex=0.1, col='red')

##Matching occurrence and climate data

occlatlon<-cbind(occ$lon,occ$lat)
presvals <- extract(predictors, occlatlon)
#presvals is the climate data for where the species is present
backvals <- values(bg)
#backvals is the climate data for the background data
bg_lonlat<-geom(bg)
lonlats<-rbind(occlatlon, bg_lonlat[,c("x","y")])
pb <- c(rep(1, nrow(presvals)), rep(0, nrow(backvals)))
#The first column of the dataset is a vector of 1s for presences and 0s for background data.
sdmdata <- data.frame(cbind(lonlats,pb, rbind(presvals, backvals)))
#here we combine the presence and background data into a single data frame

## examine how colinear predictor variables are
pairs(sdmdata[,4:7], cex=0.1)

## Fitting a species distribution model

sdmdata<-subset(sdmdata,is.na(bio_1)==F)
#here we're just removing a couple of rows where the climate data are NAs.


specdata<-as.data.frame(cbind(rep("Bradypus Linnaeus",length(sdmdata[,1])),
                              sdmdata))

names(specdata)[1:4]<-c("species","longitude","latitude","presabs")

specdata<-subset(specdata,presabs==1)

backdata<-as.data.frame(cbind(rep("background",length(sdmdata[,1])),
                              sdmdata))

names(backdata)[1:4]<-c("","longitude","latitude","presabs")

backdata<-subset(backdata,presabs==0)


write.table(specdata[,-4],paste(output_dir,"/BradypusLinnaeus_swd.csv",sep=""),col.names=T,row.names=F,sep=",")
write.table(backdata[,-4],paste(output_dir,"/background.csv",sep=""),col.names=T,row.names=F,sep=",")

model<-MaxEnt(sdmdata[,-c(1:3)],sdmdata[,3],removeDuplicates=TRUE)
#Here we've used all of the climate variables, but you could be more discerning. 
#We've also asked the model to ignore any data that comes from the same cell.
#In the maxent call we first specify the climate data - these are in all the columns except the first one.
#Next we specify the presence/background data - column1 i.e. [,3]
model
plot(model)

# look at the predicted climate suitability globally, and see how it matches where the species has been recorded
predictedocc <- predict(model, predictors, args=c("outputformat=raw")) 

par(mfrow=c(2,1))
plot(predictedocc)
plot(predictedocc)
points(occlatlon,pch=".", col='pink')

## Predicting future distributions
# we need future climate data

# download climate data
bio_fut<-cmip6_world(model='ACCESS-ESM1-5', ssp='245', time='2041-2060', var='bioc', res=10, path=output_dir)

fut_predictors<-crop(bio_fut,e)

# compare baseline to future climates for any of our variables
plot(predictors,19)
plot(fut_predictors,19) # 19 was my species most important variable

# generate a future prediction for climate suitability for my species
names(fut_predictors)<-names(predictors)

fut_predictedocc <- predict(model, fut_predictors, args=c("outputformat=raw")) 

par(mfrow=c(2,1))
plot(predictedocc,main="current")

plot(fut_predictedocc,main="2050")
# the area on the map that is light blue is smaller in 2050, suggesting species distribution is going to decline
