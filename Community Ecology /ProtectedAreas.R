# Load packages

library(tidyverse)
library(maps) # allows access to space polygons corresponding to diff countries and regions
library(sf)
library(vegan)
library(betapart)

# load csv files
bird.dat.tot <- read.csv("distributions.csv")
coords.dat <- read.csv("grid_square_coordinates_lookup.csv")

View(bird.dat.tot)
View(coords.dat)

# use maps to visualise the UK

map(database="world", regions="UK")

#superimpose coordinates 

points(coords.dat$long, coords.dat$lat)

#store polygons defining the UK and convert them to a "simple feature"
UK <- map(database="world", regions="UK", fill=TRUE)
UK.sf <- st_as_sf(UK)
UK.sf <- st_cast(UK.sf, "POLYGON") #ignore warning
plot(st_geometry(UK.sf))
points(coords.dat)

# keep grid cells with at least one corner on the main island
coords.sf <- st_as_sf(coords.dat[,1:2],coords = c("long", "lat"), crs = st_crs(UK.sf))

# isolate points intersecting the different polygons representing the UK
# only keep the points intersecting the 15th polygon, corresponding to the main island

coords.dat.island.ind <- st_intersects(UK.sf,coords.sf)[[15]]
coords.dat.island <- coords.dat[coords.dat.island.ind,]

# visualise the output to confirm you kept the correct points
plot(st_geometry(UK.sf))
points(coords.dat.island)

## SITE BY SPECIES MATRIX ##

# create one data frame for summer and for winter
##select bird data for the desired period
bird.dat <- bird.dat.tot[which(bird.dat.tot$period=="2008-11" | bird.dat.tot$period=="2007/08-10/11"),c("period","speccode","grid")]

dim(coords.dat.island)
coords.dat.island <- coords.dat.island[which(coords.dat.island$grid %in% bird.dat$grid),]
dim(coords.dat.island)
bird.dat <- bird.dat[which(bird.dat$grid %in% coords.dat.island$grid),]

bird.summer.dat <- bird.dat[which(bird.dat$period=="2008-11"),c("speccode","grid")]
bird.winter.dat <- bird.dat[which(bird.dat$period=="2007/08-10/11"),c("speccode","grid")]

# transform the dataframes into site by species data frames, using tidyverse

bird.summer.dat$presence <- 1 ##add a column with the values to populate the site-by-species data frame
bird.summer.dat.pa <- bird.summer.dat %>% 
  pivot_wider(names_from=speccode,values_from=c(presence)) ##site-by-species data frames with NAs
list0 <- as.list(rep(0,ncol(bird.summer.dat.pa))) ##values to replace the NAs
names(list0) <- names(bird.summer.dat.pa)
bird.summer.dat.pa <- as.data.frame(bird.summer.dat.pa %>% replace_na(list0)) ##replace the NAs by 0’s
row.names(bird.summer.dat.pa) <- bird.summer.dat.pa$grid ##change row names
bird.summer.dat.pa <- bird.summer.dat.pa[,-1] ##remove the first column with site names
bird.summer.dat.pa <- bird.summer.dat.pa[order(row.names(bird.summer.dat.pa)),] ##sort by grid cell names

# do it for winter
bird.winter.dat$presence <- 1
bird.winter.dat.pa <- bird.winter.dat %>% 
  pivot_wider(names_from=speccode,values_from=c(presence))
list0 <- as.list(rep(0,ncol(bird.winter.dat.pa)))
names(list0) <- names(bird.winter.dat.pa)
bird.winter.dat.pa <- as.data.frame(bird.winter.dat.pa %>% replace_na(list0))
row.names(bird.winter.dat.pa) <- bird.winter.dat.pa$grid
bird.winter.dat.pa <- bird.winter.dat.pa[,-1]
bird.winter.dat.pa <- bird.winter.dat.pa[order(row.names(bird.winter.dat.pa)),]

# confirm there are no empty cells, nor species occurring in no cell
which(colSums(bird.summer.dat.pa)==0)
which(rowSums(bird.summer.dat.pa)==0)
which(colSums(bird.winter.dat.pa)==0)
which(rowSums(bird.winter.dat.pa)==0)

## COMPUTE GENERAL PATTERNS ##
# use code from practical 1

## Gamma diversity, ncol() 
ncol(bird.summer.dat.pa)
ncol(bird.winter.dat.pa)

## Alpha diversity, mean(rowSums())
mean(rowSums(bird.summer.dat.pa))
mean(rowSums(bird.winter.dat.pa))

# plot distributions of alpha diversity values using hist()
par(mfrow=c(1,2))
hist(rowSums(bird.summer.dat.pa),main="Summer",xlab="Richness")
hist(rowSums(bird.winter.dat.pa),main="Winter",xlab="Richness")

## Species accumulation curves, specaccum() for observed SAC and poolaccum() for Chao2 estimate
SAC.summer <- specaccum(bird.summer.dat.pa)
SAC.winter <- specaccum(bird.winter.dat.pa)

plot(SAC.summer)
plot(SAC.winter)

Estim.summer <- poolaccum(bird.summer.dat.pa)
Estim.winter <- poolaccum(bird.winter.dat.pa)


par(mfrow=c(1,2))
plot(SAC.summer$richness,pch=1,lty=1,lwd=2,type="b",col="honeydew3",ylim=c(0,max(rowMeans(Estim.summer$chao))),ylab="Richness",main="Summer")
points(3:nrow(bird.summer.dat.pa),rowMeans(Estim.summer$chao),pch=2,lty=2,lwd=2,type="b",col="lavender")
plot(SAC.winter$richness,pch=1,lty=1,lwd=2,type="b",col="honeydew3",ylim=c(0,max(rowMeans(Estim.winter$chao))),ylab="Richness",main="Winter")
points(3:nrow(bird.winter.dat.pa),rowMeans(Estim.winter$chao),pch=2,lty=2,lwd=2,type="b",col="lavender")

last(rowMeans(Estim.summer$chao))/last(SAC.summer$richness) # calculates the slopes
last(rowMeans(Estim.winter$chao))/last(SAC.winter$richness)
# we shouldn't use Chao2 estimate as we have sampled all areas, cannot sample anymore to get closer to the line

## Beta diversity, beta.pair()

beta.summer <- beta.pair(bird.summer.dat.pa)
beta.winter <- beta.pair(bird.winter.dat.pa)

# compare Simpson and Sorensen dissimilarity indices

mean(beta.summer$beta.sim)
mean(beta.winter$beta.sim)

mean(beta.summer$beta.sor)
mean(beta.winter$beta.sor)

# plot their distributions
par(mfrow=c(1,2))
hist(beta.summer$beta.sim)
hist(beta.winter$beta.sim)

par(mfrow=c(1,2))
hist(beta.summer$beta.sor)
hist(beta.winter$beta.sor)

## Occupancy-frequency distribution
# compute the histogram of the species occupancy (number of grid cells they occur in)
# compute with colSums() 

freq.summer <- colSums(bird.summer.dat.pa)/nrow(bird.summer.dat.pa)
freq.winter <- colSums(bird.winter.dat.pa)/nrow(bird.winter.dat.pa)

par(mfrow=c(1,2))
hist(freq.summer,main="Summer",xlab="Occupancy",breaks = seq(0,1,0.1))
hist(freq.winter,main="Winter",xlab="Occupancy",breaks = seq(0,1,0.1))

## SPECIAL PROTECTED AREAS ##

# read shapefile
SPA <- read_sf(dsn = "GB-SPA-OSGB36-20220930")

# need to keep grid cells with at least one corner in the poolygons
coords.island.sf <- st_as_sf((coords.dat.island[,1:2]),coords = c("long", "lat"), crs = 4326) #crs - 4326for WGS84
coords.island.sf <- st_transform(coords.island.sf,crs=st_crs(SPA)) ##transform the projection to be the same as SPA
coords.dat.keep.SPA.ind <- unlist(st_intersects(SPA,coords.island.sf)) ##this time, we want to keep the points in all the polygons, not just the first one. We use unlist to collate all point indices. 
coords.dat.keep.SPA <- coords.dat.island[coords.dat.keep.SPA.ind,]

# visualise the output
plot(st_geometry(UK.sf))
points(coords.dat.keep.SPA)

# extract names of the grid points
grid.SPA <- unique(coords.dat.keep.SPA$grid)

# And we only keep the corresponding rows:

bird.summer.dat.pa.SPA <- bird.summer.dat.pa[which(row.names(bird.summer.dat.pa) %in% grid.SPA),]
bird.winter.dat.pa.SPA <- bird.winter.dat.pa[which(row.names(bird.winter.dat.pa) %in% grid.SPA),]

# Since we selected a limited number of sites, some species do not occur in the remaining cells. 
# But all cells should have species in them. You can verify this with:
  
which(colSums(bird.summer.dat.pa.SPA)==0)
which(rowSums(bird.summer.dat.pa.SPA)==0)
which(colSums(bird.winter.dat.pa.SPA)==0)
which(rowSums(bird.winter.dat.pa.SPA)==0)

# We therefore need to remove the columns corresponding to absent species, i.e. the columns whose sum is 0. 

bird.summer.dat.pa.SPA <- bird.summer.dat.pa.SPA[,-which(colSums(bird.summer.dat.pa.SPA)==0)]
bird.winter.dat.pa.SPA <- bird.winter.dat.pa.SPA[,-which(colSums(bird.winter.dat.pa.SPA)==0)]

## Compute the same patterns as before for this data

## Gamma diversity, ncol() 
ncol(bird.summer.dat.pa.SPA)
ncol(bird.winter.dat.pa.SPA)

## Alpha diversity, mean(rowSums())
mean(rowSums(bird.summer.dat.pa.SPA))
mean(rowSums(bird.winter.dat.pa.SPA))

# plot distributions of alpha diversity values using hist()
par(mfrow=c(1,2))
hist(rowSums(bird.summer.dat.pa.SPA),main="Summer SPA",xlab="Richness")
hist(rowSums(bird.winter.dat.pa.SPA),main="Winter SPA",xlab="Richness")

## Species accumulation curves, specaccum() for observed SAC and poolaccum() for Chao2 estimate
SAC.summer.SPA <- specaccum(bird.summer.dat.pa.SPA)
SAC.winter.SPA <- specaccum(bird.winter.dat.pa.SPA)

plot(SAC.summer.SPA)
plot(SAC.winter.SPA)

Estim.summer.SPA <- poolaccum(bird.summer.dat.pa.SPA)
Estim.winter.SPA <- poolaccum(bird.winter.dat.pa.SPA)


par(mfrow=c(1,2))
plot(SAC.summer.SPA$richness,pch=1,lty=1,lwd=2,type="b",col="honeydew3",ylim=c(0,max(rowMeans(Estim.summer.SPA$chao))),ylab="Richness",main="Summer SPA")
points(3:nrow(bird.summer.dat.pa.SPA),rowMeans(Estim.summer.SPA$chao),pch=2,lty=2,lwd=2,type="b",col="lavender")
plot(SAC.winter.SPA$richness,pch=1,lty=1,lwd=2,type="b",col="honeydew3",ylim=c(0,max(rowMeans(Estim.winter.SPA$chao))),ylab="Richness",main="Winter SPA")
points(3:nrow(bird.winter.dat.pa.SPA),rowMeans(Estim.winter.SPA$chao),pch=2,lty=2,lwd=2,type="b",col="lavender")

last(rowMeans(Estim.summer.SPA$chao))/last(SAC.summer.SPA$richness) # calculates the slopes
last(rowMeans(Estim.winter.SPA$chao))/last(SAC.winter.SPA$richness)
# we shouldn't use Chao2 estimate as we have sampled all areas, cannot sample anymore to get closer to the line

## Beta diversity, beta.pair()

beta.summer.SPA <- beta.pair(bird.summer.dat.pa.SPA)
beta.winter.SPA <- beta.pair(bird.winter.dat.pa.SPA)

# compare Simpson and Sorensen dissimilarity indices

mean(beta.summer.SPA$beta.sim)
mean(beta.winter.SPA$beta.sim)

mean(beta.summer.SPA$beta.sor)
mean(beta.winter.SPA$beta.sor)

# plot their distributions
par(mfrow=c(1,2))
hist(beta.summer.SPA$beta.sim)
hist(beta.winter.SPA$beta.sim)

par(mfrow=c(1,2))
hist(beta.summer.SPA$beta.sor)
hist(beta.winter.SPA$beta.sor)

## Occupancy-frequency distribution
# compute the histogram of the species occupancy (number of grid cells they occur in)
# compute with colSums() 

freq.summer.SPA <- colSums(bird.summer.dat.pa.SPA)/nrow(bird.summer.dat.pa.SPA)
freq.winter.SPA <- colSums(bird.winter.dat.pa.SPA)/nrow(bird.winter.dat.pa.SPA)

par(mfrow=c(1,2))
hist(freq.summer.SPA,main="Summer SPA",xlab="Occupancy",breaks = seq(0,1,0.1))
hist(freq.winter.SPA,main="Winter SPA",xlab="Occupancy",breaks = seq(0,1,0.1))


## COMPARE TO RANDOM DATA
#it may make more sense to compare the patterns to the same number of randomly selected grid cells. 
# You simply need to use the function sample() 

sites.summer <- sort(sample(nrow(bird.summer.dat.pa),nrow(bird.summer.dat.pa.SPA)))
sites.winter <- sort(sample(nrow(bird.winter.dat.pa),nrow(bird.winter.dat.pa.SPA)))

bird.summer.dat.pa.SS <- bird.summer.dat.pa[sites.summer,]
bird.winter.dat.pa.SS <- bird.winter.dat.pa[sites.winter,]

# remove sites with species not occurring in any of these cells
which(colSums(bird.summer.dat.pa.SS)==0)
which(rowSums(bird.summer.dat.pa.SS)==0)
which(colSums(bird.summer.dat.pa.SS)==0)
which(rowSums(bird.summer.dat.pa.SS)==0)

bird.summer.dat.pa.SS <- bird.summer.dat.pa.SS[,-which(colSums(bird.summer.dat.pa.SS)==0)]
bird.winter.dat.pa.SS <- bird.winter.dat.pa.SS[,-which(colSums(bird.winter.dat.pa.SS)==0)]

# We can then compute the patterns, and compare them to the patterns obtained for the special protected areas
