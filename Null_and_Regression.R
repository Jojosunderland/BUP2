## NULL AND REGRESSION MODELS ##

install.packages("gam")
install.packages("scam")
install.packages("gdm")
install.packages("zetadiv")

library(vegan)
library(gam)
library(scam)
library(gdm)
library(car)
library(betapart)
library(tidyverse)
library(zetadiv)

## Null models - permutation algorithms ##

# load data
load("islands_Soc_Haw_null_models_practical.RData")

# beta diversity - beta.pair and dat.soc.pa
View(dat.soc.pa)
beta.pair(dat.soc.pa)

# compute the average beta diversity across islands for each diversity index 

beta.soc.sim <- mean(beta.pair(dat.soc.pa)$beta.sim) #simpson
beta.soc.sor <- mean(beta.pair(dat.soc.pa)$beta.sor) #sorensen

beta.soc.sim
beta.soc.sor
#alternative code:
beta.soc.obs <- c(mean(beta.pair(dat.soc.pa)$beta.sor),mean(beta.pair(dat.soc.pa)$beta.sim))
beta.soc.obs

# We now need to assess if these values are different from what we can expect from chance alone. 
# we will use permutation algorithm to generate beta diversity values under randomness

?permatfull

dat.soc.pa.perm.none <- permatfull(dat.soc.pa, times = 99, mtype ="prab", fixedmar = "none") #fixedmar is the constraints on the matrix, none, fix rows, fix columns, fix both
dat.soc.pa.perm.row <- permatfull(dat.soc.pa, times = 99, mtype ="prab", fixedmar = "rows")
dat.soc.pa.perm.col <- permatfull(dat.soc.pa, times = 99, mtype ="prab", fixedmar = "columns")
dat.soc.pa.perm.both <- permatfull(dat.soc.pa, times = 99, mtype ="prab", fixedmar = "both")

View(dat.soc.pa.perm.none)
View(dat.soc.pa.perm.row)
View(dat.soc.pa.perm.col)
View(dat.soc.pa.perm.both)

# compute beta diversity indices for each permutated matrix, and then the averages

beta.mean <- function(dat){ ##this function computes the average for Sorensen and Simpson beta diversity for a givent site-by-species matrix dat
  beta.dat <- beta.pair(dat)
  return(c(mean(beta.dat$beta.sor),mean(beta.dat$beta.sim)))
}

beta.rand.soc.none <- data.frame(matrix(unlist(lapply(dat.soc.pa.perm.none$perm,beta.mean)),99,2,byrow = TRUE)) ## this applies the beta.mean() function above to each permutated site-by-species matrix
names(beta.rand.soc.none) <- c("Sorensen","Simpson") ## rename the columns of the data frame

beta.rand.soc.row <- data.frame(matrix(unlist(lapply(dat.soc.pa.perm.row$perm,beta.mean)),99,2,byrow = TRUE))
names(beta.rand.soc.row) <- c("Sorensen","Simpson")

beta.rand.soc.col <- data.frame(matrix(unlist(lapply(dat.soc.pa.perm.col$perm,beta.mean)),99,2,byrow = TRUE))
names(beta.rand.soc.col) <- c("Sorensen","Simpson")

beta.rand.soc.both <- data.frame(matrix(unlist(lapply(dat.soc.pa.perm.both$perm,beta.mean)),99,2,byrow = TRUE))
names(beta.rand.soc.both) <- c("Sorensen","Simpson")

# plot the distributions under different permutation algorithms, along with obs beta diversity values for Sor and Sim with 5/95% quartiles

par(mfrow=c(2,4))
hist(beta.rand.soc.none$Sorensen,breaks=seq(0.5,1,0.01),main="None fixed",xlab="Sorensen")
abline(v=quantile(beta.rand.soc.none$Sorensen,0.025),col="blue")
abline(v=quantile(beta.rand.soc.none$Sorensen,0.0975),col="blue")
abline(v=beta.soc.obs[1],col="red")
hist(beta.rand.soc.row$Sorensen,breaks=seq(0.5,1,0.01),main="Rows fixed",xlab="Sorensen")
abline(v=quantile(beta.rand.soc.row$Sorensen,0.025),col="blue")
abline(v=quantile(beta.rand.soc.row$Sorensen,0.0975),col="blue")
abline(v=beta.soc.obs[1],col="red")
hist(beta.rand.soc.col$Sorensen,breaks=seq(0.5,1,0.01),main="Columns fixed",xlab="Sorensen")
abline(v=quantile(beta.rand.soc.col$Sorensen,0.025),col="blue")
abline(v=quantile(beta.rand.soc.col$Sorensen,0.0975),col="blue")
abline(v=beta.soc.obs[1],col="red")
hist(beta.rand.soc.both$Sorensen,breaks=seq(0.5,1,0.01),main="Both fixed",xlab="Sorensen")
abline(v=quantile(beta.rand.soc.both$Sorensen,0.025),col="blue")
abline(v=quantile(beta.rand.soc.both$Sorensen,0.0975),col="blue")
abline(v=beta.soc.obs[1],col="red")

hist(beta.rand.soc.none$Simpson,breaks=seq(0.2,1,0.01),main="",xlab="Simpson")
abline(v=quantile(beta.rand.soc.none$Simpson,0.025),col="blue")
abline(v=quantile(beta.rand.soc.none$Simpson,0.0975),col="blue")
abline(v=beta.soc.obs[2],col="red")
hist(beta.rand.soc.row$Simpson,breaks=seq(0.2,1,0.01),main="",xlab="Simpson")
abline(v=quantile(beta.rand.soc.row$Simpson,0.025),col="blue")
abline(v=quantile(beta.rand.soc.row$Simpson,0.0975),col="blue")
abline(v=beta.soc.obs[2],col="red")
hist(beta.rand.soc.col$Simpson,breaks=seq(0.2,1,0.01),main="",xlab="Simpson")
abline(v=quantile(beta.rand.soc.col$Simpson,0.025),col="blue")
abline(v=quantile(beta.rand.soc.col$Simpson,0.0975),col="blue")
abline(v=beta.soc.obs[2],col="red")
hist(beta.rand.soc.both$Simpson,breaks=seq(0.2,1,0.01),main="",xlab="Simpson")
abline(v=quantile(beta.rand.soc.both$Simpson,0.025),col="blue")
abline(v=quantile(beta.rand.soc.both$Simpson,0.0975),col="blue")
abline(v=beta.soc.obs[2],col="red")


## do it for Hawai'i

# load lists for Hawai'i
load("permutations_hawaii.RData") 
View(dat.haw.pa)

# Simpson and sorensen averages
beta.haw.obs <- c(mean(beta.pair(dat.haw.pa)$beta.sor),mean(beta.pair(dat.haw.pa)$beta.sim))

beta.rand.haw.none <- data.frame(matrix(unlist(lapply(dat.haw.pa.perm.none$perm,beta.mean)),99,2,byrow = TRUE))
names(beta.rand.haw.none) <- c("Sorensen","Simpson")

beta.rand.haw.row <- data.frame(matrix(unlist(lapply(dat.haw.pa.perm.row$perm,beta.mean)),99,2,byrow = TRUE))
names(beta.rand.haw.row) <- c("Sorensen","Simpson")

beta.rand.haw.col <- data.frame(matrix(unlist(lapply(dat.haw.pa.perm.col$perm,beta.mean)),99,2,byrow = TRUE))
names(beta.rand.haw.col) <- c("Sorensen","Simpson")

beta.rand.haw.both <- data.frame(matrix(unlist(lapply(dat.haw.pa.perm.both$perm,beta.mean)),99,2,byrow = TRUE))
names(beta.rand.haw.both) <- c("Sorensen","Simpson")



par(mfrow=c(2,4))
hist(beta.rand.haw.none$Sorensen,breaks=seq(0,1,0.025),main="None fixed",xlab="Sorensen")
abline(v=quantile(beta.rand.haw.none$Sorensen,0.025),col="blue")
abline(v=quantile(beta.rand.haw.none$Sorensen,0.0975),col="blue")
abline(v=beta.haw.obs[1],col="red")
hist(beta.rand.haw.row$Sorensen,breaks=seq(0,1,0.025),main="Rows fixed",xlab="Sorensen")
abline(v=quantile(beta.rand.haw.row$Sorensen,0.025),col="blue")
abline(v=quantile(beta.rand.haw.row$Sorensen,0.0975),col="blue")
abline(v=beta.haw.obs[1],col="red")
hist(beta.rand.haw.col$Sorensen,breaks=seq(0,1,0.025),main="Columns fixed",xlab="Sorensen")
abline(v=quantile(beta.rand.haw.col$Sorensen,0.025),col="blue")
abline(v=quantile(beta.rand.haw.col$Sorensen,0.0975),col="blue")
abline(v=beta.haw.obs[1],col="red")
hist(beta.rand.haw.both$Sorensen,breaks=seq(0,1,0.025),main="Both fixed",xlab="Sorensen")
abline(v=quantile(beta.rand.haw.both$Sorensen,0.025),col="blue")
abline(v=quantile(beta.rand.haw.both$Sorensen,0.0975),col="blue")
abline(v=beta.haw.obs[1],col="red")

hist(beta.rand.haw.none$Simpson,breaks=seq(0,1,0.025),main="",xlab="Simpson")
abline(v=quantile(beta.rand.haw.none$Simpson,0.025),col="blue")
abline(v=quantile(beta.rand.haw.none$Simpson,0.0975),col="blue")
abline(v=beta.haw.obs[2],col="red")
hist(beta.rand.haw.row$Simpson,breaks=seq(0,1,0.025),main="",xlab="Simpson")
abline(v=quantile(beta.rand.haw.row$Simpson,0.025),col="blue")
abline(v=quantile(beta.rand.haw.row$Simpson,0.0975),col="blue")
abline(v=beta.haw.obs[2],col="red")
hist(beta.rand.haw.col$Simpson,breaks=seq(0,1,0.025),main="",xlab="Simpson")
abline(v=quantile(beta.rand.haw.col$Simpson,0.025),col="blue")
abline(v=quantile(beta.rand.haw.col$Simpson,0.0975),col="blue")
abline(v=beta.haw.obs[2],col="red")
hist(beta.rand.haw.both$Simpson,breaks=seq(0,1,0.025),main="",xlab="Simpson")
abline(v=quantile(beta.rand.haw.both$Simpson,0.025),col="blue")
abline(v=quantile(beta.rand.haw.both$Simpson,0.0975),col="blue")
abline(v=beta.haw.obs[2],col="red")

## REGRESSION MODELS ON SPECIES RICHNESS ##

# data for the whole french polynesia 
load("islands_FP_regression_models_practical.RData")

## Generalised linear models ##

##alpha diversity
islands.pred.FP$alpha <- rowSums(dat.FP.pa)

cor(islands.pred.FP[,2:7])

##GLMs
mod.glm.FP <- glm(alpha~elev_max+IslandArea+temp.mean+prec.an+temp.seas+prec.seas,data=islands.pred.FP,family = poisson())
mod.glm.FPb <- glm(alpha~elev_max+IslandArea+temp.seas+prec.an+prec.seas,data=islands.pred.FP,family = poisson())
mod.glm.FPc <- glm(alpha~elev_max+IslandArea+temp.seas+prec.an,data=islands.pred.FP,family = poisson())
mod.glm.FPd <- glm(alpha~elev_max+IslandArea,data=islands.pred.FP,family = poisson())

#inspect model outputs
summary(mod.glm.FP)
summary(mod.glm.FPb)
summary(mod.glm.FPc)
summary(mod.glm.FPd)

# We can also look at how the correlation between predictors affect the output using the variance inflation factor (vif).  
you want the vif values to be below 10. Otherwise, you need to drop some predictors with high vif values. One usually does this one predictor at a time, and re-assess the vif
vif(mod.glm.FP) # temp mean too high, remove it
vif(mod.glm.FPb)

#You can also check different combinations of predictors, and see if more complex models actually explain better the response, using the Akaike Information Criterion (AIC). 
# you only keep a more complex model if its AIC is lower than the AIC of a simpler model by a margin of 2
AIC(mod.glm.FP)
AIC(mod.glm.FPb)
AIC(mod.glm.FPc)
AIC(mod.glm.FPd)

#Finally, we can check how well the GLMs fit the data using the following code:

cor(islands.pred.FP$alpha,predict(mod.glm.FP,type="response"))


## Generalised additive models (GAM) ##

# For island area
mod.gam.FP <- gam(alpha~s(elev_max)+s(IslandArea)+s(temp.seas)+s(prec.an)+s(prec.an)+s(prec.seas),data=islands.pred.FP,family = poisson(),method = "REML")
# the s() function allows to fit the non-linear relationships. These are called “splines”


##Smoother GAMs
mod.gam.FP2 <- gam(alpha~s(elev_max,k=4)+s(IslandArea,k=4)+s(temp.seas,k=4)+s(prec.an,k=4)+s(prec.an,k=4)+s(prec.seas,k=4),data=islands.pred.FP,family = poisson(),method = "REML")

##Model summary
summary(mod.glm.FP)
summary(mod.glm.FPb)
summary(mod.glm.FPc)
summary(mod.gam.FP)
summary(mod.gam.FP2)

#plot the model outputs

##model performance - see how well GAM fits the data
cor(islands.pred.FP$alpha,predict(mod.glm.FP,type="response"))
cor(islands.pred.FP$alpha,predict(mod.glm.FPb,type="response"))
cor(islands.pred.FP$alpha,predict(mod.glm.FPc,type="response"))
cor(islands.pred.FP$alpha,predict(mod.gam.FP,type="response"))
cor(islands.pred.FP$alpha,predict(mod.gam.FP2,type="response"))

par(mfrow=c(1,5))
plot(islands.pred.FP$alpha,predict(mod.glm.FP,type="response"),type="p")
plot(islands.pred.FP$alpha,predict(mod.glm.FPb,type="response"),type="p")
plot(islands.pred.FP$alpha,predict(mod.glm.FPc,type="response"),type="p")
plot(islands.pred.FP$alpha,predict(mod.gam.FP,type="response"),type="p")
plot(islands.pred.FP$alpha,predict(mod.gam.FP2,type="response"),type="p")

##Plot GAM splines
par(mfrow=c(2,5))
plot(mod.gam.FP,ylim=c(-20,20),residuals = T,pch=1)
plot(mod.gam.FP2,ylim=c(-20,20),residuals = T,pch=1)


