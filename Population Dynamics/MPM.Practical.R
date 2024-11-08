## Matrix population modelling

# download and load packages
rm(list=ls(all=TRUE))

my_packages <- c('ggplot2', 'popbio')
new_packages <- my_packages[!(my_packages %in% installed.packages()[,'Package'])]
if(length(new_packages)) install.packages(new_packages, dependencies = T)

library(ggplot2)
library(popbio)

## Parameterising your Matrix Population Model (MPM)
# Stage-specific survival:

#plot survival rates
# first enter the values into a dataframe 
survival <- data.frame(stage=factor(c('Juvenile','Yearling','Adult'), levels=c('Juvenile','Yearling','Adult')), estimate=c(0.463, 0.510, 0.559), lcl=c(0.404, 0.445, 0.499), ucl=c(0.524, 0.574, 0.618))

# then plot by stage
ggplot(survival, aes(stage, estimate, ymin=lcl, ymax=ucl)) + 
  geom_errorbar(width=0.2) + geom_point() + ylim(0,1)

## QUESTION: Which stage has the lowest survival rate? Is this what you would expect?
# Juvenile has the lowest, this is expected as they're more vulnerable on average to predators and other external factors

#Per capita reproduction:

#estimate the number of female offspring produced by each female between censuses (using nest visitation data)
nestdata <- read.table("~/Documents/WorkingD/BUP/Population Dynamics/MPM practical/gjeroynest.txt", header = TRUE, sep = '\t')
head(nestdata)
#clutchno indicates whether it was the first, second or third etc. clutch laid in that nest in that breeding season.
#hatchingsuc indicates whether any live chicks were found for that clutch (yes = 1; no = 0).
#chickno indicates the number of chicks counted on the final visit prior to fledging.

##QUESTION: How might you estimate per capita reproduction from these data?
#calculate the average hatching success and average number of chicks using the mean() function
HatchingSuc <- mean(nestdata$hatchingsuc)
FledglingNo <- mean(nestdata$chickno)

#data frame has a row for each unique nest, column for max number of clutchno for each unique value of nestid
nests <- data.frame(nestid = sort(unique(nestdata$nestid)), numberofclutches=tapply(nestdata$clutchno, nestdata$nestid, max))
#then take the mean to be average number of clutches
ClutchNo <- mean(nests$numberofclutches)

#use these estimates to calculate per capita reproduction (Note! The mean number of chicks prior to fledging is an upwardly biased estimate of the number of fledglings, since not all will likely fledge the nest successfully)
#multiply estimates and divide by 2 as we're only modelling the female segment of the pop
(ClutchNo * HatchingSuc * FledglingNo) / 2 #Our estimate of per capita reproduction R is 1.076


## Deterministic population model

#We now have the numbers we need to parameterise a deterministic model. 
#We need a 3x3 matrix with the fertility transitions along the top row, and the survival transitions on the subsequent rows
# save our estimates of the vital rates
R <- (ClutchNo * HatchingSuc * FledglingNo) / 2
Phi.juv <- survival$estimate[survival$stage=='Juvenile'] 
Phi.yr <- survival$estimate[survival$stage=='Yearling'] 
Phi.ad <- survival$estimate[survival$stage=='Adult'] 

# remind ourselves how these relate to the transition probabilities of the matrix (see slides)
# Juvenile to Juvenile: Phi.juv * R
# Yearling to Juvenile: Phi.yr * R
# Adult to Juvenile: Phi.ad * R
# Juvenile to Yearling: Phi.juv
# Yearling to Yearling: 0 
# Adult to Yearling: 0 
# Juvenile to Adult: 0
# Yearling to Adult: Phi.yr
# Adult to Adult: Phi.ad

# put the transition probabilities into a vector 
sparrowMPM <- c(Phi.juv * R, Phi.yr * R, Phi.ad * R, Phi.juv, 0, 0, 0, Phi.yr, Phi.ad)

# save that vector as a matrix, specifying the number of rows and columns
# use the byrow=TRUE argument to tell R that the first the elements of the vector correspond to the first row of the matrix 
sparrowMPM <- matrix(sparrowMPM, nrow=3, ncol=3, byrow=T)
sparrowMPM

#We can now use the popbio package to do some analyses of our deterministic MPM
lambda(sparrowMPM) #population growth rate (lambda)

## QUESTION: What does this tell us about they dynamics of our sparrow population?
# tells us that our population is growing (at a rate of 1.03 a year)

#Projected Dynamics:
#We can project the dynamics over a given number of iterations (t) based on our matrix and a starting population (n0).
# project over 15 years
t <- 15
# start with 50 juveniles, 20 yearlings and 30 adults
n0 <- c(50,20,30)

# project dynamics 
projection <- pop.projection(sparrowMPM, n0, iterations = t)
projected <- data.frame(time=1:15, N=projection$pop.sizes)

# plot projected pop size over time
ggplot(projected, aes(time, N)) + 
  geom_line() + ylim(0,150) + ylab('Projected N') #The population is projected to increase over time since lambda>1

#Observed Dynamics:
#We can compare this with estimated population counts (N) over the 15 study years, saved in the file ‘popest.txt’
popest <- read.table("~/Documents/WorkingD/BUP/Population Dynamics/MPM practical/popest.txt", header = TRUE, sep = '\t')
head(popest)

# plot N over time
ggplot(popest, aes(year, N)) + 
  geom_line() + ylim(0,200) + ylab('Observed N')

## QUESTION: How does this population trajectory compare with our estimate of lambda?
#the trajectory suggests declines at some time intervals as well as inclines whereas lambda just projected increases
# declines are probably due to stochasticity and environmental factors

#Stable stage distribution and reproductive value:
#popbio has built-in functions for other analyses of the asymptotic dynamics of our matrix
#he stable stage distribution of our population is the long-term average relative abundance of the different stage classes
#the reproductive values of the different stage classes are the expected contribution of each individual in that stage class to future reproduction

stages <- c('Juv','Yr','Ad')
colnames(sparrowMPM) <- stages
rownames(sparrowMPM) <- stages

stable.stage(sparrowMPM)

reproductive.value(sparrowMPM)


## Pertubation analysis
#popbio also includes built-in functions to calculate the sensitivities and elasticities of the different vital rates
#these tell us about the relative importance of each vital rate (or matrix transition) in determining the population growth rate, lambda.

# REMEMBER:
#sensitivities estimate the change in lambda for an absolute change in a vital rate 
#elasticities tell us about the effect of a proportional change

# list the vital rates
sparrow.param <- list(Phi.juv = Phi.juv, Phi.yr = Phi.yr, Phi.ad = Phi.ad, R = R)

# give the matrix equation 
sparrow.equation <- expression(Phi.juv * R, Phi.yr * R, Phi.ad * R, Phi.juv, 0, 0, 0, Phi.yr, Phi.ad)

# run the sensitivity analysis
sens <- vitalsens(sparrow.equation, sparrow.param)
sens

# plot elasticity of the vital rates 
sens$vitalrate <- factor(c('Phi.juv', 'Phi.yr', 'Phi.ad', 'R'), levels = c('Phi.juv', 'Phi.yr', 'Phi.ad', 'R'))
ggplot(sens, aes(vitalrate, elasticity)) + 
  geom_bar(stat = 'identity') 

##QUESTIONa: Which vital rates are most important for population growth? Is this similar to the orca example that we saw in the lecture? 
# juvenile survival estimate is lowest, and has the highest elasticity and sensitivity
# this suggests improvements in juvenile survival can lead to a substantial increase in population growth

##QUESTIONb: Is this what you would expect based on the life-history of the species?
# The findings are consistent with life-history strategies where species with longer lifespans and lower reproductive rates (like orcas) typically have higher sensitivity to juvenile survival. 
# This suggests a life-history strategy focused on ensuring offspring survival to maintain population levels