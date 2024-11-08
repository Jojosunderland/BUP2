#download packages

my_packages <- c('dplyr', 'tidyr', 'marked', 'ggplot2', 'R2ucare')
new_packages <- my_packages[!(my_packages %in% installed.packages()[,'Package'])]
if(length(new_packages)) install.packages(new_packages) 

# load packages 
library(dplyr)
library(tidyr)
library(marked)
library(ggplot2)
library(R2ucare)

#load sparrow recapture dataset
longdata <- read.table("~/sparrowrecap.txt", header = TRUE, sep = '\t')
head(longdata)

# data exploration

# the number of unique individuals in the dataframe
length(unique(longdata$id)) 
# equal number of observations of males and females 
table(longdata$sex) 
# captures from 1998-2007
table(longdata$year) 
# at 4 different island locations
table(longdata$island) 

#marked function requires wide format so we are converting it here:

temp <- longdata[,1:2] # take the first two columns, id and year and put into a temporary dataframe
temp$detect <- 1 # add column for detection (all 1s because these represent captures) 

temp <- temp %>%
  # remove duplicates, which may occur when individuals are caught multiple times in an sampling event
  distinct() %>%
  # spread out data. The fill = 0 adds rows for combinations of id and year where individuals were not observed
  spread(year, detect, fill = 0) %>% 
  # for every individual....
  group_by(id) %>%
  # paste together 0's and 1's using unite()
  # here we are pasting the strings together from the second column (first capture event)
  # to the last capture event ("tail(names(.),1)")
  # use sep="" so there are no characters separating 0's and 1's
  unite("ch", 2:tail(names(.),1), sep = "")

sparrow <- as.data.frame(temp) # new dataframe called sparrow
head(sparrow)

##add back info on the individuals ID using match()

sparrow$island <- longdata$island[match(sparrow$id, longdata$id)] 
# this creates a new column called island in the sparrow df...
# using the entry from the island column in the longdata df... 
# where id in the sparrow df matches the id in the longdata df

sparrow$sex <- as.factor(longdata$sex[match(sparrow$id, longdata$id)])

sparrow <- droplevels(subset(sparrow, select = -id)) # remove id column so capture histories appear in first column
head(sparrow)

## Now, we are ready to run a cMR analysis

mod1 <- crm(sparrow) # capture-mark-recapture (cmr) model
mod1 # examine model and coefficient estimates

# phi = apparent survival
# p = detection probability

mod1 <- cjs.hessian(mod1) # refit model with precision estimates

# estimates of data are stored within the results section of the model under 'reals'
mod1$results$reals

#transform the data back to the data scale (from latent scale) using plogis()
plogis(mod1$results$beta$Phi)
plogis(mod1$results$beta$p)

predict(mod1, newdata=data.frame(sex = c('Female', 'Male')), se=T) 
# N.b. In this case, there are no groups or covariates in the model and so the 'newdata' argument is not used 

##THEREFORE, the pro of survivng between capture events was 0.518, 
##and the probability of detection an animal during a capture event was 0.580

## UNEQUAL SAMPLING INTERVALS

#the model assumes equal time between capture events. this assumption can be relaxed by including a vector of time intervals
mod2 <- crm(sparrow, time.intervals = c(1,2,1,1,1,1,1,3,4))
mod2$results$reals

#these models assume constant survival and detection which isn't very realistic. it's useful to include static covariates

## INCLUDING STATIC COVARIATES

#islands in the metapopulation are all different, and it is likely that probability rates differ between islands
# we can test this by addding some additional complexity to the model and allow the detection prob to vary between islands
#using the marked package

#step 1: data processing
# the outputs are lists not data frames. you can examine different components of the lists using double square brackets or $

sparrow.proc <- process.data(sparrow) # built in function for data processing
str(sparrow.proc)

head(sparrow.proc[[1]])

head(sparrow.proc$data)

#step 2: building the design matrix

sparrow.ddl <- make.design.data(sparrow.proc) # built in function for building design matrix 
str(sparrow.ddl)

head(sparrow.ddl[[1]])
head(sparrow.ddl$Phi)

# step 3: set up and execute candidate models
# specify model formulation: capture probability depends on island
p.island <- list(formula=~island) 

mod3 <- crm(sparrow.proc, 
            sparrow.ddl, 
            model.parameters = list(p = p.island), 
            accumulate=FALSE, hessian = TRUE)

mod3$results$reals


### Q: Yes detection probability varies among islands. Myken has the lowest, 
# this could be due to them also having the lowest survival
# to compare with the other more simpler model, use this code:
(mod3$results$AIC) # AIC is used for model comparison
(mod1$results$AIC)
# Lower AIC values indicate a better fitting model, so we accept the model with varying detection probabilities

## Test whether survival probabilities differ between the islands

sparrow.proc <- process.data(sparrow) 
sparrow.ddl <- make.design.data(sparrow.proc) 

Phi.island <- list(formula=~island) # survival probability depends on island
p.island <- list(formula=~island) # capture probability depends on island

mod4 <- crm(sparrow.proc, 
            sparrow.ddl, 
            model.parameters = list(Phi = Phi.island, 
                                    p = p.island), 
            accumulate=FALSE, hessian = TRUE)

mod4$results$reals

#compare models
(mod4$results$AIC)
(mod3$results$AIC)
# The AIC value of the new model is slightly higher, indicating that including this extra term does not improve model fit. 
# We therefore find no evidence that survival rates vary among islands.

## Does survival probability vary between sexes?

sparrow.proc <- process.data(sparrow)
sparrow.ddl <- make.design.data(sparrow.proc)

fit.models <- function() {
  Phi.dot <- list(formula=~1) # constant survival
  Phi.sex <- list(formula=~sex) # survival differs between sexes
  Phi.island <- list(formula=~island) # survival differs between islands
  Phi.sex.island <- list(formula=~sex+island) # survival differs between sexes and islands
  p.dot <- list(formula=~1) # constant detection
  p.sex <- list(formula=~sex) # detection probability differs between sexes
  p.island <- list(formula=~island) # detection probability differs between islands
  p.sex.island <- list(formula=~sex+island) # detection probability differs between sexes and islands
  cml <- create.model.list(c("Phi","p"))
  results <- crm.wrapper(cml, data=sparrow.proc, ddl=sparrow.ddl,
                         external=FALSE, accumulate=FALSE, hessian=TRUE)
  return(results)
}

sparrow.models <- fit.models() # run function 

sparrow.models # display model table

#plot the detection probabilities

mod5 <- sparrow.models[[2]]

ggplot(mod5$results$reals$p, aes(island, estimate, ymin=lcl, ymax=ucl)) + 
  geom_errorbar(width=0.2) + geom_point() + ylim(0,1)

#plot the sex and island differences

mod6 <- sparrow.models[[10]]
ggplot(mod6$results$reals$Phi, aes(sex, estimate, ymin=lcl, ymax=ucl)) + 
  geom_errorbar(width=0.2) + geom_point() + ylim(0,1)

mod7 <- sparrow.models[[6]]
ggplot(mod7$results$reals$Phi, aes(island, estimate, ymin=lcl, ymax=ucl)) + 
  geom_errorbar(width=0.2) + geom_point() + ylim(0,1)


## INCLUDING TIME-VARYING COVARIATES
# we cantest whether survival (or detection) probabilities differ depending on a factor that varies over time
#e.g. some aspect of weather

# test whether it being an extra cold winter affects survival
sparrow.ddl$Phi$cold <- "Cold" # new column 
sparrow.ddl$Phi$cold[sparrow.ddl$Phi$time==2 | sparrow.ddl$Phi$time==5 | sparrow.ddl$Phi$time==8] <- "VeryCold" # very cold winters between capture events 2 and 3, 5 and 6, and 8 and 9

head(sparrow.ddl$Phi)

Phi.cold <- list(formula=~cold) 
p.island <- list(formula=~island) 

mod8 <- crm(sparrow.proc, 
            sparrow.ddl, 
            model.parameters = list(Phi = Phi.cold, 
                                    p = p.island), 
            accumulate=FALSE, hessian = TRUE)

mod8$results$reals

(mod8$results$AIC)
(mod5$results$AIC)

ggplot(mod8$results$reals$Phi, aes(cold, estimate, ymin=lcl, ymax=ucl)) + 
  geom_errorbar(width=0.2) + geom_point() + ylim(0,1)

# the model including an effect of cold on survival probability has a higher AIC,
# so there is no evidence that survival varied depending on our made-up variable.


## TIME VARYING SURVIVAL AND RECAPTURE PROBABILITIES

# we can test whether survival and detection probabilities varied each year of the study period using the built-in time variable.

sparrow.proc <- process.data(sparrow)
sparrow.ddl <- make.design.data(sparrow.proc)

fit.models <- function() {
  Phi.dot <- list(formula=~1) # constant survival
  Phi.time <- list(formula=~time) # survival varies over time
  p.island <- list(formula=~island) # detection probability differs between islands
  p.time <- list(formula=~time) # detection probability varies over time
  p.island.time <- list(formula=~island+time) # detection probability varies over time
  cml <- create.model.list(c("Phi","p"))
  results <- crm.wrapper(cml, data=sparrow.proc, ddl=sparrow.ddl,
                         external=FALSE, accumulate=FALSE, hessian=TRUE)
  return(results)
}

sparrow.models <- fit.models() # run function 

sparrow.models # display model table

mod9 <- sparrow.models[[2]]

(mod9$results$AIC)

(mod5$results$AIC)

ggplot(mod9$results$reals$p, aes(time, estimate, ymin=lcl, ymax=ucl, col=island)) + 
  geom_errorbar(width=0) + geom_point() + ylim(0,1)

#It seems that detection probabilities varied over time.

## GOODNESS OF FIT TESTS
# To test whether our data violate any of the basic assumptions of the CJS model, we can use the package R2ucare
#  test whether our data violate the ‘equal detection assumption’ or the ‘equal survival assumption’.

sparrow.gof <- sparrow$ch %>%
  strsplit('') %>%
  sapply(`[`) %>%
  t() %>%
  unlist() %>%
  as.numeric %>%
  matrix(nrow = nrow(sparrow))

### R2ucare perfroms three tests, which are unhelpfully called Tests 1, 2 and 3.
#Test 1: the overall test. Overall, is there evidence that animals have equal detection probabilities and equal survival?
#Test 2: Does recapture depend on when an animal was first marked? (Tests the equal detection assumption)
#Test 3: Does marking affect survival? (Tests the equal survival assumption)

overall_CJS(sparrow.gof, rep(1,nrow(sparrow)))
##The p-value is not significant, meaning we fail to reject the null hypothesis. There is therefore no strong evidence for overall lack-of-fit.

#If we found evidence of lack-of-fit, we would want to delve deeper to find out what the problem was. First, we can look at the equal detection assumption using Test 2, which has two components:
#Test 2 CT: Is there a difference in p at t+1 between those captured and not captured at t (when animals are known to be alive because are captured later in the study)?
#Test 2 CL: Is there a difference in the expected time of next recapture between individuals captured and not captured at t when animals are known to be alive?

test2ct <- test2ct(sparrow.gof, rep(1,nrow(sparrow))) 
test2ct

test2cl <- test2cl(sparrow.gof, rep(1,nrow(sparrow)))
test2cl

#Again, we fail to reject the null hypothesis, so no evidence of a problem with the equal detection assumption.

##Finally, Test 3 which tests the equal survival assumption and also has two components:
#Test 3 SR: Do individuals with previous marks have different survival rates than first-time captures?
#Test 3 SM: For animals seen again, does when they are recaptured depend on whether they were marked on or before t?

test3sr <- test3sr(sparrow.gof, rep(1,nrow(sparrow)))
test3sr

test3sm <- test3sm(sparrow.gof, rep(1,nrow(sparrow)))
test3sm

##If your data violate the assumptions, you might need to increase the complexity of your model by adding age classes, 
# including a time-varying individual covariate, or building a multistate or multievent model.

