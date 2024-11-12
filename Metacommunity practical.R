## METACOMMUNITY MODEL ##
# Using the MCSim package

install.packages("devtools")

#install MCSim package

devtools::install_github('sokole/MCSim') 

# load packages

library(MCSim)
library(tidyverse)
library(fields)
library(vegan)
library(betapart)

#define a seed (to make simulations reproducible)
# a seed determines the starting point of random generators and allow to get the same results every time despite being 'random'
set.seed(1234) # you can remove this or set other seeds to get other results


## Model Initialisation ##

