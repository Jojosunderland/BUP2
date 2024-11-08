install.packages("geodata",dependencies=TRUE,repos="https://cloud.r-project.org")
library(geodata)
library(terra)

# data preparation
#download ring ouzel and climate data from Learn
avi_dat <- read.table('~/Documents/WorkingD/BUP/Species distribution modelling/Data_SwissBreedingBirds copy.csv', header=T, sep=',')
nrow(avi_dat)
summary(avi_dat)

# Subset the data to the columns we will be working with
ouzel_cols <- c('Turdus_torquatus', 'bio_5', 'bio_2', 'bio_14', 'blockCV_tile')

ouzel_df <- data.frame(avi_dat)[ouzel_cols]

summary(ouzel_df)

# Enter working directory for this data, download current and future climate data for switzerland
output_dir<-"~/Documents/WorkingD/BUP/Species distribution modelling"

bio_curr <-worldclim_country("Switzerland",version="2.1", var='bio', res=10, lon=5.5, lat=45.5, path=output_dir)[[c(2,5,14)]]

bio_fut <- cmip6_world(var = "bio", model = "CNRM-CM6-1-HR", ssp = "245", res = 10,  time = "2041-2060",  lon = c(5.96, 10.49),  lat = c(45.82, 47.81),path=output_dir)[[c(2,5,14)]]


# A spatial mask of Switzerland in Swiss coordinates
bg <- rast('/vsicurl/https://damariszurell.github.io/SDM-Intro/CH_mask.tif')

bio_curr <- terra::project(bio_curr, bg)
bio_fut <- terra::project(bio_fut, bg)
#we need to change the projection of our climate data to match that of the bg file.

bio_curr <- terra::resample(bio_curr, bg)
bio_fut <- terra::resample(bio_fut, bg)
#we then need to make the resolution equivalent to bg. 


bio_curr <- terra::mask(bio_curr, bg)
bio_fut <- terra::mask(bio_fut, bg)
#we then need to clip the extent to match an outline of Switzerland

names(bio_curr) <- c('bio_2', 'bio_5', 'bio_14')
names(bio_fut) <- c('bio_2', 'bio_5', 'bio_14')

# plot the climate variables and their projected change
plot(bio_curr) # current
plot(bio_fut) # future

## MODEL FITTING
# Can you code a binomial GLM with all three predictors fitted as linear and squared terms?

model <- glm( Turdus_torquatus ~ bio_2 + I(bio_2^2) + bio_5 + I(bio_5^2) + bio_14 + I(bio_14^2), family='binomial', data=ouzel_df)

summary(model) # bio 5 and 4 has a strong effect from the summary


## testing and critiquing the model
# part 1: partial effects

# understand how each climate variable affects the presence/absence of ring ouzels. to do this we will plot partial effects
#bio_2
bio_2<-seq(min(ouzel_df$bio_2),max(ouzel_df$bio_2),0.1)

newdata_bio2<-expand.grid(bio_2,mean(ouzel_df$bio_5),mean(ouzel_df$bio_14))
names(newdata_bio2)<-c("bio_2","bio_5","bio_14")
#To use predict we need to generate a new dataset. We will set the other two climate variables at their mean.

response<-predict(model,newdata=newdata_bio2,type="response")
#We've told the predict function to make a prediction on the response scale, ie. in terms of presence/absence

plot(newdata_bio2$bio_2,response,type="l")

#bio_5
bio_5<-seq(min(ouzel_df$bio_5),max(ouzel_df$bio_5),0.1)

newdata_bio5<-expand.grid(bio_5,mean(ouzel_df$bio_2),mean(ouzel_df$bio_14))
names(newdata_bio5)<-c("bio_5","bio_2","bio_14")
#To use predict we need to generate a new dataset. We will set the other two climate variables at their mean.

response<-predict(model,newdata=newdata_bio5,type="response")
#We've told the predict function to make a prediction on the response scale, ie. in terms of presence/absence

plot(newdata_bio5$bio_5,response,type="l")

#bio_14
bio_14<-seq(min(ouzel_df$bio_14),max(ouzel_df$bio_14),0.1)

newdata_bio14<-expand.grid(bio_14,mean(ouzel_df$bio_5),mean(ouzel_df$bio_2))
names(newdata_bio14)<-c("bio_14","bio_5","bio_2")
#To use predict we need to generate a new dataset. We will set the other two climate variables at their mean.

response<-predict(model,newdata=newdata_bio14,type="response")
#We've told the predict function to make a prediction on the response scale, ie. in terms of presence/absence

plot(newdata_bio14$bio_14,response,type="l")

#part 2: spatial cross validation

#first we will exclude spatial block 1.
training1<-subset(ouzel_df,blockCV_tile!=1)
#next we re-run the glm
model1<-glm( Turdus_torquatus ~ bio_2 + I(bio_2^2) + bio_5 + I(bio_5^2) + bio_14 + I(bio_14^2), family='binomial', data=training1)
#we will then subset the data for a testing dataset and we will see how well the glm fitted to the other data does in predicting presences in this testing block
testing<-subset(ouzel_df,blockCV_tile==1)
predicted<-predict(model1,testing,type="response")
#so here we have a prediction as a proportion rather than as a pres/abs.

#The next step is to take different threshold values (i.e probability values at which we count a species as present or absent)
thresholds<-seq(0,1,0.001)
tpr<-c()
fpr<-c()
#We will then use a loop to consider the true positive and false positive rate at each threshold value
for(x in 1:length(thresholds)){
  
  predicted.pres<-as.numeric(predicted>thresholds[x])
  
  present_correct<-length(which(predicted.pres*testing$Turdus_torquatus==1))
  present_incorrect<-length(which(testing$Turdus_torquatus-predicted.pres==1))
  tpr[x]<-present_correct/(present_correct+present_incorrect)
  
  absent_correct<-length(which(predicted.pres+testing$Turdus_torquatus==0))
  absent_incorrect<-length(which(testing$Turdus_torquatus-predicted.pres==-1))
  fpr[x]<-absent_incorrect/(absent_incorrect+absent_correct)
  
}
#When we've run that we can plot the receiver operating characteristic (ROC) curve
plot(fpr,tpr,xlab="false positive rate, 1-sensitivity",ylab="true positive rate, specificity",type="l")
abline(0,1,lty=2)

#Finally to calculate AUC, we can imagine lots of small rectangles. For each one we will calculate its areas and then the area under curve is the sum of those

sortedvals<-cbind(fpr,tpr)[order(fpr),]

AUC<-0
for(x in 2:length(sortedvals[,1])){
  AUC<-AUC+(sortedvals[x,1]-sortedvals[x-1,1])*sortedvals[x-1,2]}

AUC


# decide on optimal probability threshold
#Here we'll find a value that maximises true positive rate and minimises false positive rate
bestthreshold<-thresholds[which.max(tpr+(1-fpr))]

bestthreshold

## Prediction
# current distribution
bio_curr_df <- as.data.frame(bio_curr, xy = TRUE, na.rm = FALSE)

bio_curr_df$pred_glm<-predict(model,bio_curr_df,type="response")

# Make binary predictions:
bio_curr_df$bin_glm <- ifelse(bio_curr_df$pred_glm > bestthreshold, 1, 0)

#plot the current probability and predicted presence/absence.

r_pred_curr <- rast(as.matrix(bio_curr_df[,-c(3:5)]),type="xyz")
plot(r_pred_curr)

# future distribution
bio_fut_df <- as.data.frame(bio_fut, xy = TRUE, na.rm = FALSE)

bio_fut_df$pred_glm<-predict(model,bio_fut_df,type="response")

# Make binary predictions:
bio_fut_df$bin_glm <- ifelse(bio_fut_df$pred_glm > bestthreshold, 1, 0)

#plot the current probability and predicted presence/absence.

r_pred_fut <- rast(as.matrix(bio_fut_df[,-c(3:5)]),type="xyz")
plot(r_pred_fut)

