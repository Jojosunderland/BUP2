## For this practical, we will implement a DISCRETE, mechanistic model for six distinct populations:
# lichen, plants, moose forage, caribou, moose and wolves

rm(list = ls())
graphics.off()
quartz()

library(ggplot2)


#############
##Model skeleton
#############

##set up parameter values
sP <- 0.5
Kp <- 240
lP <- 2.68275
hp <- 300
sL <- 0.07
Kl <- 870
lL <- 2.68275
hl <- 400
sH <- 0.5
Kh <- 970
lH <- 3.1682
hh <- 1000
fC <- 1
eC <- 1.85
mC <- 0 #hunting
fM <- 1.5
eM <- 0.6
mM <- 0
mW <- 0.1
b <- 0.8
g <- 0.38
dC <- 460
dM <- 46



nsteps <- 200
pop.df.0 <- data.frame(time=1:nsteps,P=numeric(nsteps),L=numeric(nsteps),H=numeric(nsteps),C=numeric(nsteps),M=numeric(nsteps),W=numeric(nsteps))

## initial population sizes
pop.df.0 <- within(pop.df.0,{
  P[1] <- 240
  L[1] <- 870
  H[1] <- 970
  C[1] <- 7
  M[1] <- 25
  W[1] <- 8
})

## the model
for(t in 2:nsteps){
  pop.df.0 <- within(pop.df.0,{
    ## Primary producers
    P.birth <- sP*P[t-1]*(1-P[t-1]/Kp)
    P.death <- (C[t-1]*lP*P[t-1])/(hp+P[t-1])
    P[t] <- max(0, P[t-1] + P.birth - P.death) ##plants consumed by Caribou
    #the max is used to avoid negative values
    
    L.birth <- sL*L[t-1]*(1-L[t-1]/Kl)
    L.death <- lL*L[t-1]*C[t-1]/(hl+L[t-1])
    L[t] <- max(0,L[t-1] + L.birth - L.death) ##lichen consumed by Caribou
    
    H.birth <- sH*H[t-1]*(1-H[t-1]/Kh)
    H.death <- lH*H[t-1]*M[t-1]/(hh+H[t-1])
    H[t] <- max(0,H[t-1] + H.birth - H.death) ##plants consumed by Moose
    
    ## First trophic level
    C.growth <- C[t-1]*fC*(P[t-1]/(P[t-1]+hp))*(L[t-1]/(L[t-1]+hl))
    C.death <- C[t-1]*(1-P[t-1]/(P[t-1]+hp))*(1-L[t-1]/(L[t-1]+hl))
    C.pred <-  W[t-1]*eC*C[t-1]/(C[t-1]+dC) 
    C[t] <- max(0,C[t-1] + C.growth - C.death - C.pred) # caribou
    
    M.growth <- M[t-1]*fM*H[t-1]/(H[t-1]+hh)
    M.death <- M[t-1]*(1-H[t-1]/(H[t-1]+hh))
    M.pred <-  W[t-1]*eM*M[t-1]/(M[t-1]+dM)
    M[t] <- max(0,M[t-1] + M.growth - M.death - M.pred) # moose
    
    ## Predator
    W.growth <- (b*(C[t-1]/(dC+C[t-1])) + b*(M[t-1]/(dM+M[t-1])))*W[t-1]
    W.death <- g*W[t-1]
    W[t] <- max(0,W[t-1] + W.growth - W.death * mW) #mW is the hunting component
     #wolves
  })
}

#plot it

colors <- c("Plants"="brown","Lichen"="orange","Shrubs"="green","Caribou"="purple","Moose"="blue","Wolf"="red")
ggplot(data = pop.df.0)+
  geom_line(mapping=aes(x=time,y=P,color="Plants")) +
  geom_line(mapping=aes(x=time,y=L,color="Lichen")) +
  geom_line(mapping=aes(x=time,y=H,color="Shrubs")) +
  geom_line(mapping=aes(x=time,y=C,color="Caribou")) +
  geom_line(mapping=aes(x=time,y=M,color="Moose")) +
  geom_line(mapping=aes(x=time,y=W,color="Wolf")) +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Time", y = "Densities", color="Legend")+
  scale_color_manual(values = colors)


## TASKS ##

## limiting wolf population to 10 or 30 wolves
# added this to the for loop: W[t] <- min(10, W[t]), min(30, W[t])

for(t in 2:nsteps){
  pop.df.0 <- within(pop.df.0,{
    ## Primary producers
    P.birth <- sP*P[t-1]*(1-P[t-1]/Kp)
    P.death <- (C[t-1]*lP*P[t-1])/(hp+P[t-1])
    P[t] <- max(0, P[t-1] + P.birth - P.death) ##plants consumed by Caribou
    #the max is used to avoid negative values
    
    L.birth <- sL*L[t-1]*(1-L[t-1]/Kl)
    L.death <- lL*L[t-1]*C[t-1]/(hl+L[t-1])
    L[t] <- max(0,L[t-1] + L.birth - L.death) ##lichen consumed by Caribou
    
    H.birth <- sH*H[t-1]*(1-H[t-1]/Kh)
    H.death <- lH*H[t-1]*M[t-1]/(hh+H[t-1])
    H[t] <- max(0,H[t-1] + H.birth - H.death) ##plants consumed by Moose
    
    ## First trophic level
    C.growth <- C[t-1]*fC*(P[t-1]/(P[t-1]+hp))*(L[t-1]/(L[t-1]+hl))
    C.death <- C[t-1]*(1-P[t-1]/(P[t-1]+hp))*(1-L[t-1]/(L[t-1]+hl))
    C.pred <-  W[t-1]*eC*C[t-1]/(C[t-1]+dC) 
    C[t] <- max(0,C[t-1] + C.growth - C.death - C.pred) # caribou
    
    M.growth <- M[t-1]*fM*H[t-1]/(H[t-1]+hh)
    M.death <- M[t-1]*(1-H[t-1]/(H[t-1]+hh))
    M.pred <-  W[t-1]*eM*M[t-1]/(M[t-1]+dM)
    M[t] <- max(0,M[t-1] + M.growth - M.death - M.pred) # moose
    
    ## Predator - wolves
    W.growth <- (b*(C[t-1]/(dC+C[t-1])) + b*(M[t-1]/(dM+M[t-1])))*W[t-1]
    W.death <- g*W[t-1]
    W[t] <- min(max(0,W[t-1] + W.growth - W.death),10) #pop limited to 0-10 wolves
    #wolves
  })
}


colors <- c("Plants"="brown","Lichen"="orange","Shrubs"="green","Caribou"="purple","Moose"="blue","Wolf"="red")
ggplot(data = pop.df.0)+
  geom_line(mapping=aes(x=time,y=P,color="Plants")) +
  geom_line(mapping=aes(x=time,y=L,color="Lichen")) +
  geom_line(mapping=aes(x=time,y=H,color="Shrubs")) +
  geom_line(mapping=aes(x=time,y=C,color="Caribou")) +
  geom_line(mapping=aes(x=time,y=M,color="Moose")) +
  geom_line(mapping=aes(x=time,y=W,color="Wolf")) +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Time", y = "Densities", color="Legend")+
  scale_color_manual(values = colors)

## Wolf hunting model ##
## Added in the hunting component, changed mW to 0.1

sP <- 0.5
Kp <- 240
lP <- 2.68275
hp <- 300
sL <- 0.07
Kl <- 870
lL <- 2.68275
hl <- 400
sH <- 0.5
Kh <- 970
lH <- 3.1682
hh <- 1000
fC <- 1
eC <- 1.85
mC <- 0 #hunting
fM <- 1.5
eM <- 0.6
mM <- 0
mW <- 0.1
b <- 0.8
g <- 0.38
dC <- 460
dM <- 46



nsteps <- 200
pop.df.0 <- data.frame(time=1:nsteps,P=numeric(nsteps),L=numeric(nsteps),H=numeric(nsteps),C=numeric(nsteps),M=numeric(nsteps),W=numeric(nsteps))

## initial population sizes
pop.df.0 <- within(pop.df.0,{
  P[1] <- 240
  L[1] <- 870
  H[1] <- 970
  C[1] <- 7
  M[1] <- 25
  W[1] <- 8
})

## the model
for(t in 2:nsteps){
  pop.df.0 <- within(pop.df.0,{
    ## Primary producers
    P.birth <- sP*P[t-1]*(1-P[t-1]/Kp)
    P.death <- (C[t-1]*lP*P[t-1])/(hp+P[t-1])
    P[t] <- max(0, P[t-1] + P.birth - P.death) ##plants consumed by Caribou
    #the max is used to avoid negative values
    
    L.birth <- sL*L[t-1]*(1-L[t-1]/Kl)
    L.death <- lL*L[t-1]*C[t-1]/(hl+L[t-1])
    L[t] <- max(0,L[t-1] + L.birth - L.death) ##lichen consumed by Caribou
    
    H.birth <- sH*H[t-1]*(1-H[t-1]/Kh)
    H.death <- lH*H[t-1]*M[t-1]/(hh+H[t-1])
    H[t] <- max(0,H[t-1] + H.birth - H.death) ##plants consumed by Moose
    
    ## First trophic level
    C.growth <- C[t-1]*fC*(P[t-1]/(P[t-1]+hp))*(L[t-1]/(L[t-1]+hl))
    C.death <- C[t-1]*(1-P[t-1]/(P[t-1]+hp))*(1-L[t-1]/(L[t-1]+hl))
    C.pred <-  W[t-1]*eC*C[t-1]/(C[t-1]+dC) 
    C[t] <- max(0,C[t-1] + C.growth - C.death - C.pred) # caribou
    
    M.growth <- M[t-1]*fM*H[t-1]/(H[t-1]+hh)
    M.death <- M[t-1]*(1-H[t-1]/(H[t-1]+hh))
    M.pred <-  W[t-1]*eM*M[t-1]/(M[t-1]+dM)
    M[t] <- max(0,M[t-1] + M.growth - M.death - M.pred) # moose
    
    ## Predator
    W.growth <- (b*(C[t-1]/(dC+C[t-1])) + b*(M[t-1]/(dM+M[t-1])))*W[t-1]
    W.death <- g*W[t-1]
    W[t] <- max(0,W[t-1] + W.growth - W.death - W[t-1]*mW) #mW is the hunting component
    #wolves
  })
}

colors <- c("Plants"="brown","Lichen"="orange","Shrubs"="green","Caribou"="purple","Moose"="blue","Wolf"="red")
ggplot(data = pop.df.0)+
  geom_line(mapping=aes(x=time,y=P,color="Plants")) +
  geom_line(mapping=aes(x=time,y=L,color="Lichen")) +
  geom_line(mapping=aes(x=time,y=H,color="Shrubs")) +
  geom_line(mapping=aes(x=time,y=C,color="Caribou")) +
  geom_line(mapping=aes(x=time,y=M,color="Moose")) +
  geom_line(mapping=aes(x=time,y=W,color="Wolf")) +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Time", y = "Densities", color="Legend")+
  scale_color_manual(values = colors)


## Moose hunting ##

## Added in the hunting component, changed mM to 0.1

sP <- 0.5
Kp <- 240
lP <- 2.68275
hp <- 300
sL <- 0.07
Kl <- 870
lL <- 2.68275
hl <- 400
sH <- 0.5
Kh <- 970
lH <- 3.1682
hh <- 1000
fC <- 1
eC <- 1.85
mC <- 0 #hunting
fM <- 1.5
eM <- 0.6
mM <- 0.1 # moose huting component added
mW <- 0 # wolf hunting back to 0
b <- 0.8
g <- 0.38
dC <- 460
dM <- 46



nsteps <- 200
pop.df.0 <- data.frame(time=1:nsteps,P=numeric(nsteps),L=numeric(nsteps),H=numeric(nsteps),C=numeric(nsteps),M=numeric(nsteps),W=numeric(nsteps))

## initial population sizes
pop.df.0 <- within(pop.df.0,{
  P[1] <- 240
  L[1] <- 870
  H[1] <- 970
  C[1] <- 7
  M[1] <- 25
  W[1] <- 8
})

## the model
for(t in 2:nsteps){
  pop.df.0 <- within(pop.df.0,{
    ## Primary producers
    P.birth <- sP*P[t-1]*(1-P[t-1]/Kp)
    P.death <- (C[t-1]*lP*P[t-1])/(hp+P[t-1])
    P[t] <- max(0, P[t-1] + P.birth - P.death) ##plants consumed by Caribou
    #the max is used to avoid negative values
    
    L.birth <- sL*L[t-1]*(1-L[t-1]/Kl)
    L.death <- lL*L[t-1]*C[t-1]/(hl+L[t-1])
    L[t] <- max(0,L[t-1] + L.birth - L.death) ##lichen consumed by Caribou
    
    H.birth <- sH*H[t-1]*(1-H[t-1]/Kh)
    H.death <- lH*H[t-1]*M[t-1]/(hh+H[t-1])
    H[t] <- max(0,H[t-1] + H.birth - H.death) ##plants consumed by Moose
    
    ## First trophic level
    C.growth <- C[t-1]*fC*(P[t-1]/(P[t-1]+hp))*(L[t-1]/(L[t-1]+hl))
    C.death <- C[t-1]*(1-P[t-1]/(P[t-1]+hp))*(1-L[t-1]/(L[t-1]+hl))
    C.pred <-  W[t-1]*eC*C[t-1]/(C[t-1]+dC) 
    C[t] <- max(0,C[t-1] + C.growth - C.death - C.pred) # caribou
    
    M.growth <- M[t-1]*fM*H[t-1]/(H[t-1]+hh)
    M.death <- M[t-1]*(1-H[t-1]/(H[t-1]+hh))
    M.pred <-  W[t-1]*eM*M[t-1]/(M[t-1]+dM)
    M[t] <- max(0,M[t-1] + M.growth - M.death - M.pred - M[t-1]*mM) # moose
    # added moose hunting
    
    ## Predator
    W.growth <- (b*(C[t-1]/(dC+C[t-1])) + b*(M[t-1]/(dM+M[t-1])))*W[t-1]
    W.death <- g*W[t-1]
    W[t] <- max(0,W[t-1] + W.growth - W.death) #no hunting
    #wolves
  })
}

colors <- c("Plants"="brown","Lichen"="orange","Shrubs"="green","Caribou"="purple","Moose"="blue","Wolf"="red")
ggplot(data = pop.df.0)+
  geom_line(mapping=aes(x=time,y=P,color="Plants")) +
  geom_line(mapping=aes(x=time,y=L,color="Lichen")) +
  geom_line(mapping=aes(x=time,y=H,color="Shrubs")) +
  geom_line(mapping=aes(x=time,y=C,color="Caribou")) +
  geom_line(mapping=aes(x=time,y=M,color="Moose")) +
  geom_line(mapping=aes(x=time,y=W,color="Wolf")) +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Time", y = "Densities", color="Legend")+
  scale_color_manual(values = colors)



