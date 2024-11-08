# Introduction for modelling with R
# making population models with density-dependent logistic growth

# initialise a numeric vectir P of size 500, using numeric()

P <- numeric(length = 500)

# set the first value of P to 10 individuals, this is P at t=0
P[1] <- 10

# set parameters r and K
r <- 1.1
K <- 100

#compute P[2] at population at t = 1 from P[1]
P[2] <- r*(1-P[1]/K)*P[1] + P[1]

# compute P[3] at time =2 from P[2]
P[3] <- r *(1-P[2]/K)*P[2] + P[2]

# P[4], t=3
P[4] <- r*(1-P[3]/K)*P[3] + P[3]

# to simulate 500 time steps we will use a "for" loop
# set a loop for 499 time steps (2 to 500)

for(i in 2:500) {
  P[i] <- r*(1-P[i-1]/K)*P[i-1] + P[i-1]
}

# plot the output

plot(1:500, P, xlim=c(1,50), type = "l",
     xlab = "Time", ylab = "Population Size")

# plot P(t+1) vs P(t)

Pt <- seq(0,150,0.1)
Ptt <- r*(1-Pt/K)*Pt +Pt

plot(Pt, Ptt, col = "lightblue")

?lines()

lines(Pt, Ptt, type = "l", col = "seagreen2")

Pt <- seq(0,150,0.1) # that creates a vector from 0 to 150 with an increment of 0.1
Ptt <- r*(1-Pt/K)*Pt+Pt
lines(Pt, Ptt, col = "blue")
lines(c(1,150), c(1,150), col = "red") # that draws the diagonal, as we saw in the lecture, to show the equilibrium

## continuous model
install.packages("deSolve")
install.packages("ggplot2")
install.packages("tidyverse")
  
library(deSolve)
?ode()

# need to define parameters for ode() function

state <- c(P=10)
times <- seq(0,500, by=0.01)
parameters <- c(r=0.1, K=1000)

LG <- function(t,state,parameters){ ##logistic grown function, that takes a set of parameter values, initial conditions and a time sequence
  with(as.list(c(state, parameters)),{ ##"with" is a function that allows us to use the variable names directly - it looks for r, K and P in state and parameters
    
    dP <- r*(1-P/K)*P ##this is our logistic equation governing the rate of change of P
    
    return(list(dP)) ## return the rate of change - it needs to be a list
  }) # end with(as.list ...
}

out <- ode(y=state, times = times, func = LG, parms = parameters)
out

out.df <- data.frame(out)

plot(out.df, type = "l")

library(ggplot2)
ggplot(data = out.df)+
  geom_line(mapping=aes(x=time,y=P),color="seagreen2") +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Time", y = "P")


# model calibration using real data

#base code without the parameters
tmax <- 300
x <- numeric(tmax+1)
for(i in 2:(tmax+1)){
  x[i] <- r*x[i-1]*(1-x[i-1]/K)+x[i-1]+rnorm(1,0,x[i-1]/100)
}

# get the data

pop <- read.csv("pop_LG_simul_noise_small2.csv")

#plot

ggplot()+
  geom_line(mapping=aes(x=pop$time,y=pop$P))




