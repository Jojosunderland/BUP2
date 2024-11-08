## Predator-Prey interactions
install.packages("deSolve")
library(deSolve)
library(ggplot2)

#LV logistic growth function

LG <- function(t,state,parameters){ ##logistic grown function, that takes a set of parameter values, initial conditions and a time sequence
  with(as.list(c(state, parameters)),{ ##"with" is a function that allows us to use the variable names directly - it looks for r, K and P in state and parameters
    dP <- r*(1-P/K)*P ##this is our logistic equation governing the rate of change of P
    return(list(dP)) ## return the rate of change - it needs to be a list
  }) # end with(as.list ...
}

# LV Predator-Prey model

LV <- function(t,state,parameters){ ##lotka voltera function function, that takes a set of parameter values, initial conditions and a time sequence
  with(as.list(c(state, parameters)),{ ##"with" is a function that allows us to use the variable names directly - it looks for r, K and P in state and parameters
    dx <- alpha * x - beta * x* y ## prey population equation 
    dy <- delta * x * y - gamma * y ## predator population
    return(list(c(dx, dy))) ## return the rate of change - it needs to be a list
  }) # end with(as.list ...
}

state <- c(x=10, y = 10) ## the initial population values of prey and predator populations
parameters <- c(alpha = 0.1, beta = 0.02, delta = 0.02, gamma = 0.4) ## the equation parameters
times <-  seq(from = 1, to = 500, by = 0.1) ##a sequence of time steps – uses function seq()
out <- ode(y= state, times = times, func = LV, parms = parameters)
out.df <- data.frame(out)

#plot the output, this gives us the change in the two populations through time
ggplot(data = out.df)+
  geom_line(mapping=aes(x=time,y=x),color="blue") +
  geom_line(mapping=aes(x=time,y=y),color="red") +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Time", y = "P")

# plot the change in populations through phase space 
# plot one predator population as a function of the prey population
ggplot(data = out.df)+
  geom_path(mapping=aes(x=x,y=y),color="red") +
  xlim(0,70) +
  ylim(0,40) +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Prey", y = "Predator")

## TASK: Change the parameters and initial conditions
## Q: What happens when you increase/decrease them?
## A: When you increase initial pop size, the circle increases, when you increase the parameters it decreases


## Prey growth rate: exponential vs logistic
# in the original LV model, prey growth is exponential, lets make it logistic

LV_log <- function(t,state,parameters){ ##lotka voltera function function, that takes a set of parameter values, initial conditions and a time sequence
  with(as.list(c(state, parameters)),{ ##"with" is a function that allows us to use the variable names directly - it looks for r, K and P in state and parameters
    dx <- alpha * (1 - x/K) - beta * x* y ## prey population equation 
    dy <- delta * x * y - gamma * y ## predator population
    return(list(c(dx, dy))) ## return the rate of change - it needs to be a list
  }) # end with(as.list ...
}

# set K to 30
state <- c(x=10, y = 10) ## the initial population values of prey and predator populations
parameters <- c(alpha = 0.1, beta = 0.02, delta = 0.02, gamma = 0.4, K = 30) ## the equation parameters
times <-  seq(from = 1, to = 500, by = 0.1) ##a sequence of time steps – uses function seq()
out <- ode(y= state, times = times, func = LV_log, parms = parameters)
out.df <- data.frame(out)

## Q: How does it change  the results for the same set of parameter values listed above? Discuss (i.e. try to figure out why)



## Incorporating function response
# The rate at which predators can consume preys is called a functional response

# Type 2 functional response
x <- seq(0,50,0.1)
A <- 0.005  ## A is the hunting efficiency of the predator (high A = poor hunter)
y <- x/(1+A*x)
ggplot()+
  geom_line(mapping=aes(x=x,y=x/(1+A*x)),color="blue") +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Prey population", y = "Prey consumed")

# if predator is a bad hunter, prey reach carrying capacity much quicker as their not being killed too often

## Type 2 Functional response LV model

LV_FR2 <- function(t,state,parameters){ ##lotka voltera function function, that takes a set of parameter values, initial conditions and a time sequence
  with(as.list(c(state, parameters)),{ ##"with" is a function that allows us to use the variable names directly - it looks for r, K and P in state and parameters
    dx <- alpha * x - (beta * x* y)/(1+A*x) ## prey population equation 
    dy <- (delta * x * y)/(1+A*x) - gamma * y ## predator population
    return(list(c(dx, dy))) ## return the rate of change - it needs to be a list
  }) # end with(as.list ...
}

# Add A to parameters
state <- c(x=10, y = 10) ## the initial population values of prey and predator populations
parameters <- c(alpha = 0.1, beta = 0.02, delta = 0.02, gamma = 0.4, A = 0.02) ## the equation parameters
times <-  seq(from = 1, to = 500, by = 0.1) ##a sequence of time steps – uses function seq()
out <- ode(y= state, times = times, func = LV_FR2, parms = parameters)
out.df <- data.frame(out)

ggplot(data = out.df)+
  geom_line(mapping=aes(x=time,y=x),color="blue") +
  geom_line(mapping=aes(x=time,y=y),color="red") +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Time", y = "P")

# plot the change in populations through phase space 
# plot one predator population as a function of the prey population
ggplot(data = out.df)+
  geom_path(mapping=aes(x=x,y=y),color="red") +
  xlim(0,270) +
  ylim(0,150) +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Prey", y = "Predator")

## LV model with logistic growth and functional response type 2
x <- seq(0,30,0.1)
A <- 0.1
y <- x/(1+A*x)
ggplot()+
  geom_line(mapping=aes(x=x,y=y),color="blue") +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Prey population", y = "Prey consumed")

LV_FR2_log <- function(t,state,parameters){ ##lotka voltera function function, that takes a set of parameter values, initial conditions and a time sequence
  with(as.list(c(state, parameters)),{ ##"with" is a function that allows us to use the variable names directly - it looks for r, K and P in state and parameters
    dx <- alpha * (1-x/K) - (beta * x* y)/(1+A*x) ## prey population equation 
    dy <- (delta * x * y)/(1+A*x) - gamma * y ## predator population
    return(list(c(dx, dy))) ## return the rate of change - it needs to be a list
  }) # end with(as.list ...
}

state <- c(x=10, y = 10) ## the initial population values of prey and predator populations
parameters <- c(alpha = 0.1, beta = 0.02, delta = 0.02, gamma = 0.4, A = 0.01, K=30) ## the equation parameters
times <-  seq(from = 1, to = 500, by = 0.1) ##a sequence of time steps – uses function seq()
out <- ode(y= state, times = times, func = LV_FR2_log, parms = parameters)
out.df <- data.frame(out)

ggplot(data = out.df)+
  geom_line(mapping=aes(x=time,y=x),color="blue") +
  geom_line(mapping=aes(x=time,y=y),color="red") +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Time", y = "P")

# plot the change in populations through phase space 
# plot one predator population as a function of the prey population
ggplot(data = out.df)+
  geom_path(mapping=aes(x=x,y=y),color="red") +
  xlim(0,70) +
  ylim(0,40) +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Prey", y = "Predator")



################
# Competition
###############

## Three-species competition LV model: Limiting similarity

# 2 species:

parameters <- c(r=0.3, a12=1, a21=0.9, K = 100)
state <- c(x1 = 50, x2 = 10)

LS <- function(t,state,parameters){ ##lotka voltera function function, that takes a set of parameter values, initial conditions and a time sequence
  with(as.list(c(state, parameters)),{ ##"with" is a function that allows us to use the variable names directly - it looks for r, K and P in state and parameters
    dx1 <- r*x1*(1-(x1+a12*x2)/K) ## population 1 equation 
    dx2 <- r*x2*(1-(x2+a21*x1)/K)## population 2 population
    return(list(c(dx1, dx2))) ## return the rate of change - it needs to be a list
  }) # end with(as.list ...
}

times <- seq(0,1000,by=0.01)

out <- ode(y=state, times = times, func = LS, parms = parameters)

out.df <- data.frame(out)

ggplot(data = out.df)+
  geom_line(mapping=aes(x=time,y=x1),color="blue") +
  geom_line(mapping=aes(x=time,y=x2),color="red") +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Time", y = "P")


ggplot(data = out.df)+
  geom_path(mapping=aes(x=x1,y=x2),color="red") +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Species 1", y = "Species 2")


# add a third species:

alpha.func <- function(mu1,sig1,mu2,sig2,K1,K2,start,end){ ##this is the function to compute the alpha coefficients from the mean and standard deviations of the Gaussian niches of the species and the start and end values of the environment
  niche1 <- K1*dnorm(seq(start,end,length.out=100),mean=mu1,sd=sig1) ##dnorm() generates the values of the Gaussian. Check ?dnorm
  niche2 <- K2*dnorm(seq(start,end,length.out=100),mean=mu2,sd=sig2)
  a <- sum(niche1*niche2)/sum(niche1*niche1) ##because we have discrete values, we use a sum to approximate the integral
  return(a)
}

##Let's try different parameter values
D <- 10 ##distance between the niche optima
mu1 <- 5 ##niche optima of species 1
mu2 <- mu1+D ##niche optima of species 2
mu3 <- mu1+2*D ##niche optima of species 3
sig1 <- sig2 <- sig3 <- 10 ##all species niches have the same standard deviation for simplicity
start <- 0
end <- 30
K1 <- 200 ##carrying capacity species 1 and 3
K2 <- 250 ##carrying capacity species 2
a12 <- alpha.func(mu1,sig1,mu2,sig2,K1,K2,start,end)
a13 <- alpha.func(mu1,sig1,mu3,sig3,K1,K1,start,end)
a21 <- alpha.func(mu2,sig2,mu1,sig1,K2,K1,start,end)
a23 <- alpha.func(mu2,sig2,mu3,sig3,K2,K1,start,end)
a31 <- alpha.func(mu3,sig3,mu1,sig1,K1,K1,start,end)
a32 <- alpha.func(mu3,sig3,mu2,sig2,K1,K2,start,end)


##visualise the niches
resource <- seq(start,end,length.out=100)
niche1 <- dnorm(resource,mean=mu1,sd=sig1)*K1
niche2 <- dnorm(resource,mean=mu2,sd=sig2)*K2
niche3 <- dnorm(resource,mean=mu3,sd=sig3)*K1
ggplot()+
  geom_line(mapping=aes(x=resource,y=niche1),color="blue")+
  geom_line(mapping=aes(x=resource,y=niche2),color="red")+
  geom_line(mapping=aes(x=resource,y=niche3),color="darkgreen")

##setup and solve the system of differential equations
parameters <- c(a12=a12, a13=a13, a21=a21, a23=a23, a31=a31, a32=a32, r=0.3, K1 = K1, K2 = K2)
state <- c(x1=10, x2=10, x3=10)

LS_3 <- function(t,state,parameters){ ##lotka voltera function function, that takes a set of parameter values, initial conditions and a time sequence
  with(as.list(c(state, parameters)),{ ##"with" is a function that allows us to use the variable names directly - it looks for r, K and P in state and parameters
    dx1 <- r*x1*(1-(x1+a12*x2+a13*x3)/K1) ## population 1 equation 
    dx2 <- r*x2*(1-(x2+a21*x1+a23*x3)/K2)## population 2 equation
    dx3 <- r*x3*(1-(x3+a31*x1+a32*x2)/K1) ## population 3 equation
    return(list(c(dx1, dx2, dx3))) ## return the rate of change - it needs to be a list
  }) # end with(as.list ...
}

times <- seq(0,200,by=0.01)

out <- ode(y=state, times = times, func = LS_3, parms = parameters)

out.df <- data.frame(out)

#plot the populations
ggplot(data = out.df)+
  geom_line(mapping=aes(x=time,y=x1),color="blue") +
  geom_line(mapping=aes(x=time,y=x2),color="red") +
  geom_line(mapping=aes(x=time,y=x3), colour = 'darkgreen') +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Time", y = "P")


##Q: with these parameters which species survive?
##A: species 2 survives, species 1 and 3 go extinct (they're overlapped as they have the same K)

# if you change the distance and width of each niche of the species you'll get different outcomes
# the bigger w and the smaller d the more overlap with other niches and more competition



