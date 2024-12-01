## SIR Model for the Proof of concepts for my grant proposal

library(DSAIDE)
dsaidemenu()

vignette('DSAIDE')
help('simulate_Host_Heterogeneity_Model_ode')

library(deSolve)

# Parameters for the model
S1 <- 1 # Susceptible for host species 1
I1 <- 1 # Infected for host species 1
R1 <- 1 # Removed for host species 1

S2 <- 1 # Susceptible for host species 2
I2 <- 1 # Infected for host species 2
R2 <- 1 # Removed for host species 2

P <- 1 # Predator population

b1 <- 1 # Birth rate for host species 1
b2 <- 1 # Birth rate for host species 2
d1 <- 1 # Disease induced death rate for host species 1
d2 <- 1 # Disease induced death rate for host species 2

beta11 <- 1 # within species transmission for host species 1
beta12 <- 1 # between species transmission (2 -> 1)
beta22 <- 1 # within species transmission for host species 2
beta21 <- 1 # between species transmission (1 -> 2)

alpha1 <- 1 # predation rate on host species 1
alpha2 <- 1 # predation rate on host species 2

p <- 1 #conversion efficiency for turning prey into predator
delta <- 1 # death rate of the predator

## Anthropogenic impacts ##

gamma1 <- 1 # Impact of climate change on host species 1
gamma2 <- 1 # Impact of climate change on host species 2
theta1 <- 1 # Impact of land-use change on host species 1
theta2 <- 1 # Impact of land-use change on host species 2

b1 <- 1*gamma1*theta1 # Birth rate for host species 1
b2 <- 1*gamma2*theta2  # Birth rate for host species 2
d1 <- 1*gamma1*theta1  # Disease induced death rate for host species 1
d2 <- 1*gamma2*theta2 # Disease induced death rate for host species 2

beta11 <- 1*gamma1*theta1  # within species transmission for host species 1
beta12 <- 1*gamma2*theta2 # between species transmission (2 -> 1)
beta22 <- 1*gamma2*theta2 # within species transmission for host species 2
beta21 <- 1 *gamma1*theta1 # between species transmission (1 -> 2)


## Total population sizes for host species
N1 <- S1 + I1 + R1
N2 <- S2 + I2 + R2

## Model for host species 1
dS1 <- b1*N1 - d1*S1 - S1*beta11*I1 - S1*beta12*I2 - alpha1*P*S1
dI1 <- S1*beta11*I1 + S1*beta12*I2 - (alpha1*P - d1)*I1
dR1 <- (alpha1*P+d1)*I1

## Model for host species 2
dS2 <- b2*N2 - d2*S2 - S2*beta22*I2 - S2*beta21*I1 - alpha2*P*S2
dI2 <- S2*beta22*I2 + S2*beta21*I1 - (alpha2*P - d2)*I2
dR2 <- (alpha2*P+d2)*I2

# Model for predator population
dP <- p*alpha1*(S1+I1)*P + p*alpha2*(S2+I2)*P - delta*P

# Initial conditions for the model
initial_state <- c(S1 = 500, I1 = 10, R1 = 0, # Host species 1 (e.g., buffalo)
                   S2 = 300, I2 = 5, R2 = 0,  # Host species 2 (e.g., kudu)
                   P = 50)                    # Predator population (e.g., lions)

# Time sequence for the simulation
time <- seq(0, 365, by = 1) # Simulate for 1 year, daily steps

# Solve the differential equations
output <- ode(y = initial_state, times = time, func = multi_host_predator_model, parms = parameters)

# Convert output to a data frame
output <- as.data.frame(output)

# Plot the results
library(ggplot2)
ggplot(output, aes(x = time)) +
  geom_line(aes(y = S1, color = "S1 - Susceptible (Buffalo)")) +
  geom_line(aes(y = I1, color = "I1 - Infected (Buffalo)")) +
  geom_line(aes(y = S2, color = "S2 - Susceptible (Kudu)")) +
  geom_line(aes(y = I2, color = "I2 - Infected (Kudu)")) +
  geom_line(aes(y = P, color = "P - Predator (Lions)")) +
  labs(title = "Multi-Host SIR Model with Predator-Prey Dynamics",
       x = "Time (days)", y = "Population Size") +
  theme_bw() +
  scale_color_manual(values = c("blue", "red", "green", "purple", "orange"))
