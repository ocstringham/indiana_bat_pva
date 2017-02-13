#### try SIR model ####

## Load deSolve package
library(deSolve)

## Create an SIR function
sir <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    dS <- ( -beta * ((S * I)/N) ) + bu*S - du*S
    ##assume pups are born with WNS
    dI <-  beta * ((S * I)/N) + be*I - de*I
    
    return(list(c(dS, dI)))
  })
}

### Set parameters
## Proportion in each compartment: Susceptible 0.999999, Infected 0.000001, Recovered 0
init       <- c(S = 409*0.9, I = 40*0.1, N=409)
## beta: infection parameter; gamma: recovery parameter
parameters <- c(beta = 0.121, bu = 0.69, du = 1 - s1, be = 0.85*(lam_e/1.85), de = 1 - (lam_e/1.85) )
## Time frame
times      <- seq(0, 10, by = 1)

## Solve using ode (General Solver for Ordinary Differential Equations)
out <- ode(y = init, times = times, func = sir, parms = parameters)
## change to data frame
out <- as.data.frame(out)
## Delete time variable
out$time <- NULL
## Show data
head(out, 10)

## Plot
matplot(x = times, y = out, type = "l",
        xlab = "Time", ylab = "Susceptible and Recovered", main = "SIR Model",
        lwd = 1, lty = 1, bty = "l", col = 2:4)

## Add legend
legend(40, 0.7, c("Susceptible", "Infected", "Recovered"), pch = 1, col = 2:4, bty = "n")
