##-- Model: SIRV model --##
##

## Loading packages
library(tidyverse)
library(deSolve)



## sir() -- Function for numerical analysis of model 1

sir_model <- function(t,y,parms) {
    # The with() function gives access to the named values of parms within the
    # local environment created by the function
    with(c(as.list(y),parms),{
        N <- S + I + R + V
        dSdt <- -lambda*S - theta*S
        dIdt <- lambda*S - gamma*I
        dRdt <- gamma*I
        dVdt <- theta*S
        dc <- lambda*S
        # Note: Population size is constant, so don't need to specify dRdt
        return(list(c(dSdt,dIdt, dRdt, dVdt, dc)))
    })
}

sir_model2 <- function(t,y,parms) {
    # The with() function gives access to the named values of parms within the
    # local environment created by the function
    with(c(as.list(y),parms),{
        N <- S + I + R + V
        dSdt <- -beta*S*I/N - theta*S
        dIdt <- beta*S*I/N - gamma*I
        dRdt <- gamma*I
        dVdt <- theta*S
        dc <- beta*S*I/N
        # Note: Population size is constant, so don't need to specify dRdt
        return(list(c(dSdt,dIdt, dRdt, dVdt, dc)))
    })
}
time <- 0

N0 <- 92900000
Initial_I <- 100000
Initial_R <- 0
Initial_V <- 0.6*N0  #0.2*N0
S0 <- N0 - Initial_I - Initial_R - Initial_V

S0
N0
Initial_I
Initial_R
Initial_V


pop.SI <- c(S = S0,  # Initially 2% of the population is susceptible
            I = Initial_I,      # Suspected cases in Dec 2018
            V = Initial_V,   # 60% of the pop are vaccinated
            R = Initial_R,
            c = 0 )       


## The final input our function needs is a vector of parameter values, which we
## can create in the same way:

values <- c(beta = 26,                     # Transmission coefficient
            #lambda = 26*500/N0,
            gamma = 1/(14/30),       # 1 / infectious period = 1/(14/30) month            N = N0,            # population size (constant)
            theta = 0.6 * 0.85      # Coverage (60%) * Efficacy (85%)
            )
time.out <- seq(0, 24, by = 0.01)

#sir_example <- sir_model2(t = 0, y = pop.SI, parms = values)

sirv <- data.frame(lsoda(
    y = pop.SI,               # Initial conditions for population
    times = time.out,             # Timepoints for evaluation
    func = sir_model2,                   # Function to evaluate
    parms = values                # Vector of parameters
))

sirv

sirv_long <- sirv %>% pivot_longer(
    - time,
    
    names_to = "compartments",
    values_to = "values"
)

incidence_sim <- sirv_long %>% filter(compartments == "c") %>% 
    mutate(
        incidence = c(NA, diff(values))
    ) 
incidence_sim %>% 
    ggplot(aes(x=time, y= incidence, colour = compartments)) + # Time on the x axis, number infected (I) on the y axis
    geom_line()+
    xlab("Time in days")+
    ylab("Number infected")+
    ggtitle("Measles in New York") 



sirv_long %>% filter(compartments == "c") %>% 
    ggplot(aes(x=time, y= values, colour = compartments)) + # Time on the x axis, number infected (I) on the y axis
    geom_line()+
    xlab("Time in days")+
    ylab("Number infected")+
    ggtitle("Measles in New York")




ggplot(aes(x=time, y= values, colour = compartments), data = sirv_long) + # Time on the x axis, number infected (I) on the y axis
    geom_line()+
    xlab("Time in days")+
    ylab("Number infected")+
    ggtitle("Measles in New York")#+
    #xlim(c(0,400))


# Fitting the model

nllikelihood <- function(
        parms,
        obsDat,
         
) {
    sirv <- data.frame(lsoda(
        y = pop.SI,               # Initial conditions for population
        times = time.out,             # Timepoints for evaluation
        func = sir_model2,                   # Function to evaluate
        parms = values                # Vector of parameters
    ))
    
    incidence <- diff(sirv$c)
    obsDat_c <- obsDat[1:length(incidence)]
    # simulate an epidemic corresponding to parms, evaluated
    #simDat <- simEpidemic(tseq = obsDat$time, parms = parms)
    nlls <- -dpois(
        obsDat_c,
        lambda = incidence
        log = TRUE
    )
    return(sum(nlls))
}

optim.vals <- optim(
    par = init.pars,
    fn = nllikelihood,
    fixed.params = disease_params(),
    obsDat = myDat,
    control = list(trace = trace, maxit = 150),
    method = "SANN"
)



