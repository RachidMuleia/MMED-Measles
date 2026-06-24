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
Initial_I <- 1000
Initial_R <- 0
Initial_V <- 0.4*N0  #0.2*N0
S0 <- N0 - Initial_I - Initial_R - Initial_V
coverage <- 0.6
efficacy <- 0.85  
inf_d <- 14/30    # infectious period (in months)
trans_coef <- 6  #527*N0/(S0 * Initial_I)  # Transmission coefficient

    

S0
N0
Initial_I
Initial_R
Initial_V
trans_coef

pop.SI <- c(S = S0,  # Initially 2% of the population is susceptible
            I = Initial_I,      # Suspected cases in Dec 2018
            V = Initial_V,   # 60% of the pop are vaccinated
            R = Initial_R,
            c = 0 )       


## The final input our function needs is a vector of parameter values, which we
## can create in the same way:

values <- c(beta = trans_coef,                     
            gamma = 1/inf_d, 
            theta = - (log(1- coverage * efficacy))/12,
            N = N0                                    # population size (constant)
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
sirv_monthly <- data.frame(lsoda(
    y = pop.SI,               # Initial conditions for population
    times = seq(0, 19, by = 1),             # Timepoints for evaluation
    func = sir_model2,                   # Function to evaluate
    parms = values                # Vector of parameters
))
sirv_monthly %>% glimpse()

sirv_monthly_fit <- cases_19_8 %>% mutate(
    estimated_incidence = pmax(diff(sirv_monthly$c), 0) 
)

sirv_monthly_fit_plot <- sirv %>% mutate(
    date = min(sirv_monthly_fit$year_month) %m+%days(round(time*30)), 
    estimated_incidence = pmax(trans_coef*S*I/(S+I+R+V), 0)
)
sirv_monthly_fit %>% glimpse()
ggplot(sirv_monthly_fit, aes(x = year_month)) +
    geom_col(aes(y = measles_suspect)) +
    geom_line(data = sirv_monthly_fit_plot, aes(x = date, y = estimated_incidence))

trans_coef






sirv_long <- sirv %>% pivot_longer(
    - time,
    
    names_to = "compartments",
    values_to = "values"
)

incidence_sim <- sirv_long %>% filter(compartments == "c") %>% 
    mutate(
        incidence = c(NA, diff(values))
    ) 

sirv_est_inc <- sirv %>% mutate(
    estim_incidence <- trans_coef * S * I/(S + I + V + R)
) %>% rename(
    estim_incidence = `estim_incidence <- trans_coef * S * I/(S + I + V + R)`
)

sirv_est_inc %>% glimpse()
summary(sirv_est_inc)

sim_inc_plot <- incidence_sim %>% 
    ggplot(aes(x=time, y= incidence, colour = compartments)) + # Time on the x axis, number infected (I) on the y axis
    geom_line()+
    xlab("Time in months")+
    ylab("Number infected")+
    ggtitle("Measles in DRC: Incidence")+
    xlim(c(0,12)) 

print(cases_plot + sim_inc_plot)

ggplot() +
    aes(x = time) +
    
    geom_col(aes(x = month_2, y = measles_suspect), cases_19_8) +
    geom_line(aes(y = I, colour = "blue"), sirv_est_inc) + 
    xlim(c(0,20))
    #geom_errorbar(aes(ymin = lci, ymax = uci, color = "observed"), myDat) +
    labs(x = NULL, y = "prevalence") +
    theme_classic() +
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.5, 0.05),
        legend.justification.inside = c(0.5, 0),
        legend.direction = "horizontal"
    ) 


sirv_long %>% filter(compartments == "c") %>% 
    ggplot(aes(x=time, y= values, colour = compartments)) + # Time on the x axis, number infected (I) on the y axis
    geom_line()+
    xlab("Time in months")+
    ylab("Number infected")+
    ggtitle("Measles in DRC: Cumulative Incidence ") +
    xlim(c(0,12))




ggplot(aes(x=time, y= values, colour = compartments), data = sirv_long) + # Time on the x axis, number infected (I) on the y axis
    geom_line()+
    xlab("Time in months")+
    ylab("Number infected")+
    ggtitle("Measles in DRC")+
    xlim(c(0,12))


# Fitting the model

nllikelihood <- function(
        parms,
        obsDat
         
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
        lambda = incidence,
        log = TRUE
    )
    return(sum(nlls))
}


