##-- Model: SEIRV model for DRC measles (wave 2) --##
##
## Same as Model_Measles_script.R but with an added E (Exposed) compartment:
## people who are infected but not yet infectious (measles latent period).
##
##   S -> E -> I -> R   (+ S -> V by vaccination)
##
## Fixed:    gamma (recovery), sigma (1/latent period)
## Estimated: beta and theta (jointly, by maximum likelihood)

library(tidyverse)
library(deSolve)
library(readxl)
library(janitor)
library(lubridate)


# 1. Load data ---------------------------------------------------------------

path_data <- "/Users/rachidmuleia/Dropbox/INS/SACEMA/PMF/MMED-Measles/Data/cases.xlsx"

measles_df <- read_excel(path_data, sheet = "WEB") |>
  clean_names() |>
  filter(year > 2018, year < 2021,
         country == "Democratic Republic of the Congo") |>
  mutate(
    year = as.numeric(year),
    month = as.numeric(month),
    date = ymd(paste(year, month, "01")),
    incidence = as.numeric(measles_suspect)
  ) |>
  arrange(date)

# Keep only wave 2
fit_data <- measles_df |> filter(date >= as.Date("2019-06-01"))
fit_data |> glimpse()


# 2. Model equations (SEIRV) -------------------------------------------------
# S susceptible, E exposed (infected, not yet infectious), I infectious,
# R immune (natural), V immune (vaccination), c cumulative new cases.
# New cases are counted at symptom onset = flow out of E (sigma * E).

seirv_model <- function(t, y, parms) {
  with(c(as.list(y), parms), {
    N  <- S + E + I + R + V
    dS <- -beta * S * I / N - theta * S
    dE <-  beta * S * I / N - sigma * E
    dI <-  sigma * E - gamma * I
    dR <-  gamma * I
    dV <-  theta * S
    dc <-  sigma * E
    list(c(dS, dE, dI, dR, dV, dc))
  })
}


# 3. Population and parameters -----------------------------------------------

N0 <- 84e6                          # DRC population
p_immune <- 0.60                    # 60% already immune at t = 0
Initial_R <- 0
Initial_V <- round(N0 * p_immune)

# E and I at t = 0: scale first-month cases by latent (8 d) / infectious (14 d) periods
latent_days     <- 8
infectious_days <- 14
Initial_E <- max(fit_data$incidence[1] * (latent_days / 30), 1)
Initial_I <- max(fit_data$incidence[1] * (infectious_days / 30), 1)

S0 <- N0 - Initial_E - Initial_I - Initial_R - Initial_V

pop.SEI <- c(
  S = S0,
  E = Initial_E,
  I = Initial_I,
  R = Initial_R,
  V = Initial_V,
  c = 0
)

# Fixed parameters (not estimated)
gamma <- 30 / infectious_days        # recovery rate (monthly step)
sigma <- 30 / latent_days            # 1 / latent period (monthly step)

vaccination_coverage <- 0.60
vaccination_efficacy <- 0.85
theta <- -log(1 - vaccination_coverage * vaccination_efficacy) / 12

obsDat   <- fit_data$incidence
time.out <- seq(0, length(obsDat), by = 1)
time.plot <- seq(0, length(obsDat), by = 0.1)


# 4. Likelihood (Poisson) ----------------------------------------------------

nllikelihood <- function(parms, obsDat) {
  out <- data.frame(lsoda(
    y = pop.SEI, times = time.out, func = seirv_model, parms = parms
  ))
  incidence <- diff(out$c)
  obsDat_c  <- obsDat[seq_along(incidence)]
  incidence <- pmax(incidence, 0.01)
  -sum(dpois(obsDat_c, lambda = incidence, log = TRUE))
}


# 5. Fit beta and theta jointly (gamma, sigma fixed) -------------------------

objective <- function(par) {
  parms <- c(beta = par[1], gamma = gamma, sigma = sigma, theta = par[2])
  nllikelihood(parms, obsDat)
}

beta_start  <- 3
theta_start <- 0.059

fit <- optim(
  par     = c(beta_start, theta_start),
  fn      = objective,
  method  = "L-BFGS-B",
  lower   = c(0.1, 0),
  upper   = c(80, 1),
  hessian = TRUE,
  control = list(maxit = 500)
)

fit_beta  <- fit$par[1]
fit_theta <- fit$par[2]
fit_parms <- c(beta = fit_beta, gamma = gamma, sigma = sigma, theta = fit_theta)

cat("\n--- SEIRV fit results: beta & theta (wave 2) ---\n")
cat("beta  =", round(fit_beta, 4), "\n")
cat("theta =", round(fit_theta, 4), "\n")
cat("sigma =", round(sigma, 4), "(fixed, latent", latent_days, "d)\n")
cat("gamma =", round(gamma, 4), "(fixed, infectious", infectious_days, "d)\n")
cat("R0    =", round(fit_beta / gamma, 2), "\n")
cat("NLL   =", round(fit$value, 1), "\n")
cat("convergence =", fit$convergence, "(0 means OK)\n")


# 6. Plot observed vs fitted -------------------------------------------------
# Smooth fitted curve on the fine grid: incidence = sigma * E (onset rate).

plot_curve <- data.frame(lsoda(
  y = pop.SEI, times = time.plot, func = seirv_model, parms = fit_parms
)) |>
  mutate(
    date = min(fit_data$date) %m+% days(round(time * 30)),
    estimated_incidence = pmax(sigma * E, 0)
  )

plot_fit_SEIRV <- ggplot(fit_data, aes(x = date)) +
  geom_col(aes(y = incidence), fill = "steelblue", alpha = 0.6, width = 20) +
  geom_line(data = plot_curve, aes(x = date, y = estimated_incidence),
            color = "firebrick", linewidth = 1) +
  labs(title = "DRC measles wave 2 — SEIRV fit (beta & theta)",
       x = "Date", y = "Monthly cases") +
  theme_minimal()

ggsave(
  filename = "Model_fit_SEIRV.png",
  plot     = plot_fit_SEIRV,
  path     = "/Users/rachidmuleia/Dropbox/INS/SACEMA/PMF/MMED-Measles",
  width    = 20, height = 12, units = "cm", dpi = 300, bg = "white"
)

cat("\nSaved plot: Model_fit_SEIRV.png\n")
