##-- Model: SIRV model for DRC measles (wave 2) --##
##
## Steps:
##   1. Load data (second epidemic wave only)
##   2. Define the model (S, I, R, V compartments)
##   3. Set population size and starting conditions
##   4. Define the likelihood (Poisson)
##   6. Fit beta with optim()
##   7. Plot observed vs fitted incidence
##   8. Compare vaccination scenarios

library(tidyverse)
library(deSolve)
library(readxl)
library(janitor)
library(lubridate)


# 1. Load data ---------------------------------------------------------------

path_data <- "/Users/rachidmuleia/Dropbox/INS/SACEMA/PMF/MMED-Measles/Data/cases.xlsx"

measles_df <- read_excel(path_data, sheet = "WEB") |>
  clean_names() |>
  filter(
    year > 2018,
    year < 2021,
    country == "Democratic Republic of the Congo"
  ) |>
  mutate(
    year = as.numeric(year),
    month = as.numeric(month),
    date = ymd(paste(year, month, "01")),
    incidence = as.numeric(measles_suspect)
  ) |>
  arrange(date)

# Keep only wave 2 (from October 2019)
fit_data <- measles_df |>
  filter(date >= as.Date("2019-06-01"))

fit_data |>glimpse()


# 2. Model equations ---------------------------------------------------------
# S = susceptible, I = infectious, R = immune (natural infection)
# V = immune (vaccination), c = cumulative new cases

sir_model2 <- function(t, y, parms) {
  with(c(as.list(y), parms), {
    N <- S + I + R + V
    dS <- -beta * S * I / N - theta * S
    dI <-  beta * S * I / N - gamma * I
    dR <-  gamma * I
    dV <-  theta * S
    dc <-  beta * S * I / N
    list(c(dS, dI, dR, dV, dc))
  })
}


# 3. Population and parameters -----------------------------------------------

N0 <- 84e6                          # DRC population
p_immune <- 0.60                    # 60% already immune at t = 0
Initial_R <- 0   # in R (prior infection)
Initial_V <- round(N0 * p_immune)                     # no one vaccinated yet at t = 0

# I at t = 0: scale first month cases by infectious period (14 days)
Initial_I <- max(fit_data$incidence[1] * (14 / 30), 1)

S0 <- N0 - Initial_I - Initial_R - Initial_V

pop.SI <- c(
  S = S0,
  I = Initial_I,
  R = 0,
  V = Initial_V,
  c = 0
)

# Fixed parameters (not estimated)
gamma <- 30 / 14                    # recovery rate (14-day infectious period, monthly step)
vaccination_coverage <- 0.60
vaccination_efficacy <- 0.85
theta <- -log(1 - vaccination_coverage * vaccination_efficacy) / 12  # campaign over 12 months

# Starting guess for beta (from first-month incidence)
#beta_start <- fit_data$incidence[1] * N0 / (S0 * Initial_I)
beta_start <- 3
values <- c(
  beta = beta_start,
  gamma = gamma,
  theta = theta
)

# Monthly time points: 15 months of data -> times 0, 1, ..., 15
obsDat <- fit_data$incidence
time.out <- seq(0, length(obsDat), by = 1)


# 4. Run model once with starting values (before fitting) --------------------

sirv <- data.frame(lsoda(
  y = pop.SI,
  times = time.out,
  func = sir_model2,
  parms = values
))

sirv


# 5. Likelihood function -----------------------------------------------------
# Compare observed monthly cases to model-predicted new cases (Poisson)

nllikelihood <- function(parms, obsDat) {
  sirv <- data.frame(lsoda(
    y = pop.SI,
    times = time.out,
    func = sir_model2,
    parms = parms
  ))

  incidence <- diff(sirv$c)
  obsDat_c <- obsDat[seq_along(incidence)]
  incidence <- pmax(incidence, 0.01)

  -sum(dpois(obsDat_c, lambda = incidence, log = TRUE))
}


# optim() varies one parameter at a time, so we wrap beta in a small function
sirv_nll_beta <- function(beta) {
  parms <- c(
    beta = unname(beta),
    gamma = unname(gamma),
    theta = unname(theta)
  )
  nllikelihood(parms, obsDat)
}


# 6. Fit beta with optim() ---------------------------------------------------

fit <- optim(
  par = beta_start,
  fn = sirv_nll_beta,
  method = "L-BFGS-B",
  lower = 0.1,
  upper = 80,
  hessian = TRUE,
  control = list(maxit = 500)
)

fit_beta <- unname(fit$par)
fit_parms <- c(
  beta = fit_beta,
  gamma = gamma,
  theta = theta
)

cat("\n--- Fit results (wave 2) ---\n")
cat("beta =", round(fit_beta, 4), "\n")
cat("R0   =", round(fit_beta / gamma, 2), "\n")
cat("NLL  =", round(fit$value, 1), "\n")
cat("convergence =", fit$convergence, "(0 means OK)\n")


# 7. Plot observed vs fitted -------------------------------------------------
#
# Fitting uses monthly time.out (step = 1 month).
# For a smooth line on the plot, run the model again on a finer time grid.
# Do NOT take diff() from the fine grid and put it in fit_data (length mismatch).

# Monthly model output (same times as the fit — for checking numbers)

sirv_monthly <- data.frame(lsoda(
  y = pop.SI,
  times = time.out,
  func = sir_model2,
  parms = fit_parms
))


fit_data <- fit_data |>
  mutate(
    estimated_incidence = pmax(diff(sirv_monthly$c), 0)
  )

# Fine time grid for a smooth curve on the plot only
time.plot <- seq(0, length(obsDat), by = 0.1)

sirv_fit <- data.frame(lsoda(
  y = pop.SI,
  times = time.plot,
  func = sir_model2,
  parms = fit_parms
))

# Instantaneous incidence = dc/dt at each time point (smooth line)
plot_curve <- sirv_fit |>
  mutate(
    date = min(fit_data$date) %m+% days(round(time * 30)),
    estimated_incidence = pmax(fit_beta * S * I / (S + I + R + V), 0)
  )

plot_fit_DRC<-ggplot(fit_data, aes(x = date)) +
  geom_col(aes(y = incidence), fill = "steelblue", alpha = 0.6, width = 20) +
geom_line(
    data = plot_curve,
    aes(x = date, y = estimated_incidence),
    color = "firebrick",
    linewidth = 1
  ) +
  labs(

    x = "Date",
    y = "Monthly cases"
  ) + 
  theme_minimal()

  ggsave(
  filename = "Model_fit_DRC.png",   # .png / .pdf / .svg / .tiff ...
  plot     = plot_fit_DRC,
  path     = "/Users/rachidmuleia/Dropbox/INS/SACEMA/PMF/MMED-Measles",
  width    = 20, height = 12, units = "cm",
  dpi      = 300,
  bg       = "white"
)


# 8. Vaccination scenarios ---------------------------------------------------
#
# Each scenario sets a different vaccination rate theta.
# theta = 0 means no vaccination during the outbreak.
# Otherwise: reach (coverage x efficacy) of susceptibles over "months" months.
#
#   theta = -log(1 - coverage x efficacy) / months

scenarios <- tribble(
  ~scenario,              ~coverage, ~efficacy, ~months,
  "No vaccination",       0,         0,         NA,
  "Baseline (60%, 12 mo)", 0.60,     0.85,      12,
  "Fast (60%, 6 mo)",     0.60,      0.85,      6,
  "Slow (60%, 18 mo)",    0.60,      0.85,      18,
  "Higher (80%, 12 mo)",  0.80,      0.85,      12,
  "Lower (40%, 12 mo)",   0.40,      0.85,      12,
  "Middle (50%,12 mo)",   0.5,       0.85,      12
)

scenario_theta <- function(coverage, efficacy, months) {
  if (coverage <= 0) {
    return(0)
  }
  -log(1 - coverage * efficacy) / months
}

scenario_results <- vector("list", nrow(scenarios))

for (i in seq_len(nrow(scenarios))) {
  s <- scenarios[i, ]
  theta_i <- scenario_theta(s$coverage, s$efficacy, s$months)

  objective <- function(beta) {
    parms <- c(
      beta = unname(beta),
      gamma = unname(gamma),
      theta = unname(theta_i)
    )
    nllikelihood(parms, obsDat)
  }

  fit_i <- optim(
    par = beta_start,
    fn = objective,
    method = "L-BFGS-B",
    lower = 0.1,
    upper = 80,
    control = list(maxit = 500)
  )

  parms_i <- c(beta = unname(fit_i$par), gamma = gamma, theta = theta_i)
  sirv_i <- data.frame(lsoda(
    y = pop.SI,
    times = time.out,
    func = sir_model2,
    parms = parms_i
  ))
  est_i <- pmax(diff(sirv_i$c), 0)

  scenario_results[[i]] <- tibble(
    scenario = s$scenario,
    theta = theta_i,
    beta = unname(fit_i$par),
    R0 = unname(fit_i$par) / gamma,
    NLL = fit_i$value,
    pearson_r = cor(obsDat, est_i),
    RMSE = sqrt(mean((obsDat - est_i)^2))
  )
}

scenario_summary <- bind_rows(scenario_results) |>
  arrange(NLL) |>
  mutate(
    delta_NLL = NLL - min(NLL),
    best = NLL == min(NLL)
  )

cat("\n--- Vaccination scenarios (lower NLL = better fit) ---\n")
print(scenario_summary, n = Inf)

# Bar chart: which scenario fits best?
ggplot(scenario_summary, aes(x = reorder(scenario, -NLL), y = NLL, fill = best)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "firebrick")) +
  labs(
    title = "Model fit by vaccination scenario",
    subtitle = "Negative log-likelihood (Poisson); lower is better",
    x = NULL,
    y = "NLL"
  ) +
  theme_minimal()

# Fitted curves for each scenario (monthly)
scenario_curves <- vector("list", nrow(scenarios))

for (i in seq_len(nrow(scenarios))) {
  s <- scenarios[i, ]
  row <- scenario_summary |> filter(scenario == s$scenario)
  parms_i <- c(beta = row$beta, gamma = gamma, theta = row$theta)
  sirv_i <- data.frame(lsoda(
    y = pop.SI,
    times = time.out,
    func = sir_model2,
    parms = parms_i
  ))

  scenario_curves[[i]] <- fit_data |>
    select(date, incidence) |>
    mutate(
      scenario = s$scenario,
      estimated_incidence = pmax(diff(sirv_i$c), 0)
    )
}

scenario_plot_data <- bind_rows(scenario_curves)

ggplot(scenario_plot_data, aes(x = date)) +
  geom_col(
    data = fit_data,
    aes(y = incidence),
    fill = "grey85",
    width = 20
  ) +
  geom_line(aes(y = estimated_incidence, color = scenario), linewidth = 0.9) +
  facet_wrap(~scenario, ncol = 2) +
  labs(
    title = "Observed vs fitted — vaccination scenarios",
    subtitle = "Grey bars = observed monthly cases",
    x = "Date",
    y = "Monthly cases",
    color = "Scenario"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# Save comparison table
output_dir <- "/Users/rachidmuleia/Dropbox/INS/SACEMA/PMF/MMED-Measles/Output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write_csv(scenario_summary, file.path(output_dir, "Model_R_vaccination_scenarios.csv"))
