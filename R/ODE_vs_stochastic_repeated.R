# Load necessary libraries
library(bayesSSM)
library(ggplot2)
library(tidyr)
library(extraDistr)
library(rstan)
library(msm)

# Set up Stan options (optional but recommended)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# ---------- Functions and Fixed Settings ----------

# --- Simulation settings and true parameters ---
n_total <- 50 # Total population size
init_infected <- 3 # Initially infectious individuals
init_state <- c(n_total - init_infected, init_infected) # (s, i) at time 0
t_max <- 15 # Total number of days to simulate
true_lambda <- 0.5 # True infection parameter
true_gamma <- 0.2 # True removal parameter

# --- Gillespie simulation functions ---

epidemic_step <- function(state, t_end, lambda, gamma, n_total) {
  t <- 0
  s <- state[1]
  i <- state[2]
  while (t < t_end && i > 0) {
    rate_infection <- (lambda / n_total) * s * i
    rate_removal <- gamma * i
    rate_total <- rate_infection + rate_removal
    if (rate_total <= 0) break
    dt <- rexp(1, rate_total)
    if (t + dt > t_end) break
    t <- t + dt
    # Decide which event occurs:
    if (runif(1) < rate_infection / rate_total) {
      s <- s - 1
      i <- i + 1
    } else {
      i <- i - 1
    }
  }
  c(s, i)
}

simulate_epidemic <- function(
    n_total, init_infected, lambda, gamma, t_max
) {
  states <- matrix(0, nrow = t_max + 1, ncol = 2)
  states[1, ] <- c(n_total - init_infected, init_infected)
  state <- states[1, ]
  for (t in 1:t_max) {
    state <- epidemic_step(state, 1, lambda, gamma, n_total)
    states[t + 1, ] <- state
  }
  states
}

# --- PMCMC functions for the epidemic model ---
# Prior definitions
log_prior_lambda <- function(lambda) {
  extraDistr::dhnorm(lambda, sigma = 1, log = TRUE)
}

log_prior_gamma <- function(gamma) {
  extraDistr::dhnorm(gamma, sigma = 1, log = TRUE)
}

log_priors <- list(
  lambda = log_prior_lambda,
  gamma = log_prior_gamma
)

init_fn_epidemic <- function(particles) {
  # Return a matrix with particles rows; each row is the initial state (s, i)
  matrix(rep(init_state, each = particles), nrow = particles, byrow = FALSE)
}

transition_fn_epidemic <- function(particles, lambda, gamma, t) {
  # t=1 deterministic
  if (t == 1) {
    return(particles)
  }

  new_particles <- t(apply(particles, 1, function(state) {
    s <- state[1]
    i <- state[2]
    if (i == 0) {
      return(c(s, i))
    }
    epidemic_step(state, t_end = 1, lambda, gamma, n_total)
  }))
  new_particles
}

log_likelihood_fn_epidemic <- function(y, particles, t) {
  # particles is expected to be a matrix with columns (s, i)

  # t=1 deterministic
  if (t == 1) {
    return(rep(0, nrow(particles)))
  }

  dpois(y, lambda = particles[, 2], log = TRUE)
}

# Functions for drawing from the priors
sample_lambda <- function(n) {
  msm::rtnorm(n, mean = 0, sd = 1, lower = 0)
}

sample_gamma <- function(n) {
  msm::rtnorm(n, mean = 0, sd = 1, lower = 0)
}

# --- Stan model definition ---
stan_code <- "
functions {
  vector sir(real t, vector y, array[] real theta,
             array[] real x_r, array[] int x_i) {
      real S = y[1];
      real I = y[2];
      real R = y[3];
      real N = x_i[1];
      real lambda = theta[1];
      real gamma = theta[2];
      real dS_dt = -lambda * I * S / N;
      real dI_dt =  lambda * I * S / N - gamma * I;
      real dR_dt =  gamma * I;
      return [dS_dt, dI_dt, dR_dt]';
  }
}
data {
  int<lower=1> n_days;
  vector[3] y0;
  array[n_days] real t;
  int N;
  array[n_days] int cases;
}
transformed data {
  real t0 = 0;
  array[0] real x_r;
  array[1] int x_i = { N };
}
parameters {
  real<lower=0> lambda;
  real<lower=0> gamma;
}
transformed parameters{
  array[n_days] vector[3] y;
  {
    array[2] real theta;
    theta[1] = lambda;
    theta[2] = gamma;
    y = ode_rk45(sir, y0, t0, t, theta, x_r, x_i);
  }
}
model {
  lambda ~ normal(0, 1);
  gamma ~ normal(0, 1);
  cases ~ poisson(y[,2]);
}
"
# Compile Stan model once
stan_model_fit <- stan_model(model_code = stan_code)

# ---------- Begin the Loop over Repetitions ----------
n_sim <- 100

# Pre-allocate vectors to store mean estimates from both methods
pmcmc_lambda_est <- numeric(n_sim)
pmcmc_gamma_est <- numeric(n_sim)
stan_lambda_est <- numeric(n_sim)
stan_gamma_est <- numeric(n_sim)

for (i in 1:n_sim) {
  cat("Iteration:", i, "\n")

  # Set a different seed for each iteration to get different simulated data
  set.seed(1405 + i)

  ## SIMULATE DATA ##
  true_states <- simulate_epidemic(
    n_total, init_infected, true_lambda, true_gamma, t_max
  )
  latent_i <- true_states[, 2]
  observations <- rpois(
    length(latent_i) - 1,
    lambda = latent_i[-1]
  )

  observations <- c(init_infected, observations)

  ## MODEL FITTING: PMCMC ##
  # Generate 4 samples from the prior for 4 chains
  prior_samples <- data.frame(
    lambda = sample_lambda(4),
    gamma = sample_gamma(4)
  )

  pilot_init_params <- lapply(
    split(
      prior_samples,
      seq_len(nrow(prior_samples))
    ),
    function(row) setNames(as.numeric(row), names(prior_samples))
  )

  # Run PMCMC using pmmh
  result_stochastic_sir <- bayesSSM::pmmh(
    y = observations,
    m = 10000,
    init_fn = init_fn_epidemic,
    transition_fn = transition_fn_epidemic,
    log_likelihood_fn = log_likelihood_fn_epidemic,
    log_priors = log_priors,
    pilot_init_params = pilot_init_params,
    burn_in = 500,
    num_chains = 4,
    param_transform = list(lambda = "log", gamma = "log"),
    verbose = FALSE,
    seed = 1405 + i
  )

  # Extract the posterior means
  pmcmc_lambda_est[i] <- summary(result_stochastic_sir)["lambda", "mean"]
  pmcmc_gamma_est[i] <- summary(result_stochastic_sir)["gamma", "mean"]

  ## MODEL FITTING: Stan ODE ##
  # Prepare data for Stan model
  n_days <- t_max
  t_seq <- seq(1, n_days) # Stan expects a vector of time points
  i0 <- init_infected
  s0 <- n_total - i0
  r0 <- 0
  y0 <- c(S = s0, I = i0, R = r0)

  data_sir <- list(
    n_days = n_days,
    y0 = y0,
    t = t_seq,
    N = n_total,
    cases = observations[-1] # remove time 0 if needed
  )

  # Fit the Stan model
  fit_sir_ode_stan <- sampling(
    object = stan_model_fit,
    data = data_sir,
    seed = 1405 + i,
    iter = 6000,
    chains = 4,
    refresh = 0
  )

  fit_sir_ode_stan

  stan_summary <- summary(fit_sir_ode_stan, pars = c("lambda", "gamma"))$summary
  stan_lambda_est[i] <- stan_summary["lambda", "mean"]
  stan_gamma_est[i] <- stan_summary["gamma", "mean"]
}

# Combine to dataframe and save results
results_ode_vs_stochastic_repeated <- data.frame(
  iteration = 1:n_sim,
  pmcmc_lambda = pmcmc_lambda_est,
  pmcmc_gamma = pmcmc_gamma_est,
  stan_lambda = stan_lambda_est,
  stan_gamma = stan_gamma_est
)

saveRDS(
  results_ode_vs_stochastic_repeated,
  "results_ode_vs_stochastic_repeated.rds"
)

results_ode_vs_stochastic_repeated <- readRDS(
  "results_ode_vs_stochastic_repeated.rds"
)

pmcmc_lambda_est <- results_ode_vs_stochastic_repeated$pmcmc_lambda
pmcmc_gamma_est <- results_ode_vs_stochastic_repeated$pmcmc_gamma
stan_lambda_est <- results_ode_vs_stochastic_repeated$stan_lambda
stan_gamma_est <- results_ode_vs_stochastic_repeated$stan_gamma

# ---------- RMSE Calculation ----------
# Root Mean Squared Error for PMCMC estimates
rmse_pmcmc_lambda <- sqrt(mean((pmcmc_lambda_est - true_lambda)^2))
rmse_pmcmc_gamma <- sqrt(mean((pmcmc_gamma_est - true_gamma)^2))

# Root Mean Squared Error for Stan estimates
rmse_stan_lambda <- sqrt(mean((stan_lambda_est - true_lambda)^2))
rmse_stan_gamma <- sqrt(mean((stan_gamma_est - true_gamma)^2))

sd_pmcmc_lambda <- sd(pmcmc_lambda_est)
sd_pmcmc_gamma <- sd(pmcmc_gamma_est)
sd_stan_lambda <- sd(stan_lambda_est)
sd_stan_gamma <- sd(stan_gamma_est)

# ---------- Display Results ----------
cat("----- RMSE Comparison -----\n")
cat("PMCMC lambda RMSE:", rmse_pmcmc_lambda, "\n")
cat("PMCMC gamma RMSE: ", rmse_pmcmc_gamma, "\n")
cat("Stan lambda RMSE: ", rmse_stan_lambda, "\n")
cat("Stan gamma RMSE:  ", rmse_stan_gamma, "\n")

cat("----- Standard Deviation Comparison -----\n")
cat("PMCMC lambda SD:", sd_pmcmc_lambda, "\n")
cat("PMCMC gamma SD: ", sd_pmcmc_gamma, "\n")
cat("Stan lambda SD: ", sd_stan_lambda, "\n")
cat("Stan gamma SD:  ", sd_stan_gamma, "\n")
