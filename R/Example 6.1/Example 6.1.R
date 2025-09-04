library(bayesSSM) # My R package see https://bjarkehautop.github.io/bayesSSM/
library(ggplot2)
library(tidyr)
library(extraDistr)
set.seed(1405)

###########################################
# 1. Simulate data
###########################################

# --- Simulation settings and true parameters ---
n_total <- 500
init_infected <- 10
init_state <- c(n_total - init_infected, init_infected)
t_max <- 10 # Total number of days to simulate
true_lambda <- 0.5 # True infection parameter
true_gamma <- 0.2 # True removal parameter
true_phi <- 6 # True dispersion parameter

# --- Functions for simulating the epidemic ---

epidemic_step <- function(state, lambda, gamma, n_total) {
  t <- 0
  t_end <- 1
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
      # Infection event
      s <- s - 1
      i <- i + 1
    } else {
      # Removal event
      i <- i - 1
    }
  }
  c(s, i)
}

simulate_epidemic <- function(
    n_total, init_infected, lambda, gamma, t_max) {
  states <- matrix(0, nrow = t_max, ncol = 2)
  # initial state at t = 0
  state <- c(n_total - init_infected, init_infected)
  for (t in 1:t_max) {
    state <- epidemic_step(state, lambda, gamma, n_total)
    states[t, ] <- state
  }
  states
}

# Simulate an epidemic dataset
true_states <- simulate_epidemic(
  n_total, init_infected, true_lambda, true_gamma, t_max
)
latent_i <- true_states[, 2]

# 1st observation is deterministic.

observations <- rnbinom(
  length(latent_i),
  size = true_phi,
  mu = latent_i
)
data_sim <- data.frame(
  time = 1:t_max,
  s = true_states[, 1],
  i = true_states[, 2],
  y = observations
)
print(data_sim)

prepare_data_for_plot <- function(states, observations, t_max) {
  data <- data.frame(
    time = 1:t_max,
    s = states[, 1],
    i = states[, 2],
    y = observations
  )

  data_long <- pivot_longer(
    data,
    cols = c("s", "i", "y"),
    names_to = "state",
    values_to = "count"
  )

  data_long
}

plot_epidemic_data <- function(data_long, t_max) {
  ggplot(data_long, aes(x = time, y = count, color = state)) +
    geom_line(data = subset(data_long, state != "y"), linewidth = 1.2) +
    geom_point(data = subset(data_long, state == "y"), size = 3) +
    scale_color_manual(
      values = c("s" = "#0072B2", "i" = "#D55E00", "y" = "#009E73"),
      labels = c("Susceptible", "Infected", "Observed")
    ) +
    labs(
      x = "Time (Days)",
      y = "Count",
      title = "Susceptible, Infected, and Observed Counts",
      color = "State"
    ) +
    theme_bw() +
    theme(
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(),
      legend.text = element_text(),
      legend.position = "top",
    ) +
    scale_x_continuous(breaks = seq(1, t_max, by = 1))
}

data_long <- prepare_data_for_plot(true_states, observations, t_max)

data_long$state <- factor(data_long$state, levels = c("s", "i", "y"))


plot_epidemic_data(data_long, t_max)

ggsave(
  "outputs/epidemic_data_plot_negbin.png",
  dpi = 300,
  width = 6.27,
  height = 4,
  units = "in"
)

###########################################
# 2. Priors and Model Specification
###########################################

# --- Priors ---

# Define the log-prior for the parameters
log_prior_lambda <- function(lambda) {
  extraDistr::dhnorm(lambda, sigma = 1, log = TRUE)
}

log_prior_gamma <- function(gamma) {
  extraDistr::dhnorm(gamma, sigma = 1, log = TRUE)
}

# 1 / sqrt(phi) ~ Normal+(0, 1)
log_prior_phi <- function(phi) {
  if (phi <= 0) {
    return(-Inf)
  }

  xi <- 1 / sqrt(phi)
  # half‑Normal(0,1):
  log_prior_xi <- msm::dtnorm(
    xi,
    mean  = 0,
    sd    = 1,
    lower = 0,
    upper = Inf,
    log   = TRUE
  )

  # Jacobian of xi = phi^(–1/2):
  #   |d xi / d phi| = 1 / (2 * phi^(3/2))
  log_jacobian <- -log(2) - 1.5 * log(phi)

  log_prior_xi + log_jacobian
}

log_priors <- list(
  lambda = log_prior_lambda,
  gamma = log_prior_gamma,
  phi = log_prior_phi
)

# --- Functions for the epidemic model ---

init_fn_epidemic <- function(num_particles) {
  # Return a matrix with num_particles rows; each row is the initial state (s, i)
  n_total <- 500
  init_infected <- 10
  init_state <- c(S = n_total - init_infected, I = init_infected)
  matrix(rep(init_state, each = num_particles), nrow = num_particles, byrow = FALSE)
}

transition_fn_epidemic <- function(particles, lambda, gamma, t) {
  n_total <- 500

  epidemic_step <- function(state, lambda, gamma, n_total) {
    t <- 0
    t_end <- 1
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
        # Infection event
        s <- s - 1
        i <- i + 1
      } else {
        # Removal event
        i <- i - 1
      }
    }
    c(s, i)
  }

  new_particles <- t(apply(particles, 1, function(state) {
    s <- state[1]
    i <- state[2]
    if (i == 0) {
      return(c(s, i))
    }
    epidemic_step(state, lambda, gamma, n_total)
  }))
  new_particles
}

log_likelihood_fn_epidemic <- function(y, particles, phi, t) {
  # particles is expected to be a matrix with columns (s, i)

  dnbinom(
    y,
    size = phi,
    mu = particles[, 2],
    log = TRUE
  )
}

###########################################
# 3. Prior predictive check
###########################################
sample_gamma <- function(n) {
  rhnorm(n, sigma = 1)
}

sample_lambda <- function(n) {
  rhnorm(n, sigma = 1)
}

sample_phi <- function(n) {
  inv_sqrt_phi <- msm::rtnorm(
    n,
    mean = 0,
    sd = 1,
    lower = 0,
    upper = Inf
  )

  1 / (inv_sqrt_phi^2)
}

# Generate 100 samples and plot the prior predictive distribution along the
# observed data
prior_samples <- data.frame(
  lambda = sample_lambda(100),
  gamma = sample_gamma(100),
  phi = sample_phi(100)
)

# Simulate trajectories for the prior predictive check
simulate_trajectory <- function(
    prior_samples, n_total, init_infected, t_max) {
  prior_trajectories <- list()

  for (i in seq_len(nrow(prior_samples))) {
    lambda <- prior_samples$lambda[i]
    gamma <- prior_samples$gamma[i]
    phi <- prior_samples$phi[i]

    states <- simulate_epidemic(
      n_total, init_infected, lambda, gamma, t_max
    )

    noisy_infected <- rnbinom(
      length(latent_i),
      size = true_phi,
      mu = latent_i
    )

    prior_trajectories[[i]] <- noisy_infected
  }

  prior_trajectories
}

# Generate the prior predictive trajectories
prior_trajectories <- simulate_trajectory(
  prior_samples, n_total, init_infected, t_max
)

plot_data <- data.frame(
  time = rep(1:t_max, times = length(prior_trajectories) + 1),
  infected = c(observations, unlist(prior_trajectories)),
  type = rep(
    c(
      "Observed",
      paste0("Sampled ", seq_along(prior_trajectories))
    ),
    each = t_max
  )
)

ggplot() +
  geom_line(
    data = plot_data,
    aes(x = time, y = infected, group = type),
    color = "skyblue",
    alpha = 0.7,
    linewidth = 1.2,
    linetype = "dashed"
  ) +
  geom_line(
    data = plot_data[plot_data$type == "Observed", ],
    aes(x = time, y = infected),
    color = "darkblue",
    linewidth = 1.5
  ) +
  labs(
    x = "Time (Days)",
    y = "Observed Infected",
    title = "Prior Predictive Check"
  ) +
  scale_x_continuous(
    breaks = seq(0, 10, by = 1)
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid = element_blank()
  )

# Save the plot
ggsave(
  "outputs/prior_predictive_check_sir_example_6.1_negbin.png",
  dpi = 300,
  width = 6.27,
  height = 4,
  units = "in"
)



###########################################
# 4. PMMH
###########################################

# Convert prior_samples to a list and use those as init_params for pilot chain.

pilot_init_params <- split(prior_samples, seq(nrow(prior_samples)))
pilot_init_params <- lapply(
  pilot_init_params,
  function(row) setNames(as.numeric(row), names(prior_samples))
)

# Use first 4 as init_params for pilot chain
pilot_init_params <- pilot_init_params[1:4]
pilot_init_params

result_sir_example_6.1_negbin <- pmmh(
  y = observations,
  m = 40000,
  init_fn = init_fn_epidemic,
  transition_fn = transition_fn_epidemic,
  log_likelihood_fn = log_likelihood_fn_epidemic,
  log_priors = log_priors,
  pilot_init_params = pilot_init_params,
  burn_in = 500,
  num_chains = 4,
  param_transform = list(lambda = "log", gamma = "log", phi = "log"),
  verbose = TRUE,
  return_latent_state_est = TRUE,
  seed = 1405,
  num_cores = 4
)

saveRDS(
  result_sir_example_6.1_negbin,
  file = "outputs/result_sir_example_6.1_negbin.rds"
)

# Or load the result from a file
result_sir_example_6.1_negbin <- readRDS(
  "outputs/result_sir_example_6.1_negbin.rds"
)

chains <- result_sir_example_6.1_negbin$theta_chain

ggplot(chains, aes(x = lambda, fill = factor(chain))) +
  geom_density(alpha = 0.4) +
  theme_minimal() +
  theme(
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_line(linewidth = 0.5, color = "gray80"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Parameter Value",
    y = "Density"
  )

ggsave(
  "outputs/density_plot_example_6.1_lambda_negbin.png",
  dpi = 300,
  width = 6.27,
  height = 4,
  units = "in"
)

ggplot(chains, aes(x = gamma, fill = factor(chain))) +
  geom_density(alpha = 0.4) +
  theme_minimal() +
  theme(
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_line(linewidth = 0.5, color = "gray80"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Parameter Value",
    y = "Density"
  )

ggsave(
  "outputs/density_plot_example_6.1_gamma_negbin.png",
  dpi = 300,
  width = 6.27,
  height = 4,
  units = "in"
)

ggplot(chains, aes(x = phi, fill = factor(chain))) +
  geom_density(alpha = 0.4) +
  theme_minimal() +
  theme(
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_line(linewidth = 0.5, color = "gray80"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Parameter Value",
    y = "Density"
  )

ggsave(
  "outputs/density_plot_example_6.1_phi_negbin.png",
  dpi = 300,
  width = 6.27,
  height = 4,
  units = "in"
)

#########################################
# 5. Posterior Predictive Check
#########################################

set.seed(1405)

num_iterations <- 100
sampled_df <- chains[
  sample(nrow(chains), num_iterations, replace = FALSE),
]

posterior_trajectories <- simulate_trajectory(
  sampled_df, n_total, init_infected, t_max
)

# Prepare the data for plotting
plot_data_posterior <- data.frame(
  time = rep(1:t_max, times = length(posterior_trajectories) + 1),
  infected = c(observations, unlist(posterior_trajectories)),
  type = rep(
    c(
      "Observed",
      paste0("Sampled ", seq_along(posterior_trajectories))
    ),
    each = t_max
  )
)

# Posterior Predictive Check Plot
ggplot() +
  geom_line(
    data = subset(plot_data_posterior, type != "Observed"),
    aes(x = time, y = infected, group = type),
    color = "skyblue",
    alpha = 0.7,
    linewidth = 1,
    linetype = "dashed"
  ) +
  geom_line(
    data = subset(plot_data_posterior, type == "Observed"),
    aes(x = time, y = infected),
    color = "darkblue",
    linewidth = 1
  ) +
  labs(
    x = "Time (Days)",
    y = "Observed Infected",
    title = "Posterior Predictive Check"
  ) +
  scale_x_continuous(
    breaks = seq(0, 14, by = 1)
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid = element_blank()
  )

ggsave(
  "outputs/posterior_predictive_check_example_6.1_negbin.png",
  dpi = 300,
  width = 6.27,
  height = 4,
  units = "in"
)
