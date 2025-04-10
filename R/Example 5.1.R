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
true_phi <- 3.5 # True dispersion parameter

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
  states <- matrix(0, nrow = t_max + 1, ncol = 2)
  # initial state at t = 0
  states[1, ] <- c(n_total - init_infected, init_infected)
  state <- states[1, ]
  for (t in 1:t_max) {
    state <- epidemic_step(state, lambda, gamma, n_total)
    states[t + 1, ] <- state
  }
  states
}

# Simulate an epidemic dataset
true_states <- simulate_epidemic(
  n_total, init_infected, true_lambda, true_gamma, t_max
)
latent_i <- true_states[, 2]

observations <- rnbinom(length(latent_i), size = true_phi, mu = latent_i)

data_sim <- data.frame(
  time = 0:t_max,
  s = true_states[, 1],
  i = true_states[, 2],
  y = observations
)
print(data_sim)

prepare_data_for_plot <- function(states, observations, t_max) {
  data <- data.frame(
    time = 0:t_max,
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
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.position = "top",
    ) +
    scale_x_continuous(breaks = seq(0, t_max, by = 1))
}

data_long <- prepare_data_for_plot(true_states, observations, t_max)

data_long$state <- factor(data_long$state, levels = c("s", "i", "y"))


plot_epidemic_data(data_long, t_max)

ggsave(
  "epidemic_data_plot_negbin.png",
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

# 1 / phi ~ Normal+(0, 1)
log_prior_phi <- function(phi) {
  if (phi <= 0) return(-Inf)
  xi <- 1 / phi
  log_prior_xi <- msm::dtnorm(xi, mean = 0, sd = 1, log = TRUE)

  # log_jacobian: |d(1/phi)/dphi| = 1/phi^2
  log_jacobian <- -2 * log(phi)

  log_prior_xi + log_jacobian
}



log_priors <- list(
  lambda = log_prior_lambda,
  gamma = log_prior_gamma,
  phi = log_prior_phi
)

# --- Functions for the epidemic model ---

init_fn_epidemic <- function(particles) {
  # Return a matrix with particles rows; each row is the initial state (s, i)
  matrix(rep(init_state, each = particles), nrow = particles, byrow = FALSE)
}

transition_fn_gillespie <- function(particles, lambda, gamma) {
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

log_likelihood_fn_epidemic <- function(y, particles, phi) {
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
  inv_phi <- msm::rtnorm(n, mean = 0, sd = 1, lower = 0)
  1 / inv_phi
}

# Generate 20 samples and plot the prior predictive distribution along the
# observed data
prior_samples <- data.frame(
  lambda = sample_lambda(20),
  gamma = sample_gamma(20),
  phi = sample_phi(20)
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
      length(states[, 2]), size = phi, mu = states[, 2]
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
  time = rep(0:t_max, times = length(prior_trajectories) + 1),
  infected = c(observations, unlist(prior_trajectories)),
  type = rep(
    c(
      "Observed",
      paste0("Sampled ", seq_along(prior_trajectories))
    ),
    each = t_max + 1
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
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid = element_blank()
  )

# Save the plot
ggsave(
  "prior_predictive_check_sir_example_5.1_negbin.png",
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

# Use 4 rows of prior_samples as initial values for the pilot chain
values <- c(1:3, 14)
pilot_init_params <- pilot_init_params[values]
pilot_init_params

result_sir_example_5.1_negbin <- bayesSSM::pmmh(
  y = observations,
  m = 10000,
  init_fn = init_fn_epidemic,
  transition_fn = transition_fn_gillespie,
  log_likelihood_fn = log_likelihood_fn_epidemic,
  log_priors = log_priors,
  pilot_init_params = pilot_init_params,
  burn_in = 100,
  num_chains = 4,
  param_transform = list(lambda = "log", gamma = "log", phi = "log"),
  verbose = TRUE,
  return_latent_state_est = TRUE,
  seed = 1405,
  num_cores = 4
)

saveRDS(result_sir_example_5.1_negbin, "result_sir_example_5.1_negbin.rds")

# Or load the result from a file
result_sir_example_5.1_negbin <- readRDS("result_sir_example_5.1_negbin.rds")

chains <- result_sir_example_5.1_negbin$theta_chain

ggplot(chains, aes(x = lambda, fill = factor(chain))) +
  geom_density(alpha = 0.4) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 12),
    panel.grid.major = element_line(linewidth = 0.5, color = "gray80"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Parameter Value",
    y = "Density"
  )

ggsave(
  "density_plot_example_5.1_lambda_negbin.png",
  dpi = 300,
  width = 6.27,
  height = 4,
  units = "in"
)

ggplot(chains, aes(x = gamma, fill = factor(chain))) +
  geom_density(alpha = 0.4) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 12),
    panel.grid.major = element_line(linewidth = 0.5, color = "gray80"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Parameter Value",
    y = "Density"
  )

ggsave(
  "density_plot_example_5.1_gamma_negbin.png",
  dpi = 300,
  width = 6.27,
  height = 4,
  units = "in"
)

ggplot(chains, aes(x = phi, fill = factor(chain))) +
  geom_density(alpha = 0.4) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 12),
    panel.grid.major = element_line(linewidth = 0.5, color = "gray80"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Parameter Value",
    y = "Density"
  )

ggsave(
  "density_plot_example_5.1_phi_negbin.png",
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
  sample(nrow(chains), num_iterations, replace = TRUE),
]

posterior_trajectories <- simulate_trajectory(
  sampled_df, n_total, init_infected, t_max
)

# Prepare the data for plotting
plot_data_posterior <- data.frame(
  time = rep(0:t_max, times = length(posterior_trajectories) + 1),
  infected = c(observations, unlist(posterior_trajectories)),
  type = rep(
    c(
      "Observed",
      paste0("Sampled ", seq_along(posterior_trajectories))),
    each = t_max + 1
  )
)

# Posterior Predictive Check Plot
ggplot() +
  geom_line(
    data = subset(plot_data_posterior, type != "Observed"),
    aes(x = time, y = infected, group = type),
    color = "skyblue",
    alpha = 0.7,
    linewidth = 1.2,
    linetype = "dashed"
  ) +
  geom_line(
    data = subset(plot_data_posterior, type == "Observed"),
    aes(x = time, y = infected),
    color = "darkblue",
    linewidth = 1.5
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
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid = element_blank()
  )

ggsave(
  "posterior_predictive_check_example_5.1_negbin.png",
  dpi = 300,
  width = 6.27,
  height = 4,
  units = "in"
)

###########################################
# 6. Future inference
###########################################

# --- Estimate future trajectory of the epidemic ---

# We sample from the posterior distribution of the parameters then simulate
# the epidemic trajectory using the estimated parameters.

results_list <- list()

# Loop over each chain
for (chain_idx in seq_along(result_sir_example_5.1_negbin$latent_state_chain)) {
  chain <- result_sir_example_5.1_negbin$latent_state_chain[[chain_idx]]

  # Extract the [11,1] and [11,2] values for each matrix in the chain
  extracted_values <- do.call(
    rbind, lapply(chain, function(mat) round(mat[11, ]))
  )

  # Convert to a dataframe and add chain identifier
  df <- as.data.frame(extracted_values)
  df$chain <- chain_idx

  # Store in results list
  results_list[[chain_idx]] <- df
}

# Combine all chains into a single dataframe
latent_est_df <- do.call(rbind, results_list)

# Rename columns for clarity
colnames(latent_est_df) <- c("S", "I", "chain")


# Combine the latent state estimates with the sampled parameters
# Add latent_est_df to `chains`
chains$S <- latent_est_df$S
chains$I <- latent_est_df$I

num_iterations <- nrow(chains)

future_inference_df <- data.frame(
  iteration = 1:num_iterations,
  num_days = rep(0, num_iterations),
  not_infected = rep(0, num_iterations)
)

for (i in 1:num_iterations) {
  state <- as.numeric(chains[i, c("S", "I")])
  lambda <- chains$lambda[i]
  gamma <- chains$gamma[i]

  num_days <- 10
  while (state[2] > 0) { # Continue until the number of infected reaches 0
    num_days <- num_days + 1
    state <- epidemic_step(state, lambda = lambda, gamma = gamma, n_total)
  }
  future_inference_df$num_days[i] <- num_days
  future_inference_df$not_infected[i] <- state[1]
}

mean(future_inference_df$num_days)
mean(n_total - future_inference_df$not_infected)

# Credible interval
quantile(future_inference_df$num_days, c(0.025, 0.975))
quantile(n_total - future_inference_df$not_infected, c(0.025, 0.975))
