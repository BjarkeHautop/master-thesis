library(ggplot2)
library(ggdist)

library(tibble)
library(extraDistr)

library(bayesSSM) # My R package see https://bjarkehautop.github.io/bayesSSM/

set.seed(1405)

###########################################
# 1. Example Setup: Non-linear Gaussian SSM
###########################################
# SSM definitions:
#    X_0 ~ N(0,1)
#    X_t = phi * X_{t-1} + sin(X_{t-1}) + sigma_x * V_t,   V_t ~ N(0,1), t>=1
#    Y_t = X_t + sigma_y * W_t,              W_t ~ N(0,1), t>=1

init_fn <- function(num_particles) {
  rnorm(num_particles, mean = 0, sd = 1)
}

transition_fn <- function(particles, phi, sigma_x) {
  phi * particles + sin(particles) +
    rnorm(length(particles), mean = 0, sd = sigma_x)
}

log_likelihood_fn <- function(y, particles, sigma_y) {
  dnorm(y, mean = particles, sd = sigma_y, log = TRUE)
}

simulate_ssm <- function(t_val, phi, sigma_x, sigma_y) {
  init_state <- rnorm(1, mean = 0, sd = 1)
  x <- numeric(t_val)
  y <- numeric(t_val)
  x[1] <- phi * init_state + sin(init_state) +
    rnorm(1, mean = 0, sd = sigma_x)
  y[1] <- x[1] + rnorm(1, mean = 0, sd = sigma_y)
  for (t in 2:t_val) {
    x[t] <- phi * x[t - 1] + sin(x[t - 1]) + rnorm(1, mean = 0, sd = sigma_x)
    y[t] <- x[t] + rnorm(1, mean = 0, sd = sigma_y)
  }
  x <- c(init_state, x)
  list(x = x, y = y)
}

###########################################
# 2. Simulate data from the SSM
###########################################

t_val <- 50
phi_true <- 0.7
sigma_x_true <- 1
sigma_y_true <- 1

sim_data <- simulate_ssm(t_val, phi_true, sigma_x_true, sigma_y_true)
x_true <- sim_data$x
y_obs <- sim_data$y

###########################################
# 3. Define functions for the priors
###########################################

# Standard uniform for phi
log_prior_phi <- function(phi) {
  stats::dunif(phi, min = -1, max = 1, log = TRUE)
}

# Half-normal for sigma_x and sigma_y
log_prior_sigma_x <- function(sigma) {
  extraDistr::dhnorm(sigma, sigma = 1, log = TRUE)
}

log_prior_sigma_y <- function(sigma) {
  extraDistr::dhnorm(sigma, sigma = 1, log = TRUE)
}

log_priors <- list(
  phi = log_prior_phi,
  sigma_x = log_prior_sigma_x,
  sigma_y = log_prior_sigma_y
)


###########################################
# 4. Prior predictive check
###########################################

# Define prior sampling functions for phi, sigma_x, and sigma_y
sample_phi <- function(n) {
  runif(n, min = -1, max = 1)
}

sample_sigma_x <- function(n) {
  extraDistr::rhnorm(n)
}

sample_sigma_y <- function(n) {
  extraDistr::rhnorm(n)
}

# Generate 100 prior samples
prior_samples <- data.frame(
  phi = sample_phi(100),
  sigma_x = sample_sigma_x(100),
  sigma_y = sample_sigma_y(100)
)

# Simulate trajectories for the prior predictive check
simulate_trajectory <- function(prior_samples, t_val) {
  prior_trajectories <- list()

  for (i in seq_len(nrow(prior_samples))) {
    phi <- prior_samples$phi[i]
    sigma_x <- prior_samples$sigma_x[i]
    sigma_y <- prior_samples$sigma_y[i]

    sim <- simulate_ssm(t_val, phi, sigma_x, sigma_y)

    prior_trajectories[[i]] <- sim$y
  }

  prior_trajectories
}

# Generate the prior predictive trajectories
prior_trajectories <- simulate_trajectory(prior_samples, t_val)

# Prepare data for plotting
plot_data <- data.frame(
  time = rep(1:t_val, times = length(prior_trajectories) + 1),
  y = c(y_obs, unlist(prior_trajectories)),
  type = rep(
    c(
      "Observed",
      paste0("Sampled ", seq_along(prior_trajectories))
    ),
    each = t_val
  )
)

ggplot() +
  geom_line(
    data = plot_data,
    aes(x = time, y = y, group = type),
    color = "skyblue",
    alpha = 0.7,
    linewidth = 1.2,
    linetype = "dashed"
  ) +
  geom_line(
    data = plot_data[plot_data$type == "Observed", ],
    aes(x = time, y = y),
    color = "darkblue",
    linewidth = 1.5
  ) +
  labs(
    x = "Time",
    y = "Observed",
    title = "Prior Predictive Check"
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
  "outputs/prior_predictive_check_example_4.2.png",
  dpi = 300,
  width = 6.27,
  height = 4,
  units = "in"
)


###########################################
# 5. Run pmmh on dataset
###########################################
# Convert prior_samples to a list and use those as init_params for pilot chain.

pilot_init_params <- split(prior_samples, seq(nrow(prior_samples)))
pilot_init_params <- lapply(
  pilot_init_params,
  function(row) setNames(as.numeric(row), names(prior_samples))
)

# Use first 4 samples for the pilot chain
pilot_init_params <- pilot_init_params[1:4]

result_sir_example_4.2 <- pmmh(
  y = y_obs,
  m = 50000,
  init_fn = init_fn,
  transition_fn = transition_fn,
  log_likelihood_fn = log_likelihood_fn,
  log_priors = log_priors,
  pilot_init_params = pilot_init_params,
  burn_in = 500,
  num_chains = 4,
  algorithm = "SISAR",
  resample_fn = "stratified",
  param_transform = list(
    phi = "identity",
    sigma_x = "log",
    sigma_y = "log"
  ),
  verbose = TRUE,
  return_latent_state_est = TRUE,
  seed = 1405,
  num_cores = 1
)


saveRDS(result_sir_example_4.2, "result_sir_example_4.2.rds")

# Or load the result from a file
result_sir_example_4.2 <- readRDS("result_sir_example_4.2.rds")

###########################################
# 6. Posterior predictive check
###########################################

set.seed(1405)

# Number of posterior predictive simulations to generate:
n_ppc <- 100

chains <- result_sir_example_4.2$theta_chain

num_iterations <- 100
sampled_df <- chains[
  sample(nrow(chains), num_iterations, replace = FALSE),
]

# Simulate posterior predictive trajectories
simulate_trajectory <- function(posterior_samples, t_val) {
  ppc_trajectories <- list()

  for (i in seq_len(nrow(posterior_samples))) {
    phi <- posterior_samples$phi[i]
    sigma_x <- posterior_samples$sigma_x[i]
    sigma_y <- posterior_samples$sigma_y[i]

    sim <- simulate_ssm(t_val, phi, sigma_x, sigma_y)

    ppc_trajectories[[i]] <- sim$y
  }

  ppc_trajectories
}

# Generate the posterior predictive trajectories
ppc_trajectories <- simulate_trajectory(sampled_df, t_val)

# Prepare data for plotting
plot_data_ppc <- data.frame(
  time = rep(1:t_val, times = length(ppc_trajectories) + 1),
  y = c(y_obs, unlist(ppc_trajectories)),
  type = rep(
    c(
      "Observed",
      paste0("Sampled ", seq_along(ppc_trajectories))
    ),
    each = t_val
  )
)

# Plot the posterior predictive check
ggplot() +
  geom_line(
    data = plot_data_ppc,
    aes(x = time, y = y, group = type),
    color = "skyblue",
    alpha = 0.7,
    linewidth = 1.2,
    linetype = "dashed"
  ) +
  geom_line(
    data = plot_data_ppc[plot_data_ppc$type == "Observed", ],
    aes(x = time, y = y),
    color = "darkblue",
    linewidth = 1.5
  ) +
  labs(
    x = "Time",
    y = "Observed",
    title = "Posterior Predictive Check"
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
  "outputs/posterior_predictive_check_example_4.2.png",
  dpi = 300,
  width = 6.27,
  height = 4,
  units = "in"
)

###########################################
# 7. Density plots
###########################################

ggplot(chains, aes(x = phi, fill = chain)) +
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
  "outputs/density_plot_phi_example_4.2.png",
  dpi = 300,
  width = 6.27,
  height = 4,
  units = "in"
)

ggplot(chains, aes(x = sigma_x, fill = chain)) +
  geom_density(alpha = 0.4) +
  theme_minimal() +
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
  "outputs/density_plot_sigma_x_example_4.2.png",
  dpi = 300,
  width = 6.27,
  height = 4,
  units = "in"
)

ggplot(chains, aes(x = sigma_y, fill = chain)) +
  geom_density(alpha = 0.4) +
  theme_minimal() +
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
  "outputs/density_plot_sigma_y_example_4.2.png",
  dpi = 300,
  width = 6.27,
  height = 4,
  units = "in"
)

###########################################
# 8. Inference
###########################################

print(result_sir_example_4.2)

ggplot(chains, aes(x = phi, y = 0)) +
  stat_halfeye(
    point_interval = "mean_qi",
    .width = 0.95,
    slab_fill = "#3B83BD80",
    slab_color = "NA",
    slab_alpha = 0.6,
    interval_size = 1.2,
    point_size = 2.5,
    point_color = "black",
    interval_color = "black"
  ) +
  theme_bw() +
  theme(
    panel.border = element_rect(linewidth = 0.8, colour = "black"),
    panel.grid.major = element_line(linewidth = 0.4, colour = "grey85"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Posterior Density of ϕ",
    x = "Parameter Value",
    y = "Density"
  )

ggsave(
  "outputs/combined_posterior_density_plot_phi_example_4.2.png",
  dpi = 300,
  width = 6.27,
  height = 4,
  units = "in"
)

ggplot(chains, aes(x = sigma_x, y = 0)) +
  stat_halfeye(
    point_interval = "mean_qi",
    .width = 0.95,
    slab_fill = "#3B83BD80",
    slab_color = "NA",
    slab_alpha = 0.6,
    interval_size = 1.2,
    point_size = 2.5,
    point_color = "black",
    interval_color = "black"
  ) +
  theme_bw() +
  theme(
    panel.border = element_rect(linewidth = 0.8, colour = "black"),
    panel.grid.major = element_line(linewidth = 0.4, colour = "grey85"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Posterior Density of σₓ",
    x = "Parameter Value",
    y = "Density"
  )

ggsave(
  "outputs/combined_posterior_density_plot_sigma_x_example_4.2.png",
  dpi = 300,
  width = 6.27,
  height = 4,
  units = "in"
)


ggplot(chains, aes(x = sigma_y, y = 0)) +
  stat_halfeye(
    point_interval = "mean_qi",
    .width = 0.95,
    slab_fill = "#3B83BD80",
    slab_color = "NA",
    slab_alpha = 0.6,
    interval_size = 1.2,
    point_size = 2.5,
    point_color = "black",
    interval_color = "black"
  ) +
  theme_bw() +
  theme(
    panel.border = element_rect(linewidth = 0.8, colour = "black"),
    panel.grid.major = element_line(linewidth = 0.4, colour = "grey85"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Posterior Density of σᵧ",
    x = "Parameter Value",
    y = "Density"
  )

ggsave(
  "outputs/combined_posterior_density_plot_sigma_y_example_4.2.png",
  dpi = 300,
  width = 6.27,
  height = 4,
  units = "in"
)

###########################################
# 9. Latent state
###########################################

# Get number of outer lists and particles per list
n_lists <- length(result_sir_example_4.2$latent_state_chain)
n_particles <- length(result_sir_example_4.2$latent_state_chain[[1]])
chain_length <- length(result_sir_example_4.2$latent_state_chain[[1]][[1]])

# Initialize matrix to store all chains (rows: 40000, columns: chain_length)
all_chains_matrix <- matrix(NA, nrow = n_lists * n_particles, ncol = chain_length)

# Fill the matrix with chains
row_idx <- 1
for (l in 1:n_lists) {
  for (p in 1:n_particles) {
    all_chains_matrix[row_idx, ] <- result_sir_example_4.2$latent_state_chain[[l]][[p]]
    row_idx <- row_idx + 1
  }
}

means <- colMeans(all_chains_matrix)

rmse <- function(est, true) sqrt(mean((est - true)^2))

rmse_value <- rmse(means, x_true)

cat("RMSE of the latent state estimate:", rmse_value, "\n")
