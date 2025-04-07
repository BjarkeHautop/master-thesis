library(ggplot2)
library(tibble)
library(extraDistr)

library(bayesSSM) # My R package see https://bjarkehautop.github.io/bayesSSM/

set.seed(1405)

###########################################
# 1. Example Setup: Non-linear Gaussian SSM
###########################################
# SSM definitions:
#    X_1 ~ N(0,1)
#    X_t = phi * X_{t-1} + sin(X_{t-1}) + sigma_x * V_t,   V_t ~ N(0,1)
#    Y_t = X_t + sigma_y * W_t,              W_t ~ N(0,1)

init_fn <- function(particles) {
  rnorm(particles, mean = 0, sd = 1)
}

transition_fn <- function(particles, phi, sigma_x) {
  phi * particles + sin(particles) +
    rnorm(length(particles), mean = 0, sd = sigma_x)
}

log_likelihood_fn <- function(y, particles, sigma_y) {
  dnorm(y, mean = particles, sd = sigma_y, log = TRUE)
}

simulate_ssm <- function(t_val, phi, sigma_x, sigma_y) {
  x <- numeric(t_val)
  y <- numeric(t_val)
  x[1] <- rnorm(1, mean = 0, sd = sigma_x)
  y[1] <- rnorm(1, mean = x[1], sd = sigma_y)
  for (t in 2:t_val) {
    x[t] <- phi * x[t - 1] + sin(x[t - 1]) + rnorm(1, mean = 0, sd = sigma_x)
    y[t] <- x[t] + rnorm(1, mean = 0, sd = sigma_y)
  }
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

# Standard normal for phi
log_prior_phi <- function(phi) {
  stats::dnorm(phi, mean = 0, sd = 1, log = TRUE)
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
# 4. Prior predictice check
###########################################

# Define prior sampling functions for phi, sigma_x, and sigma_y
sample_phi <- function(n) {
  rnorm(n, mean = 0, sd = 1)
}

sample_sigma_x <- function(n) {
  extraDistr::rhnorm(n)
}

sample_sigma_y <- function(n) {
  extraDistr::rhnorm(n)
}

# Generate 4 prior samples
prior_samples <- data.frame(
  phi = sample_phi(4),
  sigma_x = sample_sigma_x(4),
  sigma_y = sample_sigma_y(4)
)

# Simulate trajectories for the prior predictive check
simulate_trajectory <- function(prior_samples, t_val) {
  prior_trajectories <- list()

  for (i in 1:nrow(prior_samples)) {
    phi <- prior_samples$phi[i]
    sigma_x <- prior_samples$sigma_x[i]
    sigma_y <- prior_samples$sigma_y[i]

    sim <- simulate_ssm(t_val, phi, sigma_x, sigma_y)

    prior_trajectories[[i]] <- sim$y
  }

  return(prior_trajectories)
}

# Generate the prior predictive trajectories
prior_trajectories <- simulate_trajectory(prior_samples, t_val)

# Prepare data for plotting
plot_data <- data.frame(
  time = rep(1:t_val, times = length(prior_trajectories) + 1),
  y = c(y_obs, unlist(prior_trajectories)),
  type = rep(c(
    "Observed",
    paste0("Sampled ", 1:length(prior_trajectories))),
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
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid = element_blank()
  )

ggsave(
  "prior_predictive_check_example_4.3.png",
  dpi = 300,
  width = 6.27,
  height = 4,
  units = "in"
)


###########################################
# 5. Run pmmh on dataset
###########################################

result_sir_example_4.3 <- pmmh(
  y = y_obs,
  m = 10000,
  init_fn = init_fn,
  transition_fn = transition_fn,
  log_likelihood_fn = log_likelihood_fn,
  log_priors = log_priors,
  init_params = c(phi = 0.5, sigma_x = 0.5, sigma_y = 0.5),
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
  num_cores = 4
)

saveRDS(result_sir_example_4.3, "result_sir_example_4.3.rds")

# Or load the result from a file
# result_sir_example_4.3 <- readRDS("result_sir_example_4.3.rds")

###########################################
# 6. Posterior predictive check
###########################################

set.seed(1405)

# Number of posterior predictive simulations to generate:
n_ppc <- 4

# Use chain 1
posterior_samples <- as.data.frame(result_sir_example_4.3$theta_chain[[1]])

# Sample posterior parameter values
params <- posterior_samples[sample(nrow(posterior_samples), n_ppc), ]

# Simulate posterior predictive trajectories
simulate_trajectory <- function(posterior_samples, t_val) {
  ppc_trajectories <- list()

  for (i in 1:nrow(posterior_samples)) {
    phi <- posterior_samples$phi[i]
    sigma_x <- posterior_samples$sigma_x[i]
    sigma_y <- posterior_samples$sigma_y[i]

    sim <- simulate_ssm(t_val, phi, sigma_x, sigma_y)

    ppc_trajectories[[i]] <- sim$y
  }

  return(ppc_trajectories)
}

# Generate the posterior predictive trajectories
ppc_trajectories <- simulate_trajectory(params, t_val)

# Prepare data for plotting
plot_data_ppc <- data.frame(
  time = rep(1:t_val, times = length(ppc_trajectories) + 1),
  y = c(y_obs, unlist(ppc_trajectories)),
  type = rep(c(
    "Observed",
    paste0("Sampled ", 1:length(ppc_trajectories))),
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
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid = element_blank()
  )

ggsave(
  "posterior_predictive_check_example_4.3.png",
  dpi = 300,
  width = 6.27,
  height = 4,
  units = "in"
)

###########################################
# 7. Density plots
###########################################

# Extract the posterior samples from the chain
chain_1 <- as.data.frame(result_sir_example_4.3$theta_chain[[1]])
chain_2 <- as.data.frame(result_sir_example_4.3$theta_chain[[2]])
chain_3 <- as.data.frame(result_sir_example_4.3$theta_chain[[3]])
chain_4 <- as.data.frame(result_sir_example_4.3$theta_chain[[4]])

phi_chains <- cbind(
  chain_1$phi, chain_2$phi,
  chain_3$phi, chain_4$phi
)

sigma_x_chains <- cbind(
  chain_1$sigma_x, chain_2$sigma_x,
  chain_3$sigma_x, chain_4$sigma_x
)

sigma_y_chains <- cbind(
  chain_1$sigma_y, chain_2$sigma_y,
  chain_3$sigma_y, chain_4$sigma_y
)

phi_df <- data.frame(
  value = c(
    phi_chains[, 1],
    phi_chains[, 2],
    phi_chains[, 3],
    phi_chains[, 4]
  ),
  chain = factor(rep(1:4, each = nrow(phi_chains)))
)

ggplot(phi_df, aes(x = value, fill = chain)) +
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

ggsave("density_plot_phi_example_4.3.png",
       dpi = 300,
       width = 6.27,
       height = 4,
       units = "in"
)

sigma_x_df <- data.frame(
  value = c(
    sigma_x_chains[, 1],
    sigma_x_chains[, 2],
    sigma_x_chains[, 3],
    sigma_x_chains[, 4]
  ),
  chain = factor(rep(1:4, each = nrow(sigma_x_chains)))
)

ggplot(sigma_x_df, aes(x = value, fill = chain)) +
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

ggsave("density_plot_sigma_x_example_4.3.png",
       dpi = 300,
       width = 6.27,
       height = 4,
       units = "in"
)

sigma_y_df <- data.frame(
  value = c(
    sigma_y_chains[, 1],
    sigma_y_chains[, 2],
    sigma_y_chains[, 3],
    sigma_y_chains[, 4]
  ),
  chain = factor(rep(1:4, each = nrow(sigma_y_chains)))
)

ggplot(sigma_y_df, aes(x = value, fill = chain)) +
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

ggsave("density_plot_sigma_y_example_4.3.png",
       dpi = 300,
       width = 6.27,
       height = 4,
       units = "in"
)

###########################################
# 8. Inference
###########################################

print(result_sir_example_4.3)
