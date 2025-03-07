library(ggplot2)
library(tibble)
library(extraDistr) # for halfnormal distribution

library(bayesSSM) # My own package.
# See https://github.com/BjarkeHautop/bayesSSM for installation instructions

set.seed(1405)

###########################################
# 1. Example Setup: Non-linear Gaussian SSM
###########################################
# SSM definitions:
#    X_1 ~ N(0,1)
#    X_t = phi * X_{t-1} + sin(X_{t-1}) + sigma_x * V_t,   V_t ~ N(0,1)
#    Y_t = X_t + sigma_y * W_t,              W_t ~ N(0,1)

init_fn_ssm <- function(particles) {
  rnorm(particles, mean = 0, sd = 1)
}

transition_fn_ssm <- function(particles, phi, sigma_x) {
  phi * particles + sin(particles) +
    rnorm(length(particles), mean = 0, sd = sigma_x)
}

log_likelihood_fn_ssm <- function(y, particles, sigma_y) {
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

# Define a function to sample from the prior
sample_prior <- function() {
  phi <- rnorm(1, mean = 0, sd = 1)
  sigma_x <- extraDistr::rhnorm(1)
  sigma_y <- extraDistr::rhnorm(1)
  list(phi = phi, sigma_x = sigma_x, sigma_y = sigma_y)
}

# Number of prior predictive simulations
n_sim <- 4

sim_list <- vector("list", n_sim)

for (i in 1:n_sim) {
  params <- sample_prior()
  sim <- simulate_ssm(t_val, params$phi, params$sigma_x, params$sigma_y)

  sim_list[[i]] <- tibble(
    time = 1:t_val,
    y = sim$y,
    sim = as.factor(i)
  )
}

sim_df <- dplyr::bind_rows(sim_list)

obs_df <- tibble(time = 1:t_val, y = y_obs)
ggplot() +
  geom_line(
    data = sim_df,
    aes(x = time, y = y, group = sim),
    color = "grey",
    alpha = 0.8,
    linewidth = 0.8
  ) +
  geom_line(
    data = obs_df,
    aes(x = time, y = y),
    color = "red",
    linewidth = 1.2
  ) +
  labs(x = "Time", y = "Observations y") +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 12),
    panel.grid.major = element_line(size = 0.5, color = "gray80"),
    panel.grid.minor = element_blank()
  )

ggsave("prior_predictive_check.png",
  dpi = 300,
  width = 6.27,
  height = 4,
  units = "in"
)


###########################################
# 5. Run pmmh on dataset
###########################################

result <- pmmh(
  y = y_obs,
  m = 15000,
  init_fn_ssm = init_fn_ssm,
  transition_fn_ssm = transition_fn_ssm,
  log_likelihood_fn_ssm = log_likelihood_fn_ssm,
  log_priors = log_priors,
  init_params = c(phi = 0.8, sigma_x = 0.6, sigma_y = 0.6),
  burn_in = 2000,
  num_chains = 1,
  algorithm = "SISAR",
  resample_fn = "stratified",
  param_transform = c("identity", "log", "log"),
  verbose = TRUE,
  seed = 1405
)

###########################################
# 6. Posterior predictive check
###########################################

# Number of posterior predictive simulations to generate:
n_ppc <- 4

# Use chain 1
posterior_samples <- as.data.frame(result$theta_chain[[1]])

ppc_list <- vector("list", n_ppc)

# Randomly sample indices from the posterior chain
params <- posterior_samples[sample(nrow(posterior_samples), 4), ]
# Loop over the selected posterior samples and simulate new datasets:
for (i in 1:n_ppc) {
  # Extract a set of parameters from the posterior sample
  param <- params[i, ]
  # Use simulate_ssm to generate replicated data using these parameters
  sim <- simulate_ssm(t_val, param$phi, param$sigma_x, param$sigma_y)

  ppc_list[[i]] <- tibble(
    time = 1:t_val,
    y = sim$y,
    sim = as.factor(i)
  )
}

# Combine all simulated datasets into one data frame
ppc_df <- dplyr::bind_rows(ppc_list)

# Plot the posterior predictive simulations along with the observed data:
obs_df <- tibble(time = 1:t_val, y = y_obs)
ggplot() +
  geom_line(
    data = ppc_df, aes(x = time, y = y, group = sim),
    color = "grey", alpha = 0.8, size = 0.8
  ) +
  geom_line(
    data = obs_df, aes(x = time, y = y),
    color = "red", size = 1.2
  ) +
  labs(x = "Time", y = "Observed y") +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 12),
    panel.grid.major = element_line(size = 0.5, color = "gray80"),
    panel.grid.minor = element_blank()
  )

ggsave("posterior_predictive_check.png",
  dpi = 300,
  width = 6.27,
  height = 4,
  units = "in"
)

###########################################
# 7. Trace plot
###########################################

# Extract the posterior samples from the chain
chain_1 <- as.data.frame(result$theta_chain[[1]])
chain_2 <- as.data.frame(result$theta_chain[2])
chain_3 <- as.data.frame(result$theta_chain[3])
chain_4 <- as.data.frame(result$theta_chain[4])

# Trace plot of chain_1
ggplot(chain_1, aes(x = seq_len(nrow(chain_1)), y = phi)) +
  geom_line() +
  labs(x = "Iteration", y = "phi") +
  theme_minimal()

ggsave("trace_plot_phi.png",
  dpi = 300,
  width = 6.27,
  height = 4,
  units = "in"
)

ggplot(chain_1, aes(x = seq_len(nrow(chain_1)), y = sigma_x)) +
  geom_line() +
  labs(x = "Iteration", y = "sigma_x") +
  theme_minimal()

ggsave("trace_plot_sigma_x.png",
  dpi = 300,
  width = 6.27,
  height = 4,
  units = "in"
)

ggplot(chain_1, aes(x = seq_len(nrow(chain_1)), y = sigma_y)) +
  geom_line() +
  labs(x = "Iteration", y = "sigma_y") +
  theme_minimal()

ggsave("trace_plot_sigma_y.png",
  dpi = 300,
  width = 6.27,
  height = 4,
  units = "in"
)

###########################################
# 8. Density plots
###########################################

# Combine chains

phi_chains <- cbind(
  chain_1$theta_chain[, 1], chain_2$theta_chain[, 1],
  chain_3$theta_chain[, 1], chain_4$theta_chain[, 1]
)

sigma_x_chains <- cbind(
  chain_1$theta_chain[, 2], chain_2$theta_chain[, 2],
  chain_3$theta_chain[, 2], chain_4$theta_chain[, 2]
)

sigma_y_chains <- cbind(
  chain_1$theta_chain[, 3], chain_2$theta_chain[, 3],
  chain_3$theta_chain[, 3], chain_4$theta_chain[, 3]
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
  theme_minimal() +
  labs(
    x = "Parameter Value",
    y = "Density"
  )

ggsave("density_plot_phi.png",
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
  theme_minimal() +
  labs(
    x = "Parameter Value",
    y = "Density"
  )

ggsave("density_plot_sigma_x.png",
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
  theme_minimal() +
  labs(
    x = "Parameter Value",
    y = "Density"
  )

ggsave("density_plot_sigma_y.png",
  dpi = 300,
  width = 6.27,
  height = 4,
  units = "in"
)

###########################################
# 9. Inference
###########################################

print(result)

###########################################
# 10. Repeat the above many times
###########################################

main_function <- function(n_rep) {
  rmse_list <- numeric(n_rep)


  init_fn_ssm <- function(particles) {
    rnorm(particles, mean = 0, sd = 1)
  }

  transition_fn_ssm <- function(particles, phi, sigma_x) {
    phi * particles + sin(particles) +
      rnorm(length(particles), mean = 0, sd = sigma_x)
  }

  log_likelihood_fn_ssm <- function(y, particles, sigma_y) {
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

  t_val <- 50
  phi_true <- 0.7
  sigma_x_true <- 1
  sigma_y_true <- 1

  result_list <- vector("list", n_rep)

  for (i in 1:n_rep) {
    print(paste("Iteration:", i))

    set.seed(1405 + i)

    # Simulate data
    sim_data <- simulate_ssm(t_val, phi_true, sigma_x_true, sigma_y_true)
    x_true <- sim_data$x
    y_obs <- sim_data$y
    init_phi <- rnorm(1, mean = 0, sd = 1)
    init_sigma_x <- extraDistr::rhnorm(1)
    init_sigma_y <- extraDistr::rhnorm(1)

    init_params <- c(
      phi = init_phi,
      sigma_x = init_sigma_x,
      sigma_y = init_sigma_y
    )
    # Run PMMH
    result <- bayesSSM::pmmh(
      y = y_obs,
      m = 15000,
      init_fn_ssm = init_fn_ssm,
      transition_fn_ssm = transition_fn_ssm,
      log_likelihood_fn_ssm = log_likelihood_fn_ssm,
      log_priors = log_priors,
      init_params = init_params,
      burn_in = 2000,
      num_chains = 4,
      algorithm = "SISAR",
      resample_fn = "stratified",
      param_transform = c("identity", "log", "log"),
      verbose = TRUE,
      seed = 1405 + i
    )

    latent_est <- result$latent_state_estimate$mean
    rmse <- sqrt(mean((x_true - latent_est)^2))
    rmse_list[i] <- rmse

    result_list[[i]] <- result
  }
  list(rmse = rmse_list, result = result_list)
}

res <- main_function(4)

print(paste0("Mean RMSE: ", mean(res$rmse)))
print(paste0("SD RMSE: ", sd(res$rmse)))
