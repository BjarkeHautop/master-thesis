library(outbreaks)
library(bayesSSM) # My R package see https://bjarkehautop.github.io/bayesSSM/
library(ggplot2)
library(msm)
library(dplyr)
library(ggdist)


set.seed(1405)

# Load the data
outbreak_data <- outbreaks::influenza_england_1978_school

# Add a time variable and create the preprocessed data with day 0 included
outbreak_data$time <- seq_len(nrow(outbreak_data))
preprocessed_data <- outbreak_data[, c("time", "in_bed")]

ggplot(preprocessed_data, aes(x = time, y = in_bed)) +
  geom_point(color = "#D55E00", size = 3) +
  labs(
    title = "Influenza Outbreak in Boarding School",
    x = "Time (Days)",
    y = "Number of Students in Bed"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none"
  ) +
  scale_x_continuous(breaks = 1:14)

ggsave(
  "influenza_outbreak_plot.png",
  dpi = 300,
  width = 6.27,
  height = 4,
  units = "in"
)

# Define variables

n_total <- 763
init_infected <- 1
init_state <- c(S = n_total - init_infected, I = init_infected)
t_max <- nrow(preprocessed_data)
observations <- preprocessed_data$in_bed


###########################################
# 1. Epidemic Simulation
###########################################

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

###########################################
# 2. Priors and Model Specification
###########################################

# --- Priors ---
# Mean infection rate 2 days (1/lambda = 2 days) and P(lambda > 1) = 0.2
# Solving for the parameters using my R package
# https://github.com/BjarkeHautop/solvetruncated
library(solvetruncated)
solvetruncated::solve_truncated_normal(0.5, 1, 0.8)

# Define the log-prior for the parameters
log_prior_lambda <- function(lambda) {
  msm::dtnorm(lambda, mean = 0, sd = 0.63, lower = 0, log = TRUE)
}

# Mean infection period 3 days (1/gamma = 3 days) and P(gamma > 1) = 0.1.
solvetruncated::solve_truncated_normal(0.33, 1, 0.9)

log_prior_gamma <- function(gamma) {
  msm::dtnorm(gamma, mean = 0, sd = 0.41, log = TRUE)
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
  n_total <- 763
  init_infected <- 1
  init_state <- c(S = n_total - init_infected, I = init_infected)
  matrix(rep(init_state, each = num_particles), nrow = num_particles, byrow = FALSE)
}

transition_fn_epidemic <- function(particles, lambda, gamma, t) {

  n_total <- 763

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
# 3. Prior Predictive Check
###########################################
sample_lambda <- function(n) {
  msm::rtnorm(n, mean = 0, sd = 0.63, lower = 0)
}

sample_gamma <- function(n) {
  msm::rtnorm(n, mean = 0, sd = 0.41, lower = 0)
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

# Simulate trajectories
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
      length(states[, 2]),
      size = phi, mu = states[, 2]
    )

    prior_trajectories[[i]] <- noisy_infected
  }

  prior_trajectories
}

# Generate the prior predictive trajectories
prior_trajectories <- simulate_trajectory(
  prior_samples, n_total, init_infected, t_max
)

# Prepare the data for plotting
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
    linewidth = 1,
    linetype = "dashed"
  ) +
  geom_line(
    data = plot_data[plot_data$type == "Observed", ],
    aes(x = time, y = infected),
    color = "darkblue",
    linewidth = 1
  ) +
  labs(
    x = "Time (Days)",
    y = "Observed Infected",
    title = "Prior Predictive Check"
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

# Save the plot
ggsave(
  "prior_predictive_check_influenza_sir_negbin.png",
  dpi = 300,
  width = 6.27,
  height = 4,
  units = "in"
)

###########################################
# 4. PMMH
###########################################

pilot_init_params <- split(prior_samples, seq_len(nrow(prior_samples)))
pilot_init_params <- lapply(
  pilot_init_params,
  function(row) setNames(as.numeric(row), names(prior_samples))
)

# Use 1st 4 rows of prior_samples as initial values for the pilot chain
pilot_init_params <- pilot_init_params[1:4]
pilot_init_params

result_influenza_sir_negbin <- pmmh(
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

saveRDS(result_influenza_sir_negbin, file = "result_influenza_sir_negbin.rds")

# Or load the result from the saved file
result_influenza_sir_negbin <- readRDS("result_influenza_sir_negbin.rds")

chains <- result_influenza_sir_negbin$theta_chain

ggplot(chains, aes(x = lambda, fill = chain)) +
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
  "density_plot_lambda_influenza_sir_negbin.png",
  dpi = 300,
  width = 6.27,
  height = 4,
  units = "in"
)

ggplot(chains, aes(x = gamma, fill = chain)) +
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
  "density_plot_gamma_influenza_sir_negbin.png",
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
  "density_plot_phi_influenza_sir_negbin.png",
  dpi = 300,
  width = 6.27,
  height = 4,
  units = "in"
)

### Combined plots

ggplot(chains, aes(x = lambda, y = 0)) +
  stat_halfeye(
    point_interval = "mean_qi",
    .width = 0.95,
    slab_fill         = "#3B83BD80",
    slab_color        = "NA",
    slab_alpha        = 0.6,
    interval_size     = 1.2,
    point_size        = 2.5,
    point_color       = "black",
    interval_color    = "black"
  ) +
  theme_bw() +
  theme(
    panel.border        = element_rect(linewidth = 0.8, colour = "black"),
    panel.grid.major    = element_line(linewidth = 0.4, colour = "grey85"),
    panel.grid.minor    = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    plot.title          = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Posterior Density of Transmission Rate λ",
    x = "Parameter Value",
    y = "Density"
  )

ggsave(
  "combined_density_lambda_influenza_sir_negbin.png",
  dpi = 300,
  width = 6.27,
  height = 4,
  units = "in"
)

ggplot(chains, aes(x = gamma, y = 0)) +
  stat_halfeye(
    point_interval = "mean_qi",
    .width = 0.95,
    slab_fill         = "#3B83BD80",
    slab_color        = "NA",
    slab_alpha        = 0.6,
    interval_size     = 1.2,
    point_size        = 2.5,
    point_color       = "black",
    interval_color    = "black"
  ) +
  theme_bw() +
  theme(
    panel.border        = element_rect(linewidth = 0.8, colour = "black"),
    panel.grid.major    = element_line(linewidth = 0.4, colour = "grey85"),
    panel.grid.minor    = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    plot.title          = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Posterior Density of Recovery Rate γ",
    x = "Parameter Value",
    y = "Density"
  )

ggsave(
  "combined_density_gamma_influenza_sir_negbin.png",
  dpi = 300,
  width = 6.27,
  height = 4,
  units = "in"
)

###########################################
# 5. Inference on R_0 and recovery_time
###########################################

chains$r0 <- chains$lambda / chains$gamma

r0_mean <- mean(chains$r0)

r0_ci <- quantile(chains$r0, probs = c(0.025, 0.975))
r0_matrix <- matrix(chains$r0, nrow = nrow(chains) / 4, ncol = 4)

ess(r0_matrix)
rhat(r0_matrix)

r0_mean
r0_ci

chains$recovery_time <- 1 / chains$gamma

recovery_time_mean <- mean(chains$recovery_time)

recovery_time_ci <- quantile(chains$recovery_time, probs = c(0.025, 0.975))
recovery_time_matrix <- matrix(
  chains$recovery_time,
  nrow = nrow(chains) / 4, ncol = 4
)

ess(recovery_time_matrix)
rhat(recovery_time_matrix)

recovery_time_mean
recovery_time_ci


# Combined density plot for R0 and recovery time

ggplot(chains, aes(x = r0, y = 0)) +
  stat_halfeye(
    point_interval = "mean_qi",
    .width = 0.95,
    slab_fill         = "#3B83BD80",
    slab_color        = "NA",
    slab_alpha        = 0.6,
    interval_size     = 1.2,
    point_size        = 2.5,
    point_color       = "black",
    interval_color    = "black"
  ) +
  theme_bw() +
  theme(
    panel.border        = element_rect(linewidth = 0.8, colour = "black"),
    panel.grid.major    = element_line(linewidth = 0.4, colour = "grey85"),
    panel.grid.minor    = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    plot.title          = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Posterior Density of Basic Reproduction Number R₀",
    x = "Parameter Value",
    y = "Density"
  )

ggsave(
  "combined_density_r0_influenza_sir_negbin.png",
  dpi = 300,
  width = 6.27,
  height = 4,
  units = "in"
)

ggplot(chains, aes(x = recovery_time, y = 0)) +
  stat_halfeye(
    point_interval = "mean_qi",
    .width = 0.95,
    slab_fill         = "#3B83BD80",
    slab_color        = "NA",
    slab_alpha        = 0.6,
    interval_size     = 1.2,
    point_size        = 2.5,
    point_color       = "black",
    interval_color    = "black"
  ) +
  theme_bw() +
  theme(
    panel.border        = element_rect(linewidth = 0.8, colour = "black"),
    panel.grid.major    = element_line(linewidth = 0.4, colour = "grey85"),
    panel.grid.minor    = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    plot.title          = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Posterior Density of Mean Recovery Time 1/γ",
    x = "Parameter Value",
    y = "Density"
  )

ggsave(
  "combined_density_recovery_time_influenza_sir_negbin.png",
  dpi = 300,
  width = 6.27,
  height = 4,
  units = "in"
)


###########################################
# 6. Posterior Predictive Check
###########################################

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
  "posterior_predictive_check_influenza_sir_negbin.png",
  dpi = 300,
  width = 6.27,
  height = 4,
  units = "in"
)

###############################################
# 7. Estimating latent state (number of infected)
###############################################

s_i_est <- result_influenza_sir_negbin$latent_state_chain
array_s_i_est <- simplify2array(lapply(s_i_est, function(x) simplify2array(x)))

# Compute quantiles for each entry in the 15x2 matrix
quantiles_s_i <- apply(
  array_s_i_est, c(1, 2),
  function(x) quantile(x, probs = c(0.025, 0.5, 0.975))
)

# Convert to a data frame for plotting
quantiles_s_i <- as.data.frame(quantiles_s_i)

# Select infected quantiles
infected_quantiles <- quantiles_s_i[17:30]
infected_quantiles <- t(infected_quantiles)

infected_quantiles_df <- as.data.frame(infected_quantiles)
infected_quantiles_df$time <- 1:(nrow(infected_quantiles_df))

rownames(infected_quantiles_df) <- NULL

# Remove % in the col names

colnames(infected_quantiles_df) <- make.names(colnames(infected_quantiles_df))

# Make ggplot
ggplot(infected_quantiles_df, mapping = aes(x = time)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "green", alpha = 0.35) +
  geom_line(mapping = aes(x = time, y = X50.), color = "blue", linewidth = 1) +
  geom_point(
    data = preprocessed_data,
    aes(x = time, y = in_bed),
    color = "#D55E00",
    size = 3,
    shape = 16
  ) +
  labs(
    x = "Day",
    y = "Number of Infected",
    title = "Estimated and Observed Number of Infected",
  ) +
  scale_x_continuous(
    breaks = seq(1, 14, by = 1)
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

ggsave(
  "latent_state_estimation_influenza_sir_negbin.png",
  dpi = 300,
  width = 6.27,
  height = 4,
  units = "in"
)

###############################################
# 8. Prior sensitivity analysis
###############################################

# Old priors
log_prior_lambda_old <- function(lambda) {
  msm::dtnorm(lambda, mean = 0, sd = 0.63, lower = 0, log = TRUE)
}

log_prior_gamma_old <- function(gamma) {
  msm::dtnorm(gamma, mean = 0, sd = 0.41, lower = 0, log = TRUE)
}

# New priors
log_prior_lambda_new <- function(lambda) {
  msm::dtnorm(lambda, mean = 0, sd = 1, lower = 0, log = TRUE)
}

log_prior_gamma_new <- function(gamma) {
  msm::dtnorm(gamma, mean = 0, sd = 1, lower = 0, log = TRUE)
}

# Function to compute importance weights
importance_weights <- function(lambda, gamma) {
  log_w <- log_prior_lambda_new(lambda) + log_prior_gamma_new(gamma) -
    log_prior_lambda_old(lambda) - log_prior_gamma_old(gamma)
  exp(log_w) # Convert to normal scale
}

# Extract posterior samples
lambda_samples <- result_influenza_sir_negbin$theta_chain$lambda
gamma_samples <- result_influenza_sir_negbin$theta_chain$gamma
phi_samples <- result_influenza_sir_negbin$theta_chain$phi
r0_samples <- lambda_samples / gamma_samples
recovery_time_samples <- 1 / gamma_samples

# Compute importance sampling weights
weights <- importance_weights(lambda_samples, gamma_samples)

# Normalize weights
weights <- weights / sum(weights)

# Compute weighted mean and credible intervals
weighted_mean_lambda <- sum(weights * lambda_samples)
weighted_mean_gamma <- sum(weights * gamma_samples)
weighted_mean_r0 <- sum(weights * r0_samples)
weighted_mean_recovery_time <- sum(weights * recovery_time_samples)


weighted_quantile <- function(x, w, probs) {
  ord <- order(x)
  x_sorted <- x[ord]
  w_sorted <- w[ord]
  cum_w <- cumsum(w_sorted) / sum(w_sorted)
  approx(cum_w, x_sorted, xout = probs, method = "linear", rule = 2)$y
}

# Compute weighted credible intervals
credible_interval_lambda <- weighted_quantile(
  lambda_samples, weights, c(0.025, 0.975)
)
credible_interval_gamma <- weighted_quantile(
  gamma_samples, weights, c(0.025, 0.975)
)
credible_interval_r0 <- weighted_quantile(
  r0_samples, weights, c(0.025, 0.975)
)
credible_interval_recovery_time <- weighted_quantile(
  recovery_time_samples, weights, c(0.025, 0.975)
)

weighted_mean_lambda
weighted_mean_gamma
weighted_mean_r0
weighted_mean_recovery_time

credible_interval_lambda
credible_interval_gamma
credible_interval_r0
credible_interval_recovery_time


##############################################
# 9. MLE using iterated filtering 2
##############################################

# Stratified resampling function
resample_stratified <- function(indices, weights) {
  n <- length(indices)
  positions <- (runif(n) + seq_len(n) - 1) / n
  cumulative_sum <- cumsum(weights)
  indices_resampled <- numeric(n)
  i <- 1
  j <- 1
  while (i <= n) {
    if (positions[i] < cumulative_sum[j]) {
      indices_resampled[i] <- j
      i <- i + 1
    } else {
      j <- j + 1
    }
  }
  indices[indices_resampled]
}


if2_epidemic <- function(
    y,
    start_params,
    Np = 500,
    Nmif = 100,
    rw_sd = c(lambda = 0.1, gamma = 0.1, phi = 0.1),
    cooling_fraction = 0.95,
    return_trace = FALSE
) {
  p <- length(start_params)
  params <- matrix(start_params, nrow = Np, ncol = p, byrow = TRUE)
  colnames(params) <- names(start_params)

  mean_params <- matrix(NA, nrow = Nmif, ncol = p)
  colnames(mean_params) <- names(start_params)

  for (m in seq_len(Nmif)) {
    if (m %% 10 == 0) cat("IF2 iteration", m, "of", Nmif, "\n")
    sd_m <- rw_sd * cooling_fraction^(m - 1)

    # Add noise to parameter swarm
    Theta_F <- params + matrix(
      rnorm(Np * p, mean = 0, sd = sd_m), nrow = Np, ncol = p
    )
    colnames(Theta_F) <- names(start_params)
    X_F <- init_fn_epidemic(Np)

    for (n in seq_along(y)) {
      Theta_P <- Theta_F + matrix(
        rnorm(Np * p, mean = 0, sd = sd_m), nrow = Np, ncol = p
      )
      colnames(Theta_P) <- names(start_params)

      X_P <- matrix(NA, nrow = Np, ncol = ncol(X_F))
      for (i in seq_len(Np)) {
        X_P[i, ] <- transition_fn_epidemic(
          particles = matrix(X_F[i, ], nrow = 1),
          lambda = Theta_P[i, "lambda"],
          gamma  = Theta_P[i, "gamma"],
          t = n
        )
      }

      log_w <- log_likelihood_fn_epidemic(
        y = y[n],
        particles = X_P,
        phi = Theta_P[,"phi"],
        t = n
      )

      w <- exp(log_w - max(log_w))
      w <- w / sum(w)
      idx <- resample_stratified(1:Np, w)
      Theta_F <- Theta_P[idx, , drop = FALSE]
      X_F     <- X_P[idx, , drop = FALSE]
    }

    mean_params[m, ] <- colMeans(Theta_F)
    params <- Theta_F
  }

  theta_est <- colMeans(params)
  names(theta_est) <- names(start_params)

  if (return_trace) {
    return(list(estimates = theta_est, trace = mean_params))
  } else {
    return(theta_est)
  }
}

set.seed(1405)
mle_obs <- if2_epidemic(
  y = observations,
  start_params = c(
    lambda = median(chains$lambda),
    gamma = median(chains$gamma),
    phi = median(chains$phi)
  ),
  Np = 200,
  Nmif = 500,
  rw_sd = c(lambda = 0.1, gamma = 0.1, phi = 0.1),
  cooling_fraction = 0.99,
  return_trace = FALSE
)
mle_obs


# Bootstrap
B <- 200
boot.mat <- matrix(NA, nrow = B, ncol = 3)
for (b in 1:B) {
  cat("Bootstrap iteration", b, "of", B, "\n")
  set.seed(1405 + b)
  latent_i <- rep(0, length(observations))

  # If disease dies out in the first week make a new sample
  while (any(latent_i[1:7] == 0)) {
    sim_states <- simulate_epidemic(
      n_total, init_infected, mle_obs[["lambda"]], mle_obs[["gamma"]], t_max
    )
    latent_i <- sim_states[, 2]
    yb <- rnbinom(length(latent_i), size = mle_obs[["phi"]], mu = latent_i)
  }

  fit <- if2_epidemic(
    y = yb,
    start_params = mle_obs,
    Np = 500,
    Nmif = 100,
    rw_sd = c(lambda = 0.1, gamma = 0.1, phi = 0.1),
    cooling_fraction = 0.95
  )
  boot.mat[b, ] <- fit
  print(boot.mat[b, ])
}

# Save as rds
saveRDS(boot.mat, file = "boot_mat.rds")

# Load the boot.mat if needed
boot.mat <- readRDS("boot_mat.rds")
boot.mat

# Compute the mean and 95% CI
lambda_ci <- quantile(boot.mat[, 1], probs = c(0.025, 0.975))

gamma_ci <- quantile(boot.mat[, 2], probs = c(0.025, 0.975))

# Compute R0 and recovery time
boot.mat <- as.data.frame(boot.mat)
boot.mat$r0 <- boot.mat[, 1] / boot.mat[, 2]
boot.mat$recovery_time <- 1 / boot.mat[, 2]

r0_ci <- quantile(boot.mat$r0, probs = c(0.025, 0.975))

recovery_time_ci <- quantile(boot.mat$recovery_time, probs = c(0.025, 0.975))

# Print results
cat("95% CI for lambda:", lambda_ci, "\n")
cat("95% CI for gamma:", gamma_ci, "\n")
cat("95% CI for R0:", r0_ci, "\n")
cat("95% CI for recovery time:", recovery_time_ci, "\n")
