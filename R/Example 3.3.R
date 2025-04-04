library(ggplot2)
library(tibble)
library(tidyr)
set.seed(1405)
###########################################
# 1. RESAMPLING FUNCTIONS
###########################################

# Multinomial resampling: sample indices according to the weights.
resample_multinomial <- function(particles, weights) {
  n <- length(weights)
  indices <- sample(1:n, size = n, replace = TRUE, prob = weights)

  particles[indices]
}

# Stratified resampling: divides [0,1] into n strata and samples one point
# per stratum.
resample_stratified <- function(particles, weights) {
  n <- length(weights)
  positions <- (runif(1) + 0:(n - 1)) / n
  cumulative_sum <- cumsum(weights)
  indices <- numeric(n)
  i <- 1
  j <- 1
  while (i <= n) {
    if (positions[i] < cumulative_sum[j]) {
      indices[i] <- j
      i <- i + 1
    } else {
      j <- j + 1
    }
  }

  particles[indices]
}

# Systematic resampling: similar to stratified sampling but uses one
# random start.
resample_systematic <- function(particles, weights) {
  n <- length(weights)
  u0 <- runif(1, 0, 1 / n)
  positions <- u0 + (0:(n - 1)) / n
  cumulative_sum <- cumsum(weights)
  indices <- numeric(n)
  j <- 1
  for (i in 1:n) {
    while (positions[i] > cumulative_sum[j]) {
      j <- j + 1
    }
    indices[i] <- j
  }

  particles[indices]
}

###########################################
# 2. PARTICLE FILTER FUNCTION WITH FILTERING WEIGHTS HISTORY
###########################################

particle_filter <- function(
  y, n, init_fn, transition_fn, likelihood_fn,
  algorithm = c("SIS", "SISR", "SISAR"),
  resample_fn = c("multinomial", "stratified", "systematic"),
  threshold = NULL, return_particles = TRUE, ...
) {

  algorithm <- match.arg(algorithm)
  resample_fn <- match.arg(resample_fn)

  if (resample_fn == "multinomial") {
    resample_func <- resample_multinomial
  } else if (resample_fn == "stratified") {
    resample_func <- resample_stratified
  } else if (resample_fn == "systematic") {
    resample_func <- resample_systematic
  }

  t_len <- length(y)
  state_est <- numeric(t_len)
  ess_vec <- numeric(t_len)
  if (return_particles) {
    particles_history <- matrix(NA, nrow = t_len, ncol = n)
  }
  # NEW: store the filtering weights at each time step.
  weights_history <- matrix(NA, nrow = t_len, ncol = n)

  # Initialize particles and weights at time 1.
  particles <- init_fn(n, ...)
  weights <- likelihood_fn(y[1], particles, t = 1, ...)
  weights <- weights / sum(weights)

  state_est[1] <- sum(particles * weights)
  ess_vec[1] <- 1 / sum(weights^2)
  if (return_particles) {
    particles_history[1, ] <- particles
  }
  weights_history[1, ] <- weights

  # Filtering recursion
  for (t in 2:t_len) {
    particles <- transition_fn(particles, t, ...)
    likelihoods <- likelihood_fn(y[t], particles, t, ...)
    weights <- weights * likelihoods
    weights <- weights / sum(weights)

    ess_current <- 1 / sum(weights^2)
    ess_vec[t] <- ess_current

    if (algorithm == "SISR") {
      particles <- resample_func(particles, weights)
      weights <- rep(1 / n, n)
      ess_vec[t] <- n # reset ESS after resampling
    } else if (algorithm == "SISAR") {
      if (is.null(threshold)) {
        threshold <- n / 2
      }
      if (ess_current < threshold) {
        particles <- resample_func(particles, weights)
        weights <- rep(1 / n, n)
        ess_vec[t] <- n
      }
    }
    state_est[t] <- sum(particles * weights)

    if (return_particles) {
      particles_history[t, ] <- particles
    }
    weights_history[t, ] <- weights
  }

  result <- list(state_est = state_est, ESS = ess_vec, algorithm = algorithm)
  if (return_particles) {
    result$particles_history <- particles_history
    result$weights_history <- weights_history
  }
  return(result)
}

###########################################
# 3. STATE SPACE MODEL FUNCTIONS FOR SSM
###########################################
# SSM definitions:
#   X_1 ~ N(0,1)
#   X_t = phi * X_{t-1} + sin(X_{t-1}) + sigma_x * V_t,   V_t ~ N(0,1)
#   Y_t = X_t + sigma_y * W_t,              W_t ~ N(0,1)

init_fn_ssm <- function(n, ...) {
  rnorm(n, mean = 0, sd = 1)
}

transition_fn_ssm <- function(particles, t, phi, sigma_x, ...) {
  phi * particles + sin(particles) +
    rnorm(length(particles), mean = 0, sd = sigma_x)
}

likelihood_fn_ssm <- function(y, particles, t, sigma_y, ...) {
  dnorm(y, mean = particles, sd = sigma_y)
}

simulate_ssm <- function(t, phi, sigma_x, sigma_y) {
  x <- numeric(t)
  y <- numeric(t)
  x[1] <- rnorm(1, mean = 0, sd = 1)
  y[1] <- rnorm(1, mean = x[1], sd = sigma_y)
  for (t in 2:t) {
    x[t] <- phi * x[t - 1] + sin(x[t - 1]) + rnorm(1, mean = 0, sd = sigma_x)
    y[t] <- x[t] + rnorm(1, mean = 0, sd = sigma_y)
  }
  list(x = x, y = y)
}

###########################################
# 4. TRANSITION DENSITY FUNCTION FOR SMOOTHING
###########################################
transition_density_ssm <- function(x_current, x_prev, phi, sigma_x) {
  # Density for: X_t = phi * x_prev + sin(x_prev) + N(0, sigma_x^2)
  dnorm(x_current, mean = phi * x_prev + sin(x_prev), sd = sigma_x)
}

###########################################
# 5. BACKWARD SMOOTHING FUNCTION (FFBSm)
###########################################
backward_smoothing <- function(
  particles_history, weights_history, phi, sigma_x
) {
  t_len <- nrow(particles_history)
  n <- ncol(particles_history)

  # Initialize smoothing weights at time t as the filtering weights at t.
  smooth_weights <- matrix(0, nrow = t_len, ncol = n)
  smooth_weights[t_len, ] <- weights_history[t_len, ]

  # Backward recursion: for t = t-1 down to 1.
  for (t in (t_len - 1):1) {
    for (i in 1:n) {
      sum_val <- 0
      for (j in 1:n) {
        trans_prob <- transition_density_ssm(
          particles_history[t + 1, j],
          particles_history[t, i],
          phi, sigma_x
        )
        sum_val <- sum_val + trans_prob * smooth_weights[t + 1, j]
      }
      smooth_weights[t, i] <- weights_history[t, i] * sum_val
    }
    # Normalize the smoothing weights at time t.
    smooth_weights[t, ] <- smooth_weights[t, ] / sum(smooth_weights[t, ])
  }

  # Compute the smoothed state estimate for each time point.
  smoothed_est <- numeric(t_len)
  for (t in 1:t_len) {
    smoothed_est[t] <- sum(particles_history[t, ] * smooth_weights[t, ])
  }

  list(smoothed_est = smoothed_est, smooth_weights = smooth_weights)
}

###########################################
# 6. PARAMETERS AND SIMULATION
###########################################

t_val <- 50
n_particles <- 1000
phi_val <- 0.7
sigma_x_val <- 1
sigma_y_val <- 1

# Simulate data from the SSM
sim_data <- simulate_ssm(t_val, phi_val, sigma_x_val, sigma_y_val)
x_true <- sim_data$x
y_obs <- sim_data$y

# Plot the true state and observations
df <- tibble(
  time = seq_along(x_true),
  x_true = x_true,
  y_obs = y_obs
)

ggplot(df, aes(x = time)) +
  geom_line(aes(y = x_true, color = "True State"), linewidth = 1) +
  geom_point(aes(y = x_true, color = "True State"), size = 2, alpha = 0.7) +
  geom_point(aes(y = y_obs, color = "Observed"), size = 2, alpha = 0.7) +
  scale_color_manual(values = c("True State" = "blue", "Observed" = "red")) +
  labs(x = "Time", y = "Value", color = "Legend") +
  theme_minimal()

###########################################
# 7. RUN PARTICLE FILTER AND BACKWARD SMOOTHING
###########################################
# Run the particle filter
# (here, using the SISR algorithm with stratified resampling)
result_sisar <- particle_filter(
  y = y_obs, n = n_particles,
  init_fn = init_fn_ssm,
  transition_fn = transition_fn_ssm,
  likelihood_fn = likelihood_fn_ssm,
  algorithm = "SISAR",
  resample_fn = "stratified",
  phi = phi_val, sigma_x = sigma_x_val, sigma_y = sigma_y_val
)

# Run backward smoothing using the stored particle and weight histories
smooth_result <- backward_smoothing(result_sisar$particles_history,
  result_sisar$weights_history,
  phi = phi_val, sigma_x = sigma_x_val
)

###########################################
# 8. PLOT FILTERING AND SMOOTHED ESTIMATES VS TRUE STATE
###########################################
time <- 1:t_val
df_state <- data.frame(
  time = time,
  Latent_state = x_true,
  Filter_Est = result_sisar$state_est,
  Smooth_Est = smooth_result$smoothed_est
)

# Pivot data to long format
df_long <- pivot_longer(
  df_state,
  cols = -time,
  names_to = "Method",
  values_to = "State"
)

df_long$Method <- recode(
  df_long$Method,
  "Latent_state" = "Latent State",
  "Filter_Est" = "Filtering",
  "Smooth_Est" = "Smoothing"
)

ggplot(df_long, aes(x = time, y = State, color = Method, linetype = Method)) +
  geom_line(size = 1.2, alpha = 0.8) +
  labs(title = "Latent State, Filtering & Smoothing Estimates", y = "State") +
  scale_color_viridis_d(
    name = "Estimation Method",
    option = "plasma"
  ) +
  scale_linetype_manual(
    values = c("solid", "solid", "solid")
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.position = "top",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    panel.grid.minor = element_line(color = "gray90", linetype = "dotted")
  ) +
  guides(linetype = "none")

ggsave("Example_3.3_estimates.png", width = 6.27, height = 4, dpi = 300)



rmse_filtering <- sqrt(mean((result_sisar$state_est - x_true)^2))
rmse_smoothed <- sqrt(mean((smooth_result$smoothed_est - x_true)^2))

cat("RMSE for Filtering Estimates: ", rmse_filtering, "\n")
cat("RMSE for Smoothed Estimates: ", rmse_smoothed, "\n")


###########################################
# 9. Repeat many times
###########################################

# Function to compute RMSE
to_rmse <- function(true_vals, estimates) {
  sqrt(mean((true_vals - estimates)^2))
}

library(foreach)
library(doParallel)
library(progressr)
library(parallel)

# Number of simulations to run in parallel (only ran 10 times). It is ran
# 10.000 times in Julia, see Example 3.3.jl for details.

n_simulations <- 10

# Use all available cores minus one (to leave one core free)
n_cores <- parallel::detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Enable progress bar
handlers(global = TRUE)
handlers("txtprogressbar")

rmse_values <- with_progress({
  p <- progressor(along = 1:n_simulations)

  foreach(
    i = 1:n_simulations,
    .combine = c,
    .export = c(
      "simulate_ssm", "init_fn_ssm", "transition_fn_ssm",
      "likelihood_fn_ssm", "particle_filter", "backward_smoothing",
      "to_rmse", "phi_val", "sigma_x_val", "sigma_y_val",
      "t_val", "n_particles"
    )
  ) %dopar% {
    # Simulate data
    sim_data <- simulate_ssm(t_val, phi_val, sigma_x_val, sigma_y_val)
    x_true <- sim_data$x
    y_obs <- sim_data$y

    # Run the particle filter
    result_sisar <- particle_filter(
      y = y_obs, n = n_particles,
      init_fn = init_fn_ssm,
      transition_fn = transition_fn_ssm,
      likelihood_fn = likelihood_fn_ssm,
      algorithm = "SISAR",
      resample_fn = "stratified",
      phi = phi_val, sigma_x = sigma_x_val, sigma_y = sigma_y_val
    )

    # Run backward smoothing
    smooth_result <- backward_smoothing(result_sisar$particles_history,
      result_sisar$weights_history,
      phi = phi_val, sigma_x = sigma_x_val
    )

    # Update progress bar
    p()

    # Compute and return the RMSE for this simulation
    to_rmse(x_true, smooth_result$smoothed_est)
  }
})

stopCluster(cl)

mean_rmse <- mean(rmse_values)
sd_rmse <- sd(rmse_values)

cat("Mean RMSE:", mean_rmse, "\n")
cat("SD RMSE:", sd_rmse, "\n")
