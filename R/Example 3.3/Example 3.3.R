library(ggplot2)
library(tibble)
library(tidyr)
library(dplyr)
library(bayesSSM)
set.seed(1405)

###########################################
# 1. EXAMPLE SETUP: NON-LINEAR GAUSSIAN SSM
###########################################
# SSM definition:
#    X_0 ~ N(0,1)
#    X_t = phi * X_{t-1} + sin(X_{t-1}) + sigma_x * V_t,   V_t ~ N(0,1), t>=1
#    Y_t = X_t + sigma_y * W_t,              W_t ~ N(0,1), t>=1

init_fn <- function(num_particles) {
  rnorm(num_particles, mean = 0, sd = 1)
}

transition_fn <- function(particles, t, phi, sigma_x) {
  phi * particles + sin(particles) +
    rnorm(length(particles), mean = 0, sd = sigma_x)
}

log_likelihood_fn <- function(y, particles, t, sigma_y) {
  dnorm(y, mean = particles, sd = sigma_y, log = TRUE)
}

simulate_ssm <- function(t_val, phi, sigma_x, sigma_y) {
  init_state <- rnorm(1, mean = 0, sd = sigma_x)
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
# 2. PARAMETERS
###########################################

t_val <- 50
n_particles <- 1000
phi_val <- 0.7
sigma_x_val <- 1
sigma_y_val <- 1

###########################################
# 3. ONE SIMULATION
###########################################

sim_data <- simulate_ssm(t_val, phi_val, sigma_x_val, sigma_y_val)
x_true <- sim_data$x
y_obs  <- sim_data$y

###########################################
# 4. TRANSITION DENSITY FUNCTION FOR SMOOTHING
###########################################
transition_density_ssm <- function(x_current, x_prev, phi, sigma_x) {
  # Density for: X_t = phi *  X_{t-1} + sin(x_prev) + N(0, sigma_x^2)
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
# 6. RUN PARTICLE FILTER AND BACKWARD SMOOTHING
###########################################
result_sisar <- particle_filter(
  y = y_obs,
  n = n_particles,
  init_fn = init_fn,
  transition_fn = transition_fn,
  log_likelihood_fn = log_likelihood_fn,
  algorithm = "SISAR",
  resample_fn = "stratified",
  phi = phi_val,
  sigma_x = sigma_x_val,
  sigma_y = sigma_y_val
)

smooth_result <- backward_smoothing(
  particles_history = result_sisar$particles_history,
  weights_history   = result_sisar$weights_history,
  phi               = phi_val,
  sigma_x           = sigma_x_val
)



###########################################
# 7. PLOT FILTERING AND SMOOTHED ESTIMATES VS TRUE STATE
###########################################
time <- 0:t_val
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
  geom_line(linewidth = 1.2, alpha = 0.8) +
  labs(
    title = "Latent State, Filtering & Smoothing Estimates",
    x = "Time",
    y = "State"
  ) +
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
    panel.grid = element_blank()
  ) +
  guides(linetype = "none")

ggsave("Example_3.3_estimates.png", width = 6.27, height = 4, dpi = 300)



rmse_filtering <- sqrt(mean((result_sisar$state_est - x_true)^2))
rmse_smoothed <- sqrt(mean((smooth_result$smoothed_est - x_true)^2))

cat("RMSE for Filtering Estimates: ", rmse_filtering, "\n")
cat("RMSE for Smoothed Estimates: ", rmse_smoothed, "\n")


###########################################
# 8. Repeat many times
###########################################

# It is ran 10.000 times in Julia, see Example 3.3.jl for details.

