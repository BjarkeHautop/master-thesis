library(ggplot2)
library(tibble)
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
# 2. MAIN PARTICLE FILTER FUNCTION
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

  particles <- init_fn(n, ...)
  weights <- likelihood_fn(y[1], particles, t = 1, ...)
  weights <- weights / sum(weights)

  state_est[1] <- sum(particles * weights)
  ess_vec[1] <- 1 / sum(weights^2)
  if (return_particles) {
    particles_history[1, ] <- particles
  }

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
      ess_vec[t] <- n  # reset ESS after resampling
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
  }

  result <- list(state_est = state_est, ESS = ess_vec, algorithm = algorithm)
  if (return_particles) result$particles_history <- particles_history
  return(result)
}

###########################################
# 3. EXAMPLE SETUP: NON-LINEAR GAUSSIAN SSM
###########################################
# SSM definitions:
#    X_1 ~ N(0,1)
#    X_t = phi * X_{t-1} + sin(X_{t-1}) + sigma_x * V_t,   V_t ~ N(0,1)
#    Y_t = X_t + sigma_y * W_t,              W_t ~ N(0,1)

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
# 4. PARAMETERS
###########################################

t_val <- 50
n_particles <- 1000
phi_val <- 0.7
sigma_x_val <- 1
sigma_y_val <- 1

###########################################
# 5. ONE SIMULATION
###########################################

sim_data <- simulate_ssm(t_val, phi_val, sigma_x_val, sigma_y_val)
x_true <- sim_data$x
y_obs  <- sim_data$y

df <- tibble(
  time = seq_along(x_true),
  x_true = x_true,
  y_obs = y_obs
)

# Plot
ggplot(df, aes(x = time)) +
  geom_line(aes(y = x_true, color = "X_t"), linewidth = 1.2) +
  geom_point(aes(y = y_obs, color = "Y_t"), size = 1.8, alpha = 0.5) +
  scale_color_viridis_d(option = "plasma", begin = 0.2, end = 0.8) +
  labs(title = "Latent and Observed States Over Time",
       x = "Time", y = "State Value", color = "State Type") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "none",
    panel.grid = element_blank()
  )

ggsave("example_3.1.png", width = 6.27, height = 4, dpi = 300)


result_sis <- particle_filter(
  y = y_obs,
  n = n_particles,
  init_fn = init_fn_ssm,
  transition_fn = transition_fn_ssm,
  likelihood_fn = likelihood_fn_ssm,
  algorithm = "SIS",
  phi = phi_val,
  sigma_x = sigma_x_val,
  sigma_y = sigma_y_val
)

result_sisr <- particle_filter(
  y = y_obs,
  n = n_particles,
  init_fn = init_fn_ssm,
  transition_fn = transition_fn_ssm,
  likelihood_fn = likelihood_fn_ssm,
  algorithm = "SISR",
  resample_fn = "stratified",
  phi = phi_val,
  sigma_x = sigma_x_val,
  sigma_y = sigma_y_val
)

result_sisar <- particle_filter(
  y = y_obs,
  n = n_particles,
  init_fn = init_fn_ssm,
  transition_fn = transition_fn_ssm,
  likelihood_fn = likelihood_fn_ssm,
  algorithm = "SISAR",
  resample_fn = "stratified",
  phi = phi_val,
  sigma_x = sigma_x_val,
  sigma_y = sigma_y_val
)

rmse <- function(est, true) sqrt(mean((est - true)^2))

rmse_sis   <- rmse(result_sis$state_est, x_true)
rmse_sisr  <- rmse(result_sisr$state_est, x_true)
rmse_sisar <- rmse(result_sisar$state_est, x_true)

cat("Single simulation RMSEs:\n")
cat("  SIS   :", rmse_sis, "\n")
cat("  SISR  :", rmse_sisr, "\n")
cat("  SISAR :", rmse_sisar, "\n\n")

time <- 1:t_val
df_state <- data.frame(
  Time = rep(time, 4),
  State = c(x_true,
            result_sis$state_est,
            result_sisr$state_est,
            result_sisar$state_est),
  Method = factor(rep(c("Latent State", "SIS", "SISR", "SISAR"),
                      each = t_val),
                  levels = c("Latent State", "SIS", "SISR", "SISAR"))
)

ggplot(df_state, aes(x = Time, y = State, color = Method, linetype = Method)) +
  geom_line(size = 1.2, alpha = 0.8) +
  labs(title = "Latent States and Particle Filter Estimates", y = "State") +
  scale_color_viridis_d(option = "plasma") +
  scale_linetype_manual(values = c("solid", "solid", "solid", "solid")) +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 12),
    legend.position = "top",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

ggsave("example_3.1_estimates.png", width = 6.27, height = 4, dpi = 300)


###########################################
# 6. REPEATED SIMULATIONS
###########################################

n_reps <- 10000
rmse_sis_all   <- numeric(n_reps)
rmse_sisr_all  <- numeric(n_reps)
rmse_sisar_all <- numeric(n_reps)

for (i in 1:n_reps) {
  sim_data <- simulate_ssm(t_val, phi_val, sigma_x_val, sigma_y_val)
  x_true <- sim_data$x
  y_obs  <- sim_data$y

  result_sis <- particle_filter(
    y = y_obs,
    n = n_particles,
    init_fn = init_fn_ssm,
    transition_fn = transition_fn_ssm,
    likelihood_fn = likelihood_fn_ssm,
    algorithm = "SIS",
    phi = phi_val,
    sigma_x = sigma_x_val,
    sigma_y = sigma_y_val
  )

  result_sisr <- particle_filter(
    y = y_obs,
    n = n_particles,
    init_fn = init_fn_ssm,
    transition_fn = transition_fn_ssm,
    likelihood_fn = likelihood_fn_ssm,
    algorithm = "SISR",
    resample_fn = "stratified",
    phi = phi_val,
    sigma_x = sigma_x_val,
    sigma_y = sigma_y_val
  )

  result_sisar <- particle_filter(
    y = y_obs,
    n = n_particles,
    init_fn = init_fn_ssm,
    transition_fn = transition_fn_ssm,
    likelihood_fn = likelihood_fn_ssm,
    algorithm = "SISAR",
    resample_fn = "stratified",
    phi = phi_val,
    sigma_x = sigma_x_val,
    sigma_y = sigma_y_val
  )

  rmse_sis_all[i]   <- rmse(result_sis$state_est, x_true)
  rmse_sisr_all[i]  <- rmse(result_sisr$state_est, x_true)
  rmse_sisar_all[i] <- rmse(result_sisar$state_est, x_true)
}

rmse_sis_mean   <- mean(rmse_sis_all)
rmse_sis_sd     <- sd(rmse_sis_all)
rmse_sisr_mean  <- mean(rmse_sisr_all)
rmse_sisr_sd    <- sd(rmse_sisr_all)
rmse_sisar_mean <- mean(rmse_sisar_all)
rmse_sisar_sd   <- sd(rmse_sisar_all)

cat("RMSE over", n_reps, "replications:\n")
cat("  SIS   : Mean =", rmse_sis_mean, "SD =", rmse_sis_sd, "\n")
cat("  SISR  : Mean =", rmse_sisr_mean, "SD =", rmse_sisr_sd, "\n")
cat("  SISAR : Mean =", rmse_sisar_mean, "SD =", rmse_sisar_sd, "\n\n")

df_rmse <- data.frame(
  Method = factor(
    c("SIS", "SISR", "SISAR"),
    levels = c("SIS", "SISR", "SISAR")
  ),
  rmse_mean = c(rmse_sis_mean, rmse_sisr_mean, rmse_sisar_mean),
  rmse_sd   = c(rmse_sis_sd, rmse_sisr_sd, rmse_sisar_sd)
)

p_rmse <- ggplot(df_rmse, aes(x = Method, y = rmse_mean, fill = Method)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_errorbar(aes(ymin = rmse_mean - rmse_sd, ymax = rmse_mean + rmse_sd),
                width = 0.2) +
  theme_minimal() +
  labs(
    title = "RMSE of Particle Filter Methods (10000 Replications)",
    y = "RMSE"
  ) +
  scale_fill_manual(values = c("red", "blue", "green"))
print(p_rmse)
