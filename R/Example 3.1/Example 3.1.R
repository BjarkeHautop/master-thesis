library(ggplot2)
library(tibble)
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
y_obs <- sim_data$y

df <- tibble(
  time = seq_along(x_true[-1]),
  x_true = x_true[-1],
  y_obs = y_obs
)

# Plot
ggplot(df, aes(x = time)) +
  geom_line(aes(y = x_true, color = "X_t"), linewidth = 1) +
  geom_point(aes(y = y_obs, color = "Y_t"), size = 1.8, alpha = 0.5) +
  scale_color_viridis_d(option = "plasma", begin = 0.2, end = 0.8) +
  labs(
    title = "Latent and Observed States Over Time",
    x = "Time", y = "State Value", color = "State Type"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none",
    panel.grid = element_blank()
  )

ggsave("outputs/example_3.1.png", width = 6.27, height = 4, dpi = 300)

result_sis <- particle_filter(
  y = y_obs,
  n = n_particles,
  init_fn = init_fn,
  transition_fn = transition_fn,
  log_likelihood_fn = log_likelihood_fn,
  algorithm = "SIS",
  resample_fn = "stratified",
  phi = phi_val,
  sigma_x = sigma_x_val,
  sigma_y = sigma_y_val
)

result_sisr <- particle_filter(
  y = y_obs,
  n = n_particles,
  init_fn = init_fn,
  transition_fn = transition_fn,
  log_likelihood_fn = log_likelihood_fn,
  algorithm = "SISR",
  resample_fn = "stratified",
  phi = phi_val,
  sigma_x = sigma_x_val,
  sigma_y = sigma_y_val
)

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

rmse <- function(est, true) sqrt(mean((est - true)^2))

rmse_sis <- rmse(result_sis$state_est, x_true)
rmse_sisr <- rmse(result_sisr$state_est, x_true)
rmse_sisar <- rmse(result_sisar$state_est, x_true)

cat("Single simulation RMSEs (my package):\n")
cat("  SIS   :", rmse_sis, "\n")
cat("  SISR  :", rmse_sisr, "\n")
cat("  SISAR :", rmse_sisar, "\n\n")

time <- 0:t_val
df_state <- data.frame(
  Time = rep(time, 4),
  State = c(
    x_true,
    result_sis$state_est,
    result_sisr$state_est,
    result_sisar$state_est
  ),
  Method = factor(
    rep(c("Latent State", "SIS", "SISR", "SISAR"),
      each = length(time)
    ),
    levels = c("Latent State", "SIS", "SISR", "SISAR")
  )
)

ggplot(df_state, aes(x = Time, y = State, color = Method, linetype = Method)) +
  geom_line(linewidth = 1, alpha = 0.8) +
  labs(title = "Latent States and Particle Filter Estimates", y = "State Value") +
  scale_color_viridis_d(option = "plasma") +
  scale_linetype_manual(values = c("solid", "solid", "solid", "solid")) +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold"),
    legend.position = "top",
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid = element_blank()
  )

ggsave("example_3.1_estimates.png", width = 6.27, height = 4, dpi = 300)


###########################################
# 4. REPEATED SIMULATIONS
###########################################

n_reps <- 10000
rmse_sis_all <- numeric(n_reps)
rmse_sisr_all <- numeric(n_reps)
rmse_sisar_all <- numeric(n_reps)

for (i in 1:n_reps) {
  set.seed(1405 + i)
  if (i %% 500 == 0) {
    cat("Iteration", i, "of", n_reps, "\n")
  }
  sim_data <- simulate_ssm(t_val, phi_val, sigma_x_val, sigma_y_val)
  x_true <- sim_data$x
  y_obs <- sim_data$y

  result_sis <- particle_filter(
    y = y_obs,
    n = n_particles,
    init_fn = init_fn,
    transition_fn = transition_fn,
    log_likelihood_fn = log_likelihood_fn,
    algorithm = "SIS",
    phi = phi_val,
    sigma_x = sigma_x_val,
    sigma_y = sigma_y_val
  )

  result_sisr <- particle_filter(
    y = y_obs,
    n = n_particles,
    init_fn = init_fn,
    transition_fn = transition_fn,
    log_likelihood_fn = log_likelihood_fn,
    algorithm = "SISR",
    resample_fn = "stratified",
    phi = phi_val,
    sigma_x = sigma_x_val,
    sigma_y = sigma_y_val
  )

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

  rmse_sis_all[i] <- rmse(result_sis$state_est, x_true)
  rmse_sisr_all[i] <- rmse(result_sisr$state_est, x_true)
  rmse_sisar_all[i] <- rmse(result_sisar$state_est, x_true)
}

rmse_sis_mean <- mean(rmse_sis_all)
rmse_sis_sd <- sd(rmse_sis_all)
rmse_sisr_mean <- mean(rmse_sisr_all)
rmse_sisr_sd <- sd(rmse_sisr_all)
rmse_sisar_mean <- mean(rmse_sisar_all)
rmse_sisar_sd <- sd(rmse_sisar_all)

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
  rmse_sd = c(rmse_sis_sd, rmse_sisr_sd, rmse_sisar_sd)
)

p_rmse <- ggplot(df_rmse, aes(x = Method, y = rmse_mean, fill = Method)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_errorbar(aes(ymin = rmse_mean - rmse_sd, ymax = rmse_mean + rmse_sd),
    width = 0.2
  ) +
  theme_minimal() +
  labs(
    title = "RMSE of Particle Filter Methods (10000 Replications)",
    y = "RMSE"
  ) +
  scale_fill_manual(values = c("red", "blue", "green"))
print(p_rmse)
