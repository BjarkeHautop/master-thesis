library(ggplot2)
set.seed(1405)

###########################################
# 1. RESAMPLING FUNCTIONS
###########################################

# Multinomial resampling: sample indices according to the weights.
resample_multinomial <- function(particles, weights) {
  N <- length(weights)
  indices <- sample(1:N, size = N, replace = TRUE, prob = weights)
  return(particles[indices])
}

# Stratified resampling: divides [0,1] into N strata and samples one point per stratum.
resample_stratified <- function(particles, weights) {
  N <- length(weights)
  positions <- (runif(1) + 0:(N - 1)) / N
  cumulative_sum <- cumsum(weights)
  indices <- numeric(N)
  i <- 1; j <- 1
  while (i <= N) {
    if (positions[i] < cumulative_sum[j]) {
      indices[i] <- j
      i <- i + 1
    } else {
      j <- j + 1
    }
  }
  return(particles[indices])
}

# Systematic resampling: similar to stratified sampling but uses one random start.
resample_systematic <- function(particles, weights) {
  N <- length(weights)
  u0 <- runif(1, 0, 1/N)
  positions <- u0 + (0:(N - 1)) / N
  cumulative_sum <- cumsum(weights)
  indices <- numeric(N)
  j <- 1
  for (i in 1:N) {
    while (positions[i] > cumulative_sum[j]) {
      j <- j + 1
    }
    indices[i] <- j
  }
  return(particles[indices])
}

###########################################
# 2. MAIN PARTICLE FILTER FUNCTION
###########################################

particle_filter <- function(y, N, init_fn, transition_fn, likelihood_fn,
                            algorithm = c("SIS", "SISR", "SISAR"),
                            resample_fn = c("multinomial", "stratified", "systematic"),
                            threshold = NULL, return_particles = TRUE, ...) {

  algorithm <- match.arg(algorithm)
  resample_fn <- match.arg(resample_fn)

  if (resample_fn == "multinomial") {
    resample_func <- resample_multinomial
  } else if (resample_fn == "stratified") {
    resample_func <- resample_stratified
  } else if (resample_fn == "systematic") {
    resample_func <- resample_systematic
  }

  T_len <- length(y)
  state_est <- numeric(T_len)
  ESS_vec   <- numeric(T_len)
  if (return_particles) {
    particles_history <- matrix(NA, nrow = T_len, ncol = N)
  }

  particles <- init_fn(N, ...)
  weights <- likelihood_fn(y[1], particles, t = 1, ...)
  weights <- weights / sum(weights)

  state_est[1] <- sum(particles * weights)
  ESS_vec[1] <- 1 / sum(weights^2)
  if (return_particles) {
    particles_history[1, ] <- particles
  }

  for (t in 2:T_len) {
    particles <- transition_fn(particles, t, ...)
    likelihoods <- likelihood_fn(y[t], particles, t, ...)
    weights <- weights * likelihoods
    weights <- weights / sum(weights)

    ESS_current <- 1 / sum(weights^2)
    ESS_vec[t] <- ESS_current

    if (algorithm == "SISR") {
      particles <- resample_func(particles, weights)
      weights <- rep(1/N, N)
      ESS_vec[t] <- N  # reset ESS after resampling
    } else if (algorithm == "SISAR") {
      if (is.null(threshold)) {
        threshold <- N / 2
      }
      if (ESS_current < threshold) {
        particles <- resample_func(particles, weights)
        weights <- rep(1/N, N)
        ESS_vec[t] <- N
      }
    }
    state_est[t] <- sum(particles * weights)

    if (return_particles) {
      particles_history[t, ] <- particles
    }
  }

  result <- list(state_est = state_est, ESS = ESS_vec, algorithm = algorithm)
  if (return_particles) result$particles_history <- particles_history
  return(result)
}

###########################################
# 3. EXAMPLE SETUP: LINEAR GAUSSIAN SSM
###########################################
# SSM definitions:
#    X_t = phi * X_{t-1} + sigma_x * V_t,   V_t ~ N(0,1)
#    Y_t = X_t + sigma_y * W_t,              W_t ~ N(0,1)

init_fn_ssm <- function(N, phi, sigma_x, ...) {
  rnorm(N, mean = 0, sd = sigma_x / sqrt(1 - phi^2))
}

transition_fn_ssm <- function(particles, t, phi, sigma_x, ...) {
  phi * particles + rnorm(length(particles), mean = 0, sd = sigma_x)
}

likelihood_fn_ssm <- function(y, particles, t, sigma_y, ...) {
  dnorm(y, mean = particles, sd = sigma_y)
}

simulate_ssm <- function(T, phi, sigma_x, sigma_y) {
  x <- numeric(T)
  y <- numeric(T)
  x[1] <- rnorm(1, mean = 0, sd = sigma_x / sqrt(1 - phi^2))
  y[1] <- rnorm(1, mean = x[1], sd = sigma_y)
  for (t in 2:T) {
    x[t] <- phi * x[t - 1] + rnorm(1, mean = 0, sd = sigma_x)
    y[t] <- x[t] + rnorm(1, mean = 0, sd = sigma_y)
  }
  list(x = x, y = y)
}

###########################################
# 4. PARAMETERS
###########################################

T_val <- 50
N_particles <- 1000
phi_val <- 0.7
sigma_x_val <- 1
sigma_y_val <- 1

###########################################
# 5. ONE SIMULATION
###########################################

sim_data <- simulate_ssm(T_val, phi_val, sigma_x_val, sigma_y_val)
x_true <- sim_data$x
y_obs  <- sim_data$y

result_SIS <- particle_filter(y = y_obs, N = N_particles,
                              init_fn = init_fn_ssm,
                              transition_fn = transition_fn_ssm,
                              likelihood_fn = likelihood_fn_ssm,
                              algorithm = "SIS",
                              phi = phi_val, sigma_x = sigma_x_val, sigma_y = sigma_y_val)

result_SISR <- particle_filter(y = y_obs, N = N_particles,
                               init_fn = init_fn_ssm,
                               transition_fn = transition_fn_ssm,
                               likelihood_fn = likelihood_fn_ssm,
                               algorithm = "SISR",
                               resample_fn = "stratified",
                               phi = phi_val, sigma_x = sigma_x_val, sigma_y = sigma_y_val)

result_SISAR <- particle_filter(y = y_obs, N = N_particles,
                                init_fn = init_fn_ssm,
                                transition_fn = transition_fn_ssm,
                                likelihood_fn = likelihood_fn_ssm,
                                algorithm = "SISAR",
                                resample_fn = "stratified",
                                phi = phi_val, sigma_x = sigma_x_val, sigma_y = sigma_y_val)

rmse <- function(est, true) sqrt(mean((est - true)^2))

rmse_SIS   <- rmse(result_SIS$state_est, x_true)
rmse_SISR  <- rmse(result_SISR$state_est, x_true)
rmse_SISAR <- rmse(result_SISAR$state_est, x_true)

cat("Single simulation RMSEs:\n")
cat("  SIS   :", rmse_SIS, "\n")
cat("  SISR  :", rmse_SISR, "\n")
cat("  SISAR :", rmse_SISAR, "\n\n")

time <- 1:T_val
df_state <- data.frame(
  Time = rep(time, 4),
  State = c(x_true,
            result_SIS$state_est,
            result_SISR$state_est,
            result_SISAR$state_est),
  Method = factor(rep(c("True State", "SIS", "SISR", "SISAR"),
                      each = T_val),
                  levels = c("True State", "SIS", "SISR", "SISAR"))
)

p1 <- ggplot(df_state, aes(x = Time, y = State, color = Method, linetype = Method)) +
  geom_line(size = 1.2) +
  theme_minimal() +
  labs(title = "True State vs Particle Filter Estimates", y = "State") +
  scale_color_manual(values = c("black", "red", "blue", "green"))
print(p1)

df_ess <- data.frame(
  Time = rep(time, 3),
  ESS = c(result_SIS$ESS, result_SISR$ESS, result_SISAR$ESS),
  Method = factor(rep(c("SIS", "SISR", "SISAR"),
                      each = T_val),
                  levels = c("SIS", "SISR", "SISAR"))
)

p2 <- ggplot(df_ess, aes(x = Time, y = ESS, color = Method)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  geom_hline(yintercept = N_particles/2, linetype = "dashed", color = "gray") +
  theme_minimal()
print(p2)

###########################################
# 6. REPEATED SIMULATIONS
###########################################

n_reps <- 10000
rmse_SIS_all   <- numeric(n_reps)
rmse_SISR_all  <- numeric(n_reps)
rmse_SISAR_all <- numeric(n_reps)

for (i in 1:n_reps) {
  sim_data <- simulate_ssm(T_val, phi_val, sigma_x_val, sigma_y_val)
  x_true <- sim_data$x
  y_obs  <- sim_data$y

  result_SIS   <- particle_filter(y = y_obs, N = N_particles,
                                  init_fn = init_fn_ssm,
                                  transition_fn = transition_fn_ssm,
                                  likelihood_fn = likelihood_fn_ssm,
                                  algorithm = "SIS",
                                  phi = phi_val, sigma_x = sigma_x_val, sigma_y = sigma_y_val)

  result_SISR  <- particle_filter(y = y_obs, N = N_particles,
                                  init_fn = init_fn_ssm,
                                  transition_fn = transition_fn_ssm,
                                  likelihood_fn = likelihood_fn_ssm,
                                  algorithm = "SISR",
                                  resample_fn = "stratified",
                                  phi = phi_val, sigma_x = sigma_x_val, sigma_y = sigma_y_val)

  result_SISAR <- particle_filter(y = y_obs, N = N_particles,
                                  init_fn = init_fn_ssm,
                                  transition_fn = transition_fn_ssm,
                                  likelihood_fn = likelihood_fn_ssm,
                                  algorithm = "SISAR",
                                  resample_fn = "stratified",
                                  phi = phi_val, sigma_x = sigma_x_val, sigma_y = sigma_y_val)

  rmse_SIS_all[i]   <- rmse(result_SIS$state_est, x_true)
  rmse_SISR_all[i]  <- rmse(result_SISR$state_est, x_true)
  rmse_SISAR_all[i] <- rmse(result_SISAR$state_est, x_true)
}

rmse_SIS_mean   <- mean(rmse_SIS_all)
rmse_SIS_sd     <- sd(rmse_SIS_all)
rmse_SISR_mean  <- mean(rmse_SISR_all)
rmse_SISR_sd    <- sd(rmse_SISR_all)
rmse_SISAR_mean <- mean(rmse_SISAR_all)
rmse_SISAR_sd   <- sd(rmse_SISAR_all)

cat("RMSE over", n_reps, "replications:\n")
cat("  SIS   : Mean =", rmse_SIS_mean, "SD =", rmse_SIS_sd, "\n")
cat("  SISR  : Mean =", rmse_SISR_mean, "SD =", rmse_SISR_sd, "\n")
cat("  SISAR : Mean =", rmse_SISAR_mean, "SD =", rmse_SISAR_sd, "\n\n")

df_rmse <- data.frame(
  Method = factor(c("SIS", "SISR", "SISAR"), levels = c("SIS", "SISR", "SISAR")),
  RMSE_mean = c(rmse_SIS_mean, rmse_SISR_mean, rmse_SISAR_mean),
  RMSE_sd   = c(rmse_SIS_sd, rmse_SISR_sd, rmse_SISAR_sd)
)

p_rmse <- ggplot(df_rmse, aes(x = Method, y = RMSE_mean, fill = Method)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_errorbar(aes(ymin = RMSE_mean - RMSE_sd, ymax = RMSE_mean + RMSE_sd),
                width = 0.2) +
  theme_minimal() +
  labs(title = "RMSE of Particle Filter Methods (10000 Replications)", y = "RMSE") +
  scale_fill_manual(values = c("red", "blue", "green"))
print(p_rmse)
