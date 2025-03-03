library(ggplot2)
library(tibble)
library(extraDistr) # for half-normal functions
library(mvtnorm) # for multivariate normal sampling

source("Rhat and ESS.R")
set.seed(1405)

###########################################
# 1. RESAMPLING FUNCTIONS
###########################################

# Multinomial resampling: sample indices according to the weights.
resample_multinomial <- function(particles, weights) {
  N <- length(weights)
  indices <- sample(1:N, size = N, replace = TRUE, prob = weights)

  particles[indices]
}

# Stratified resampling: divides [0,1] into N strata and samples one point per
# stratum.
resample_stratified <- function(particles, weights) {
  N <- length(weights)
  positions <- (runif(1) + 0:(N - 1)) / N
  cumulative_sum <- cumsum(weights)
  indices <- findInterval(positions, cumulative_sum) + 1

  particles[indices]
}

# Systematic resampling: similar to stratified sampling but uses one random
# start, the rest are just shifted by a constant.
resample_systematic <- function(particles, weights) {
  N <- length(weights)
  u0 <- runif(1, 0, 1 / N)
  positions <- u0 + (0:(N - 1)) / N
  cumulative_sum <- cumsum(weights)
  indices <- findInterval(positions, cumulative_sum) + 1

  particles[indices]
}

###########################################
# 2. Particle Filter (returns a likelihood estimate)
###########################################

particle_filter_loglike <- function(
    y, N, init_fn, transition_fn, likelihood_fn,
    algorithm = c("SIS", "SISR", "SISAR"),
    resample_fn = c("multinomial", "stratified", "systematic"),
    threshold = NULL, return_particles = TRUE, ...) {
  algorithm <- match.arg(algorithm)
  resample_fn <- match.arg(resample_fn)

  if (resample_fn == "multinomial") {
    resample_func <- function(particles, weights) {
      indices <- sample(
        seq_along(weights),
        size = length(weights),
        replace = TRUE,
        prob = weights
      )
      particles[indices]
    }
  } else if (resample_fn == "stratified") {
    resample_func <- resample_stratified
  } else if (resample_fn == "systematic") {
    resample_func <- resample_systematic
  }

  T_len <- length(y)
  state_est <- numeric(T_len)
  ESS_vec <- numeric(T_len)
  loglike <- 0 # log-likelihood accumulator
  if (return_particles) {
    particles_history <- matrix(NA, nrow = T_len, ncol = N)
  }

  # Helper function: log-sum-exp trick for numerical stability
  logsumexp <- function(lw) {
    max_lw <- max(lw)
    max_lw + log(sum(exp(lw - max_lw)))
  }

  # t = 1: Initialize particles and compute log weights
  particles <- init_fn(N, ...)
  log_weights <- likelihood_fn(y[1], particles, t = 1, ...)

  # Compute incremental likelihood L_1 = (1/N)*sum(exp(log_weights))
  log_L_t <- logsumexp(log_weights) - log(N)
  loglike <- log_L_t

  # Normalize the log weights and store in probability scale
  log_normalizer <- logsumexp(log_weights)
  log_weights <- log_weights - log_normalizer
  weights <- exp(log_weights)

  state_est[1] <- sum(particles * weights)
  ESS_vec[1] <- 1 / sum(weights^2)
  if (return_particles) particles_history[1, ] <- particles

  # Loop over time steps
  for (t in 2:T_len) {
    particles <- transition_fn(particles, t, ...)
    log_likelihoods <- likelihood_fn(y[t], particles, t, ...)

    # Update log weights: add new log likelihoods to the log of the previous
    # (normalized) weights.
    log_weights_old <- log(weights)
    log_w <- log_weights_old + log_likelihoods

    # Compute the incremental log likelihood L_t
    log_L_t <- logsumexp(log_w) - log(N)
    loglike <- loglike + log_L_t

    # Normalize weights using log-sum-exp trick
    log_normalizer <- logsumexp(log_w)
    log_weights <- log_w - log_normalizer
    weights <- exp(log_weights)

    ESS_current <- 1 / sum(weights^2)
    ESS_vec[t] <- ESS_current

    if (algorithm == "SISR") {
      particles <- resample_func(particles, weights)
      weights <- rep(1 / N, N)
      ESS_vec[t] <- N
    } else if (algorithm == "SISAR") {
      if (is.null(threshold)) threshold <- N / 2
      if (ESS_current < threshold) {
        particles <- resample_func(particles, weights)
        weights <- rep(1 / N, N)
        ESS_vec[t] <- N
      }
    }
    state_est[t] <- sum(particles * weights)
    if (return_particles) particles_history[t, ] <- particles
  }

  result <- list(state_est = state_est, ESS = ESS_vec, loglike = loglike, algorithm = algorithm)
  if (return_particles) result$particles_history <- particles_history

  result
}

###########################################
# 3. Example Setup: Non-linear Gaussian SSM
###########################################
# SSM definitions:
#    X_1 ~ N(0,1)
#    X_t = phi * X_{t-1} + sin(X_{t-1}) + sigma_x * V_t,   V_t ~ N(0,1)
#    Y_t = X_t + sigma_y * W_t,              W_t ~ N(0,1)

init_fn_ssm <- function(N, ...) {
  rnorm(N, mean = 0, sd = 1)
}

transition_fn_ssm <- function(particles, t, phi, sigma_x, ...) {
  # X_t = phi*X_{t-1} + sin(X_{t-1}) + sigma_x*V_t,  V_t ~ N(0,1)
  phi * particles + sin(particles) +
    rnorm(length(particles), mean = 0, sd = sigma_x)
}

likelihood_fn_ssm <- function(y, particles, t, sigma_y, ...) {
  # Y_t ~ N(X_t, sigma_y^2)
  dnorm(y, mean = particles, sd = sigma_y, log = TRUE)
}

simulate_ssm <- function(T_val, phi, sigma_x, sigma_y) {
  x <- numeric(T_val)
  y <- numeric(T_val)
  x[1] <- rnorm(1, mean = 0, sd = 1)
  y[1] <- rnorm(1, mean = x[1], sd = sigma_y)
  for (t in 2:T) {
    x[t] <- phi * x[t - 1] + sin(x[t - 1]) + rnorm(1, mean = 0, sd = sigma_x)
    y[t] <- x[t] + rnorm(1, mean = 0, sd = sigma_y)
  }
  list(x = x, y = y)
}

###########################################
# 4. Simulate data from the SSM
###########################################

T_val <- 50
N_particles <- 100
phi_true <- 0.7
sigma_x_true <- 1
sigma_y_true <- 1

sim_data <- simulate_ssm(T_val, phi_true, sigma_x_true, sigma_y_true)
x_true <- sim_data$x
y_obs <- sim_data$y

###########################################
# 5. Define functions for the priors and proposals
###########################################

# Prior for phi ~ N(0,1)
log_prior_phi <- function(phi) {
  dnorm(phi, mean = 0, sd = 1, log = TRUE)
}

# For sigma_x and sigma_y we use a half-normal prior with sd=1.
log_prior_sigma <- function(sigma) {
  dhnorm(sigma, sigma = 1, log = TRUE)
}


###########################################
# 6. Prior predictice check
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
  sim <- simulate_ssm(T_val, params$phi, params$sigma_x, params$sigma_y)

  sim_list[[i]] <- tibble(
    time = 1:T_val,
    y = sim$y,
    sim = as.factor(i)
  )
}

sim_df <- dplyr::bind_rows(sim_list)

obs_df <- tibble(time = 1:T_val, y = y_obs)
ggplot() +
  geom_line(
    data = sim_df,
    aes(x = time, y = y, group = sim),
    color = "grey",
    alpha = 0.6
  ) +
  geom_line(
    data = obs_df,
    aes(x = time, y = y),
    color = "red",
    size = 1.2
  ) +
  labs(x = "Time", y = "y") +
  theme_minimal()


###########################################
# 7. Pilot Run Function to Choose N (for likelihood variance)
###########################################

pilot_run <- function(y, pilot_N, pilot_reps, theta, ...) {
  pilot_loglikes <- numeric(pilot_reps)
  for (i in 1:pilot_reps) {
    pf_result <- particle_filter_loglike(y, pilot_N,
      init_fn = init_fn_ssm,
      transition_fn = transition_fn_ssm,
      likelihood_fn = likelihood_fn_ssm,
      algorithm = "SISAR",
      resample_fn = "stratified",
      phi = theta[1],
      sigma_x = theta[2],
      sigma_y = theta[3],
      return_particles = FALSE
    )
    pilot_loglikes[i] <- pf_result$loglike
  }
  variance_estimate <- var(pilot_loglikes)
  target_N <- ceiling(pilot_N * variance_estimate)
  target_N <- max(target_N, 100) # Ensure a minimum number of particles
  list(
    variance_estimate = variance_estimate,
    target_N = target_N,
    pilot_loglikes = pilot_loglikes
  )
}



###########################################
# 8. Pilot chain function for PMMH
###########################################

run_pilot_chain <- function(y, pilot_M, pilot_N, init_theta, prop_sd) {
  cat("Running pilot chain to estimate posterior mean and covariance...\n")
  pilot_theta_chain <- matrix(NA, nrow = pilot_M, ncol = 3)
  pilot_loglike_chain <- numeric(pilot_M)
  current_theta <- init_theta
  pilot_theta_chain[1, ] <- current_theta
  pf_result <- particle_filter_loglike(y, pilot_N,
    init_fn = init_fn_ssm,
    transition_fn = transition_fn_ssm,
    likelihood_fn = likelihood_fn_ssm,
    algorithm = "SISAR",
    resample_fn = "stratified",
    phi = current_theta[1],
    sigma_x = current_theta[2],
    sigma_y = current_theta[3]
  )
  current_loglike <- pf_result$loglike
  pilot_loglike_chain[1] <- current_loglike

  for (m in 2:pilot_M) {
    prop_phi <- current_theta[1] + rnorm(1, mean = 0, sd = prop_sd[1])
    prop_log_sigma_x <- log(current_theta[2]) +
      rnorm(1, mean = 0, sd = prop_sd[2])
    prop_log_sigma_y <- log(current_theta[3]) +
      rnorm(1, mean = 0, sd = prop_sd[3])
    proposed_theta <- c(prop_phi, exp(prop_log_sigma_x), exp(prop_log_sigma_y))

    log_prior_current <- log_prior_phi(current_theta[1]) +
      log_prior_sigma(current_theta[2]) +
      log_prior_sigma(current_theta[3])
    log_prior_proposed <- log_prior_phi(proposed_theta[1]) +
      log_prior_sigma(proposed_theta[2]) +
      log_prior_sigma(proposed_theta[3])
    log_jacobian_current <- log(current_theta[2]) + log(current_theta[3])
    log_jacobian_proposed <- log(proposed_theta[2]) + log(proposed_theta[3])

    pf_prop <- particle_filter_loglike(y, pilot_N,
      init_fn = init_fn_ssm,
      transition_fn = transition_fn_ssm,
      likelihood_fn = likelihood_fn_ssm,
      algorithm = "SISAR",
      resample_fn = "stratified",
      phi = proposed_theta[1],
      sigma_x = proposed_theta[2],
      sigma_y = proposed_theta[3]
    )
    proposed_loglike <- pf_prop$loglike

    log_accept_ratio <- (log_prior_proposed + proposed_loglike +
      log_jacobian_proposed) -
      (log_prior_current + current_loglike + log_jacobian_current)

    if (log(runif(1)) < log_accept_ratio) {
      current_theta <- proposed_theta
      current_loglike <- proposed_loglike
    }
    pilot_theta_chain[m, ] <- current_theta
    pilot_loglike_chain[m] <- current_loglike
  }
  burn_in <- floor(pilot_M / 2)
  pilot_theta_mean <- colMeans(pilot_theta_chain[(burn_in + 1):pilot_M, ])
  pilot_theta_cov <- cov(pilot_theta_chain[(burn_in + 1):pilot_M, ])

  cat("Pilot chain posterior mean:\n")
  print(pilot_theta_mean)
  cat("Pilot chain posterior covariance:\n")
  print(pilot_theta_cov)

  # Estimate the target number of particles using the pilot run
  pilot_result <- pilot_run(y, pilot_N, pilot_reps = 10, pilot_theta_mean)
  target_N <- pilot_result$target_N
  cat("Estimated target number of particles for PMMH:", target_N, "\n")

  list(
    pilot_theta_mean = pilot_theta_mean,
    pilot_theta_cov = pilot_theta_cov,
    target_N = target_N,
    pilot_theta_chain = pilot_theta_chain,
    pilot_loglike_chain = pilot_loglike_chain
  )
}



###########################################
# 9. PMMH sampler function
###########################################

PMMH <- function(y, M, N, init_theta, proposal_cov, burn_in = 1000) {
  current_theta <- init_theta
  theta_chain <- matrix(NA, nrow = M, ncol = 3)
  loglike_chain <- numeric(M)
  state_est_chain <- vector("list", M)

  theta_chain[1, ] <- current_theta

  pf_result <- particle_filter_loglike(y, N,
    init_fn = init_fn_ssm,
    transition_fn = transition_fn_ssm,
    likelihood_fn = likelihood_fn_ssm,
    algorithm = "SISAR",
    resample_fn = "stratified",
    phi = current_theta[1],
    sigma_x = current_theta[2],
    sigma_y = current_theta[3]
  )
  current_loglike <- pf_result$loglike
  current_state_est <- pf_result$state_est
  state_est_chain[[1]] <- current_state_est

  for (m in 2:M) {
    psi_current <- c(
      current_theta[1],
      log(current_theta[2]),
      log(current_theta[3])
    )
    psi_proposed <- as.numeric(
      rmvnorm(1, mean = psi_current, sigma = proposal_cov)
    )
    proposed_theta <- c(
      psi_proposed[1],
      exp(psi_proposed[2]),
      exp(psi_proposed[3])
    )

    # Evaluate log-priors:
    log_prior_current <- log_prior_phi(current_theta[1]) +
      log_prior_sigma(current_theta[2]) +
      log_prior_sigma(current_theta[3])
    log_prior_proposed <- log_prior_phi(proposed_theta[1]) +
      log_prior_sigma(proposed_theta[2]) +
      log_prior_sigma(proposed_theta[3])

    # Jacobian corrections for the transformation sigma = exp(log_sigma)
    log_jacobian_current <- log(current_theta[2]) + log(current_theta[3])
    log_jacobian_proposed <- log(proposed_theta[2]) + log(proposed_theta[3])

    # Run particle filter for proposed parameters:
    pf_prop <- particle_filter_loglike(y, N,
      init_fn = init_fn_ssm,
      transition_fn = transition_fn_ssm,
      likelihood_fn = likelihood_fn_ssm,
      algorithm = "SISAR",
      resample_fn = "stratified",
      phi = proposed_theta[1],
      sigma_x = proposed_theta[2],
      sigma_y = proposed_theta[3]
    )
    proposed_loglike <- pf_prop$loglike

    # Compute log acceptance ratio:
    log_num <- log_prior_proposed + proposed_loglike + log_jacobian_proposed
    log_denom <- log_prior_current + current_loglike + log_jacobian_current
    log_accept_ratio <- log_num - log_denom

    if (log(runif(1)) < log_accept_ratio) {
      current_theta <- proposed_theta
      current_loglike <- proposed_loglike
      current_state_est <- pf_prop$state_est
    }
    theta_chain[m, ] <- current_theta
    loglike_chain[m] <- current_loglike
    state_est_chain[[m]] <- current_state_est
  }

  # Remove burn-in samples:
  theta_chain_burned <- theta_chain[(burn_in + 1):M, , drop = FALSE]
  loglike_chain_burned <- loglike_chain[(burn_in + 1):M]

  # Compute posterior estimates:
  phi_estimate <- mean(theta_chain_burned[, 1])
  sigma_x_estimate <- mean(theta_chain_burned[, 2])
  sigma_y_estimate <- mean(theta_chain_burned[, 3])

  latent_state_matrix <- do.call(rbind, state_est_chain[(burn_in + 1):M])
  latent_state_estimate <- colMeans(latent_state_matrix)

  list(
    phi_estimate = phi_estimate,
    sigma_x_estimate = sigma_x_estimate,
    sigma_y_estimate = sigma_y_estimate,
    theta_chain = theta_chain_burned,
    loglike_chain = loglike_chain_burned,
    latent_state_chain = state_est_chain,
    latent_state_estimate = latent_state_estimate
  )
}


###########################################
# 10. Run PMMH on the simulated data
###########################################

# Initialize parameters from the prior:
sample_phi <- function() rnorm(1, mean = 0, sd = 1)
sample_sigma <- function() rht(1, nu = 3, sigma = 2.5) # half-t sample

init_theta <- c(sample_phi(), sample_sigma(), sample_sigma())
# prop_sd is needed for the pilot chain
prop_sd <- c(0.1, 0.1, 0.1)

# First, estimate the proposal covariance using a pilot run.
pilot_M <- 2000
pilot_N <- 100

pilot_chain <- run_pilot_chain(y_obs, pilot_M, pilot_N, init_theta, prop_sd)

init_theta <- pilot_chain$pilot_theta_mean
proposal_cov_est <- pilot_chain$pilot_theta_cov
N_particles <- pilot_chain$target_N

M_iter <- 15000
burn_in <- 2000

# Now run the main PMMH using the multivariate normal proposal.
pmmh_out <- PMMH(y_obs, M_iter, N_particles, init_theta, proposal_cov_est,
  burn_in = burn_in
)

###########################################
# 11. Posterior predictive check
###########################################

# Number of posterior predictive simulations to generate:
n_ppc <- 4

posterior_samples <- pmmh_out$theta_chain

ppc_list <- vector("list", n_ppc)

# Randomly sample indices from the posterior chain
sample_indices <- sample(seq_len(posterior_samples), n_ppc)

# Loop over the selected posterior samples and simulate new datasets:
for (i in 1:n_ppc) {
  # Extract a set of parameters from the posterior sample
  params <- posterior_samples[sample_indices[i], ]
  # Use simulate_ssm to generate replicated data using these parameters
  sim <- simulate_ssm(T_val, params[1], params[2], params[3])

  ppc_list[[i]] <- tibble(
    time = 1:T_val,
    y = sim$y,
    sim = as.factor(i)
  )
}

# Combine all simulated datasets into one data frame
ppc_df <- dplyr::bind_rows(ppc_list)

# Plot the posterior predictive simulations along with the observed data:
obs_df <- tibble(time = 1:T_val, y = y_obs)
ggplot() +
  geom_line(
    data = ppc_df, aes(x = time, y = y, group = sim),
    color = "grey", alpha = 0.6
  ) +
  geom_line(
    data = obs_df, aes(x = time, y = y),
    color = "red", size = 1.2
  ) +
  labs(x = "Time", y = "Observed y") +
  theme_minimal()


###########################################
# 12. Inference
###########################################

# Posterior estimates for phi, sigma_x, and sigma_y after burn-in
cat("Posterior estimates after burn-in:\n")
cat("phi estimate: ", pmmh_out$phi_estimate, "\n")
cat("sigma_x estimate: ", pmmh_out$sigma_x_estimate, "\n")
cat("sigma_y estimate: ", pmmh_out$sigma_y_estimate, "\n")

# Plot traces of the parameter chains:
par(mfrow = c(3, 1))
plot(
  pmmh_out$theta_chain[, 1],
  type = "l",
  main = expression(phi),
  ylab = "phi",
  xlab = "Iteration"
)
plot(
  pmmh_out$theta_chain[, 2],
  type = "l",
  main = expression(sigma[x]),
  ylab = expression(sigma[x]),
  xlab = "Iteration"
)
plot(
  pmmh_out$theta_chain[, 3],
  type = "l",
  main = expression(sigma[y]),
  ylab = expression(sigma[y]),
  xlab = "Iteration"
)


###########################################
# 13. Run 4 chains of PMMH
###########################################

M_iter <- 15000
burn_in <- 2000

# Initialize empty lists to store the chains
chains <- list()

# Run the PMMH 4 times and store each result in the chains list
for (i in 1:4) {
  print(paste("Running chain", i, "..."))
  init_theta <- c(sample_phi(), sample_sigma(), sample_sigma())
  chains[[paste0("chain_", i)]] <- PMMH(
    y_obs,
    M_iter,
    N_particles,
    init_theta,
    proposal_cov_est,
    burn_in = burn_in
  )
}

# Access individual chains
chain_1 <- chains$chain_1
chain_2 <- chains$chain_2
chain_3 <- chains$chain_3
chain_4 <- chains$chain_4

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


###########################################
# 14. Density plot of the 4 chains
###########################################

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


###########################################
# 15. ESS and split_rhat for the 4 chains
###########################################

# Run the split_rhat function for each parameter
phi_rhat <- split_rhat(phi_chains)
sigma_x_rhat <- split_rhat(sigma_x_chains)
sigma_y_rhat <- split_rhat(sigma_y_chains)

# Print the Rhat values
print(paste("Rhat for phi:", phi_rhat))
print(paste("Rhat for sigma_x:", sigma_x_rhat))
print(paste("Rhat for sigma_y:", sigma_y_rhat))


# Run the effectiveSize function for each parameter
phi_ess <- effectiveSize(phi_chains)
sigma_x_ess <- effectiveSize(sigma_x_chains)
sigma_y_ess <- effectiveSize(sigma_y_chains)


# Print the ESS values
print(paste("ESS for phi:", phi_ess))
print(paste("ESS for sigma_x:", sigma_x_ess))
print(paste("ESS for sigma_y:", sigma_y_ess))

###########################################
# 16.Inference
############################################

# Combine chains to estimate parameter means
phi_estimates <- mean(colMeans(phi_chains))
sigma_x_estimates <- mean(colMeans(sigma_x_chains))
sigma_y_estimates <- mean(colMeans(sigma_y_chains))

# Credible interval
phi_ci <- quantile(phi_chains, c(0.025, 0.975))
sigma_x_ci <- quantile(sigma_x_chains, c(0.025, 0.975))
sigma_y_ci <- quantile(sigma_y_chains, c(0.025, 0.975))

# Print
print(paste("Mean phi estimate:", round(phi_estimates, 2)))
print(paste("Mean sigma_x estimate:", round(sigma_x_estimates, 2)))
print(paste("Mean sigma_y estimate:", round(sigma_y_estimates, 2)))

print(paste("95% CI for phi:", round(phi_ci, 2)))
print(paste("95% CI for sigma_x:", round(sigma_x_ci, 2)))
print(paste("95% CI for sigma_y:", round(sigma_y_ci, 2)))


# RMSE of the latent state
chains_latent <- cbind(
  chain_1$latent_state_estimate,
  chain_2$latent_state_estimate,
  chain_3$latent_state_estimate,
  chain_4$latent_state_estimate
)

latent_mean <- apply(chains_latent, 1, mean)

rmse <- function(est, true) sqrt(mean((est - true)^2))
cat("RMSE of latent state estimate:", rmse(latent_mean, x_true), "\n")
