library(outbreaks)
library(bayesSSM)  # My R package see https://bjarkehautop.github.io/bayesSSM/
library(ggplot2)
library(msm)
library(dplyr)

set.seed(1405)

# Load the data
outbreak_data <- outbreaks::influenza_england_1978_school

# Add a time variable and create the preprocessed data with day 0 included
outbreak_data$time <- 1:nrow(outbreak_data)
preprocessed_data <- outbreak_data[, c("time", "in_bed")]
preprocessed_data <- rbind(data.frame(time = 0, in_bed = 0), preprocessed_data)

ggplot(preprocessed_data, aes(x = time, y = in_bed)) +
  #geom_line(color = "#0072B2", size = 1.2) +
  geom_point(color = "#D55E00", size = 3) +
  labs(
    title = "Influenza Outbreak in Boarding School",
    x = "Time (Days)",
    y = "Number of Students in Bed"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "none"
  ) +
  scale_x_continuous(breaks = 0:14)
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
t_max <- nrow(preprocessed_data)-1
observations <- preprocessed_data$in_bed


###########################################
# 1. Epidemic Simulation
###########################################

epidemic_step <- function(state, t_end, lambda, gamma, n_total) {
  t <- 0
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
  states <- matrix(0, nrow = t_max + 1, ncol = 2)
  # initial state at t = 0
  states[1, ] <- c(n_total - init_infected, init_infected)
  state <- states[1, ]
  for (t in 1:t_max) {
    state <- epidemic_step(state, 1, lambda, gamma, n_total)
    states[t + 1, ] <- state
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

log_priors <- list(
  lambda = log_prior_lambda,
  gamma = log_prior_gamma
)

# --- Functions for the epidemic model ---

init_fn_epidemic <- function(particles) {
  # Return a matrix with particles rows; each row is the initial state (s, i)
  n_total <- 763
  init_infected <- 1
  init_state <- c(S = n_total - init_infected, I = init_infected)
  matrix(rep(init_state, each = particles), nrow = particles, byrow = FALSE)
}

transition_fn_epidemic <- function(particles, lambda, gamma) {
  n_total <- 763

  epidemic_step <- function(state, t_end, lambda, gamma, n_total) {
    t <- 0
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
    epidemic_step(state, t_end = 1, lambda, gamma, n_total)
  }))
  new_particles
}

log_likelihood_fn_epidemic <- function(y, particles) {
  # particles is expected to be a matrix with columns (s, i)
  dpois(y, lambda = particles[, 2], log = TRUE)

  # Negative binomial
  # dnbinom(y, size = phi_true, mu = particles[, 2], log = TRUE)
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

# Generate 100 samples and plot the prior predictive distribution along the
# observed data
prior_samples <- data.frame(
  lambda = sample_lambda(100),
  gamma = sample_gamma(100)
)

# Simulate trajectories
simulate_trajectory <- function(
    prior_samples, n_total, init_infected, t_max
) {
  prior_trajectories <- list()

  for (i in 1:nrow(prior_samples)) {
    lambda <- prior_samples$lambda[i]
    gamma <- prior_samples$gamma[i]

    # Simulate epidemic for each sample
    states <- simulate_epidemic(
      n_total, init_infected, lambda, gamma, t_max
    )

    # Add Poisson noise to the infected counts (same as observation process)
    noisy_infected <- rpois(length(states[, 2]), lambda = states[, 2])

    prior_trajectories[[i]] <- noisy_infected
  }

  return(prior_trajectories)
}

# Generate the prior predictive trajectories
prior_trajectories <- simulate_trajectory(
  prior_samples, n_total, init_infected, t_max
)

# Prepare the data for plotting
plot_data <- data.frame(
  time = rep(0:t_max, times = length(prior_trajectories) + 1),
  infected = c(observations, unlist(prior_trajectories)),
  type = rep(c(
    "Observed",
    paste0("Sampled ", 1:length(prior_trajectories))),
    each = t_max + 1
  )
)

ggplot() +
  geom_line(
    data = plot_data,
    aes(x = time, y = infected, group = type),
    color = "skyblue",
    alpha = 0.7,
    linewidth = 1.2,
    linetype = "dashed"
  ) +
  geom_line(
    data = plot_data[plot_data$type == "Observed", ],
    aes(x = time, y = infected),
    color = "darkblue",
    linewidth = 1.5
  ) +
  labs(
    x = "Time (Days)",
    y = "Observed Infected",
    title = "Prior Predictive Check"
  ) +
  scale_x_continuous(
    breaks = seq(0, 10, by = 1)
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

# Save the plot
ggsave(
  "prior_predictive_check_influenza.png",
  dpi = 300,
  width = 6.27,
  height = 4,
  units = "in"
)

###########################################
# 4. PMMH
###########################################

result_influenza <- bayesSSM::pmmh(
  y = observations,
  m = 10000,
  init_fn = init_fn_epidemic,
  transition_fn = transition_fn_epidemic,
  log_likelihood_fn = log_likelihood_fn_epidemic,
  log_priors = log_priors,
  init_params = c(lambda = 0.5, gamma = 0.5),
  burn_in = 500,
  num_chains = 4,
  param_transform = list(lambda = "log", gamma = "log"),
  verbose = TRUE,
  return_latent_state_est = TRUE,
  seed = 1405,
  num_cores = 4
)

saveRDS(result_influenza, file = "result_influenza.rds")

# Or load the result from the saved file
# result_influenza <- readRDS("result_influenza.rds")

chains <- result_influenza$theta_chain

chain_1 <- chains[[1]]
chain_2 <- chains[[2]]
chain_3 <- chains[[3]]
chain_4 <- chains[[4]]

lambda_chains <- cbind(
  chain_1$lambda, chain_2$lambda, chain_3$lambda, chain_4$lambda
)

lambda_df <- data.frame(
  value = c(
    lambda_chains
  ),
  chain = factor(rep(1:4, each = nrow(lambda_chains)))
)

ggplot(lambda_df, aes(x = value, fill = chain)) +
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

ggsave("density_plot_lambda_influenza.png",
       dpi = 300,
       width = 6.27,
       height = 4,
       units = "in"
)

gamma_chains <- cbind(
  chain_1$gamma, chain_2$gamma, chain_3$gamma, chain_4$gamma
)

gamma_df <- data.frame(
  value = c(
    gamma_chains
  ),
  chain = factor(rep(1:4, each = nrow(gamma_chains)))
)

ggplot(gamma_df, aes(x = value, fill = chain)) +
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

ggsave("density_plot_gamma_influenza.png",
       dpi = 300,
       width = 6.27,
       height = 4,
       units = "in"
)

###########################################
# 5. Inference on R_0
###########################################

combined_df_with_r0 <- data.frame(
  iteration = rep(1:nrow(lambda_chains), times = 4),
  chain = factor(rep(1:4, each = nrow(lambda_chains))),
  lambda = c(lambda_chains),
  gamma = c(gamma_chains),
  r0 = c(lambda_chains/gamma_chains)
)

m <- nrow(lambda_chains)
k <- 4

r0_matrix <- matrix(combined_df_with_r0$r0, nrow = m, ncol = k, byrow = FALSE)

r0_mean <- mean(combined_df_with_r0$r0)

r0_ci <- quantile(combined_df_with_r0$r0, probs = c(0.025, 0.975))

ess(r0_matrix)
rhat(r0_matrix)

r0_mean
r0_ci


###########################################
# 6. Posterior Predictive Check
###########################################

set.seed(1405)

combined_df <- data.frame(
  iteration = rep(1:nrow(lambda_chains), times = 4),
  chain = factor(rep(1:4, each = nrow(lambda_chains))),
  lambda = c(lambda_chains),
  gamma = c(gamma_chains)
)

num_iterations <- 100
sampled_df <- combined_df[
  sample(nrow(combined_df), num_iterations, replace = TRUE),
]

posterior_trajectories <- simulate_trajectory(
  sampled_df, n_total, init_infected, t_max
)

# Prepare the data for plotting
plot_data_posterior <- data.frame(
  time = rep(0:t_max, times = length(posterior_trajectories) + 1),
  infected = c(observations, unlist(posterior_trajectories)),
  type = rep(c(
    "Observed",
    paste0("Sampled ", 1:length(posterior_trajectories))),
    each = t_max + 1
  )
)

# Posterior Predictive Check Plot
ggplot() +
  geom_line(
    data = subset(plot_data_posterior, type != "Observed"),
    aes(x = time, y = infected, group = type),
    color = "skyblue",
    alpha = 0.7,
    linewidth = 1.2,
    linetype = "dashed"
  ) +
  geom_line(
    data = subset(plot_data_posterior, type == "Observed"),
    aes(x = time, y = infected),
    color = "darkblue",
    linewidth = 1.5
  ) +
  labs(
    x = "Time (Days)",
    y = "Observed Infected",
    title = "Posterior Predictive Check"
  ) +
  scale_x_continuous(
    breaks = seq(0, 10, by = 1)
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

ggsave("posterior_predictive_check_influenza.png",
       dpi = 300,
       width = 6.27,
       height = 4,
       units = "in"
)
