function resample_multinomial(particles, weights::AbstractVector)
    n = size(particles, 1)
    @assert n == length(weights) "Number of particles must match the length of weights"

    norm_weights = weights / sum(weights)
    cumw = cumsum(norm_weights)
    inds = [searchsortedfirst(cumw, rand()) for _ = 1:n]

    if isa(particles, AbstractMatrix)
        return particles[inds, :]
    elseif isa(particles, AbstractVector)
        return particles[inds]
    else
        throw(ArgumentError("particles must be a vector or a matrix"))
    end
end

function resample_stratified(particles, weights::AbstractVector)
    n = size(particles, 1)
    @assert n == length(weights) "Number of particles must match the length of weights"
    u0 = rand()
    positions = (u0 .+ (0:(n-1))) ./ n

    norm_weights = weights / sum(weights)
    cumw = cumsum(norm_weights)

    inds = Vector{Int}(undef, n)
    j = 1
    for i = 1:n
        while positions[i] > cumw[j]
            j += 1
        end
        inds[i] = j
    end

    if isa(particles, AbstractMatrix)
        return particles[inds, :]
    elseif isa(particles, AbstractVector)
        return particles[inds]
    else
        throw(ArgumentError("particles must be a vector or a matrix"))
    end
end

function resample_systematic(particles::AbstractMatrix, weights::AbstractVector)
    n, _ = size(particles)
    @assert n == length(weights) "Number of particles must match length of weights"
    start = rand() / n
    positions = start .+ (0:(n-1)) ./ n
    cumw = cumsum(weights)
    inds = Vector{Int}(undef, n)
    j = 1
    for i = 1:n
        while positions[i] > cumw[j]
            j += 1
        end
        inds[i] = j
    end

    if isa(particles, AbstractMatrix)
        return particles[inds, :]
    elseif isa(particles, AbstractVector)
        return particles[inds]
    else
        throw(ArgumentError("particles must be a vector or a matrix"))
    end
end


function particle_filter(y, N, init_fn, transition_fn, likelihood_fn;
    algorithm::String = "SIS",
    resample_fn::String = "multinomial",
    threshold = nothing,
    return_particles::Bool = true,
    kwargs...)
    # Select resampling function
    resample_func = resample_fn == "multinomial" ? resample_multinomial :
    resample_fn == "stratified"   ? resample_stratified :
    resample_fn == "systematic"   ? resample_systematic :
    error("Unknown resample function: $resample_fn")

    T = length(y)
    # We have T+1 time points: t=0 (prior) through t=T
    state_est = zeros(T+1)
    ESS_vec   = zeros(T+1)
    particles_history = return_particles ? Array{Float64}(undef, T+1, N) : nothing
    weights_history   = Array{Float64}(undef, T+1, N)

    #---- Initialization at t = 0 ----#
    particles = init_fn(N; kwargs...)
    # initial weights uniform
    weights = fill(1/N, N)

    # store initial estimate and ESS
    state_est[1] = sum(particles .* weights)
    ESS_vec[1]   = 1 / sum(weights .^ 2)
    if return_particles
        particles_history[1, :] = particles
    end
    weights_history[1, :] = weights

    # Set default threshold for SISAR if not provided
    if algorithm == "SISAR" && threshold === nothing
        threshold = N/2
    end

    #---- Filtering recursion for t = 1:T ----#
    for t in 1:T
        # propagate particles
        particles = transition_fn(particles, t; kwargs...)

        # compute and normalize weights
        lik = likelihood_fn(y[t], particles, t; kwargs...)
        weights .= weights .* lik
        weights ./= sum(weights)

        # compute ESS
        ess = 1 / sum(weights .^ 2)

        # resampling if needed
        if algorithm == "SISR"
            particles = resample_func(particles, weights)
            weights .= 1/N
            ess = N
        elseif algorithm == "SISAR" && ess < threshold
            particles = resample_func(particles, weights)
            weights .= 1/N
            ess = N
        end

        # time index in storage is t+1
        idx = t + 1
        state_est[idx] = sum(particles .* weights)
        ESS_vec[idx]   = ess
        if return_particles
            particles_history[idx, :] = particles
        end
        weights_history[idx, :] = weights
    end

    return (
        state_est = state_est,
        ESS       = ESS_vec,
        algorithm = algorithm,
        particles_history = particles_history,
        weights_history   = weights_history
    )
end


# Example usage: Non-linear Gaussian SSM
using Distributions, Random, Plots, Statistics

# Model functions
init_fn(N; phi, sigma_x, sigma_y ) = rand(Normal(0,1), N)
transition_fn(particles, t; phi, sigma_x, sigma_y) = phi .* particles .+ sin.(particles) .+ rand(Normal(0, sigma_x), length(particles))
likelihood_fn(y, particles, t; phi, sigma_x, sigma_y) = pdf.(Normal.(particles, sigma_y), y)

# Simulate SSM
function simulate_ssm(t_val, phi, sigma_x, sigma_y)
    x = Vector{Float64}(undef, t_val+1)
    y = Vector{Float64}(undef, t_val)
    x[1] = rand(Normal(0, sigma_x))
    for t in 2:t_val+1
        x[t] = phi * x[t-1] + sin(x[t-1]) + rand(Normal(0, sigma_x))
        if t > 1
            y[t-1] = x[t] + rand(Normal(0, sigma_y))
        end
    end
    return x, y
end

# Parameters
t_val = 50
n_particles = 1000
phi = 0.7
sigma_x = 1.0
sigma_y = 1.0

# Run simulation
Random.seed!(1405)
x_true, y_obs = simulate_ssm(t_val, phi, sigma_x, sigma_y)

# Run particle filters
res_sis   = particle_filter(y_obs,   n_particles, init_fn, transition_fn, likelihood_fn;
                             algorithm="SIS",   resample_fn="stratified",
                             phi=phi, sigma_x=sigma_x, sigma_y=sigma_y)
res_sisr  = particle_filter(y_obs,   n_particles, init_fn, transition_fn, likelihood_fn;
                             algorithm="SISR",  resample_fn="stratified",
                             phi=phi, sigma_x=sigma_x, sigma_y=sigma_y)
res_sisar = particle_filter(y_obs,   n_particles, init_fn, transition_fn, likelihood_fn;
                             algorithm="SISAR", resample_fn="stratified",
                             phi=phi, sigma_x=sigma_x, sigma_y=sigma_y)

# Compute RMSE
rmse(est, true_val) = sqrt(mean((est .- true_val).^2))
println("RMSE SIS:   ", rmse(res_sis.state_est, x_true))
println("RMSE SISR:  ", rmse(res_sisr.state_est, x_true))
println("RMSE SISAR: ", rmse(res_sisar.state_est, x_true))

# Plot results
time = 1:t_val
p = plot(time, x_true[2:end], label="X_t", lw=2)
scatter!(p, time, y_obs, label="Y_t", ms=3, alpha=0.5)
title!(p, "Latent and Observed States Over Time")
xlabel!(p, "Time")
ylabel!(p, "State Value")



# Transistion density
function transition_density_ssm(x_curr::Real, x_prev::Real, phi::Real, sigma_x::Real)
    mean = phi * x_prev + sin(x_prev)
    return pdf(Normal(mean, sigma_x), x_curr)
end


# FFBSm

function backward_smoothing(
    particles_history::AbstractMatrix{<:Real},
    weights_history::AbstractMatrix{<:Real},
    phi::Real,
    sigma_x::Real
)
    T1, N = size(particles_history)     # T1 = T+1
    smooth_weights = zeros(eltype(weights_history), T1, N)

    # Initialize at final time T
    smooth_weights[T1, :] .= weights_history[T1, :]

    # Backward recursion
    for t in (T1-1):-1:1
        # Preallocate an array for the sums over future particles
        sums = zeros(eltype(weights_history), N)
        for i in 1:N
            # For each particle at time t, compute sum_j p(x_{t+1}^j | x_t^i) * w_{t+1}^j
            for j in 1:N
                probs = transition_density_ssm(
                    particles_history[t+1, j],
                    particles_history[t,   i],
                    phi, sigma_x
                )
                sums[i] += probs * smooth_weights[t+1, j]
            end
            # Weight by the filtering weight at time t
            smooth_weights[t, i] = weights_history[t, i] * sums[i]
        end
        # Normalize smoothing weights at time t
        smooth_weights[t, :] ./= sum(smooth_weights[t, :])
    end

    # Compute smoothed state estimates
    smoothed_est = similar(sum(weights_history; dims=2))
    for t in 1:T1
        smoothed_est[t] = sum(particles_history[t, :] .* smooth_weights[t, :])
    end

    return (smoothed_est = smoothed_est, smooth_weights = smooth_weights)
end


# --- after your res_sisar = particle_filter(...) call ---

# 1. Extract the stored particles & weights
particles_hist = res_sisar.particles_history   # (T+1) × N
weights_hist   = res_sisar.weights_history     # (T+1) × N

# 2. Run backward smoothing
sm_res = backward_smoothing(particles_hist, weights_hist, phi, sigma_x)
smoothed_est = sm_res.smoothed_est             # length T+1

# 3. Plot true, filtered, and smoothed trajectories
time = 0:t_val   # recall state_est is indexed 1→t=0, 2→t=1, …, T+1→t=T
filtered_est = res_sisar.state_est
println("RMSE smoothing: ", rmse(smoothed_est, x_true))

using Plots
p = plot(time, x_true,               label="True state",       lw=2)
plot!(p, time, filtered_est,         label="Filtered mean",    lw=2, ls=:dash)
plot!(p, time, smoothed_est,         label="Smoothed mean",     lw=2, ls=:dot)
xlabel!(p, "Time")
ylabel!(p, "State")
title!(p, "Filtering vs. Smoothing")



# Repeat 10,000 times
times = 10000
smooth_res = Array{Float64}(undef, times)
for i in 1:times
    if i % 500 == 0
        println("Iteration ", i)
    end
    Random.seed!(1405+i)
    x_true, y_obs = simulate_ssm(t_val, phi, sigma_x, sigma_y)

    # Run particle filters
    res_sisar = particle_filter(y_obs,   n_particles, init_fn, transition_fn, likelihood_fn;
                                algorithm="SISAR", resample_fn="stratified",
                                phi=phi, sigma_x=sigma_x, sigma_y=sigma_y)

    # 1. Extract the stored particles & weights
    particles_hist = res_sisar.particles_history   # (T+1) × N
    weights_hist   = res_sisar.weights_history     # (T+1) × N

    # 2. Run backward smoothing
    sm_res = backward_smoothing(particles_hist, weights_hist, phi, sigma_x)
    smooth_res[i] = rmse(sm_res.smoothed_est, x_true)
end

using DelimitedFiles

writedlm("smooth_res.csv", smooth_res, ',')
smooth_res = readdlm("smooth_res.csv", ',')[:]

mean_smooth_res = mean(smooth_res)
sd_smooth_res = std(smooth_res)

println("mean RMSE smoothing: ", mean_smooth_res)
println("sd RMSE smoothing: ", sd_smooth_res)
