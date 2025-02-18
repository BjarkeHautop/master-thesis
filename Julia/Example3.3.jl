using Distributed
# Add extra workers if not already added.
if nprocs() == 1
    addprocs(max(Sys.CPU_THREADS - 1, 1))
end

@everywhere begin
    using Random, Distributions, StatsBase, Statistics

    #############################################################
    # 1. RESAMPLING FUNCTIONS
    #############################################################
    
    # Multinomial resampling: sample indices according to the weights.
    function resample_multinomial(particles, weights)
        N = length(weights)
        indices = sample(1:N, Weights(weights), N, replace=true)
        return particles[indices]
    end

    # Stratified resampling: divides [0,1] into N strata and samples one point per stratum.
    function resample_stratified(particles, weights)
        N = length(weights)
        positions = ((rand() .+ (0:N-1)) ./ N)
        cumulative_sum = cumsum(weights)
        indices = zeros(Int, N)
        i = 1
        j = 1
        while i ≤ N
            if positions[i] < cumulative_sum[j]
                indices[i] = j
                i += 1
            else
                j += 1
            end
        end
        return particles[indices]
    end

    # Systematic resampling: similar to stratified sampling but uses one random start.
    function resample_systematic(particles, weights)
        N = length(weights)
        u0 = rand() / N
        positions = u0 .+ (0:N-1) ./ N
        cumulative_sum = cumsum(weights)
        indices = zeros(Int, N)
        j = 1
        for i in 1:N
            while positions[i] > cumulative_sum[j]
                j += 1
            end
            indices[i] = j
        end
        return particles[indices]
    end

    #############################################################
    # 2. PARTICLE FILTER FUNCTION WITH FILTERING WEIGHTS HISTORY
    #############################################################
    
    function particle_filter(y, N, init_fn, transition_fn, likelihood_fn;
                             algorithm::String="SIS",
                             resample_fn::String="multinomial",
                             threshold=nothing,
                             return_particles::Bool=true,
                             kwargs...)
        resample_func = resample_fn == "multinomial" ? resample_multinomial :
                        resample_fn == "stratified"   ? resample_stratified :
                        resample_fn == "systematic"   ? resample_systematic :
                        error("Unknown resample function")
        
        T_len = length(y)
        state_est = zeros(T_len)
        ESS_vec   = zeros(T_len)
        particles_history = return_particles ? Array{Float64}(undef, T_len, N) : nothing
        weights_history   = Array{Float64}(undef, T_len, N)
        
        # Initialization at time 1.
        particles = init_fn(N; kwargs...)
        weights   = likelihood_fn(y[1], particles, 1; kwargs...)
        weights   = weights ./ sum(weights)
        
        state_est[1] = sum(particles .* weights)
        ESS_vec[1]   = 1 / sum(weights.^2)
        if return_particles
            particles_history[1, :] = particles
        end
        weights_history[1, :] = weights

        # Filtering recursion.
        for t in 2:T_len
            particles = transition_fn(particles, t; kwargs...)
            likelihoods = likelihood_fn(y[t], particles, t; kwargs...)
            weights = weights .* likelihoods
            weights = weights ./ sum(weights)
            
            ESS_current = 1 / sum(weights.^2)
            ESS_vec[t] = ESS_current
            
            if algorithm == "SISR"
                particles = resample_func(particles, weights)
                weights = fill(1/N, N)
                ESS_vec[t] = N
            elseif algorithm == "SISAR"
                if threshold === nothing
                    threshold = N / 2
                end
                if ESS_current < threshold
                    particles = resample_func(particles, weights)
                    weights = fill(1/N, N)
                    ESS_vec[t] = N
                end
            end
            
            state_est[t] = sum(particles .* weights)
            
            if return_particles
                particles_history[t, :] = particles
            end
            weights_history[t, :] = weights
        end
        
        return (
            state_est = state_est,
            ESS = ESS_vec,
            algorithm = algorithm,
            particles_history = particles_history,
            weights_history = weights_history
        )
    end

    #############################################################
    # 3. STATE SPACE MODEL FUNCTIONS FOR SSM
    #############################################################
    
    # Initialization: X₀ ~ N(0,1)
    function init_fn_ssm(N; kwargs...)
        return randn(N)
    end

    # Transition function: Xₜ = φ*Xₜ₋₁ + sin(Xₜ₋₁) + σₓ * N(0,1)
    function transition_fn_ssm(particles, t; phi, sigma_x, kwargs...)
        return phi .* particles .+ sin.(particles) .+ sigma_x .* randn(length(particles))
    end

    # Likelihood function: Yₜ = Xₜ + σ_y * N(0,1)
    function likelihood_fn_ssm(y, particles, t; sigma_y, kwargs...)
        return pdf.(Normal.(particles, sigma_y), y)
    end

    # Simulation of SSM data.
    function simulate_ssm(T, phi, sigma_x, sigma_y)
        x = zeros(T)
        y = zeros(T)
        x[1] = randn()  # X₀ ~ N(0,1)
        y[1] = rand(Normal(x[1], sigma_y))
        for t in 2:T
            x[t] = phi * x[t - 1] + sin(x[t - 1]) + sigma_x * randn()
            y[t] = x[t] + sigma_y * randn()
        end
        return (x = x, y = y)
    end

    #############################################################
    # 4. TRANSITION DENSITY FUNCTION FOR SMOOTHING
    #############################################################
    
    function transition_density_ssm(x_current, x_prev; phi, sigma_x)
        # Density for: Xₜ = φ*x_prev + sin(x_prev) + N(0, σₓ²)
        return pdf(Normal(phi * x_prev + sin(x_prev), sigma_x), x_current)
    end

    #############################################################
    # 5. BACKWARD SMOOTHING FUNCTION (FFBSm)
    #############################################################
    
    function backward_smoothing(particles_history, weights_history; phi, sigma_x)
        T_len, N = size(particles_history)
        smooth_weights = zeros(T_len, N)
        smooth_weights[T_len, :] = weights_history[T_len, :]
        for t in (T_len-1):-1:1
            for i in 1:N
                sum_val = 0.0
                for j in 1:N
                    trans_prob = transition_density_ssm(particles_history[t+1, j],
                                                          particles_history[t, i];
                                                          phi=phi, sigma_x=sigma_x)
                    sum_val += trans_prob * smooth_weights[t+1, j]
                end
                smooth_weights[t, i] = weights_history[t, i] * sum_val
            end
            smooth_weights[t, :] /= sum(smooth_weights[t, :])
        end
        smoothed_est = [sum(particles_history[t, :] .* smooth_weights[t, :]) for t in 1:T_len]
        return (smoothed_est = smoothed_est, smooth_weights = smooth_weights)
    end

    #############################################################
    # 6. HELPER FUNCTION: RMSE
    #############################################################
    
    function to_rmse(true_vals, estimates)
        return sqrt(mean((true_vals .- estimates).^2))
    end
end  # end @everywhere block

###############################################################################
# MAIN CODE (runs on the master process)
###############################################################################

using DataFrames, Plots

# 6. PARAMETERS AND SIMULATION
T_val = 50
N_particles = 1000
phi_val = 0.7
sigma_x_val = 1.0
sigma_y_val = 1.0

# Simulate data from the SSM
sim_data = simulate_ssm(T_val, phi_val, sigma_x_val, sigma_y_val)
x_true = sim_data.x
y_obs  = sim_data.y

# Plot the true state and observations.
df = DataFrame(time=1:length(x_true), x_true=x_true, y_obs=y_obs)
plt1 = plot(df.time, df.x_true, label="True State", lw=2)
scatter!(df.time, df.x_true, label="True State", markersize=4, alpha=0.7)
scatter!(df.time, df.y_obs, label="Observed", markersize=4, alpha=0.7)
xlabel!("Time")
ylabel!("Value")
title!("True State and Observations")
display(plt1)

# 7. RUN PARTICLE FILTER AND BACKWARD SMOOTHING
result_SISAR = particle_filter(y_obs, N_particles,
                               init_fn_ssm,
                               transition_fn_ssm,
                               likelihood_fn_ssm;
                               algorithm="SISAR",
                               resample_fn="stratified",
                               phi=phi_val, sigma_x=sigma_x_val, sigma_y=sigma_y_val)

smooth_result = backward_smoothing(result_SISAR.particles_history,
                                   result_SISAR.weights_history;
                                   phi=phi_val, sigma_x=sigma_x_val)

# 8. PLOT FILTERING AND SMOOTHED ESTIMATES VS TRUE STATE
time_vec = 1:T_val
df_state = DataFrame(Time = time_vec,
                     True_State = x_true,
                     Filter_Est = result_SISAR.state_est,
                     Smooth_Est = smooth_result.smoothed_est)

plt2 = plot(df_state.Time, df_state.True_State, label="True State", lw=2, color=:black)
plot!(df_state.Time, df_state.Filter_Est, label="Filter Est", lw=2, color=:blue)
plot!(df_state.Time, df_state.Smooth_Est, label="Smooth Est", lw=2, color=:red)
xlabel!("Time")
ylabel!("State")
title!("True State, Filtering & Smoothed Estimates")
display(plt2)

rmse_filtering = sqrt(mean((result_SISAR.state_est .- x_true).^2))
rmse_smoothed  = sqrt(mean((smooth_result.smoothed_est .- x_true).^2))

println("RMSE for Filtering Estimates: ", rmse_filtering)
println("RMSE for Smoothed Estimates: ", rmse_smoothed)

###############################################################################
# 9. REPEAT MANY TIMES (PARALLEL SIMULATIONS)
###############################################################################

@everywhere function single_simulation(T_val, N_particles, phi_val, sigma_x_val, sigma_y_val)
    sim_data = simulate_ssm(T_val, phi_val, sigma_x_val, sigma_y_val)
    x_true = sim_data.x
    y_obs  = sim_data.y

    result_SISAR = particle_filter(y_obs, N_particles,
                                   init_fn_ssm,
                                   transition_fn_ssm,
                                   likelihood_fn_ssm;
                                   algorithm="SISAR",
                                   resample_fn="stratified",
                                   phi=phi_val, sigma_x=sigma_x_val, sigma_y=sigma_y_val)
    
    smooth_result = backward_smoothing(result_SISAR.particles_history,
                                       result_SISAR.weights_history;
                                       phi=phi_val, sigma_x=sigma_x_val)
    
    return to_rmse(x_true, smooth_result.smoothed_est)
end

n_simulations = 10000

# Start timing
time_before = time()

# Run parallel simulations with progress tracking
rmse_values = pmap(1:n_simulations) do i
    if i % 100 == 0
        println("Iteration: ", i)
    end
    result = single_simulation(T_val, N_particles, phi_val, sigma_x_val, sigma_y_val)
    return result
end

# Stop timing
time_after = time()

mean_rmse = mean(rmse_values)
sd_rmse   = std(rmse_values)

println("Number of simulations: ", n_simulations)
println("Mean RMSE: ", mean_rmse)
println("SD RMSE: ", sd_rmse)
println("Time taken: ", time_after - time_before, " seconds")
