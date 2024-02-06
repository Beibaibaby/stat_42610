using Distributions
using Plots
using Statistics
using Measures
using QuadGK

# Analytical solution for ρ(τ)
function rho_analytical(τ, λ, τ_abs, τ_rel)
    if τ <= τ_abs
        return 0.0
    else
        H_τ = λ * (1 - exp(-(τ - τ_abs) / τ_rel))
        S_τ = exp(-λ * (τ - τ_abs - τ_rel + τ_rel * exp(-(τ - τ_abs) / τ_rel)))
        return H_τ * S_τ
    end
end

function simulate_spike_times(λ, τ_abs, τ_rel, T, dt)
    t = 0.0
    spike_times = []
    while t < T
        # Calculate hazard rate at current time, considering λ in Hz and converting dt to seconds for this calculation
        rate = λ * (1 - exp(-(t - (isempty(spike_times) ? 0 : spike_times[end])-τ_abs) / τ_rel))
        p_spike = 1 - exp(-rate * dt / 1000)  # Convert dt to seconds for this calculation because λ is in Hz
        if rand() < p_spike
            push!(spike_times, t)
            t += τ_abs  # Ensure τ_abs is in ms, consistent with t and dt
            #to the next iteration
            continue
        end
        t += dt
    end
    return spike_times
end




# Corrected function for numerical estimation and plotting
function estimate_rho(spike_times, τ_abs, τ_rel, λ)
    isi = diff(spike_times)
    # Assuming continuous intervals for plotting
    τ_vals = minimum(isi):0.001:maximum(isi)
    ρ_analytical_vals = [rho_analytical(τ, λ/1000, τ_abs, τ_rel) for τ in τ_vals]

    # Plot numerical estimation as histogram
    histogram(isi, bins=100, normed=true, alpha=0.5, label="Numerical τ_abs=$τ_abs, τ_rel=$τ_rel",left_margin=20mm,bottom_margin=20mm,size=(800,600))

    # Overlay analytical solution
    plot!(τ_vals, ρ_analytical_vals, label="Analytical τ_abs=$τ_abs, τ_rel=$τ_rel", lw=2)
end


# Function to compute CV from spike times
function compute_cv(spike_times)
    isi = diff(spike_times)
    return std(isi) / mean(isi)
end

function rho(tau, lambda, tau_rel, tau_abs)
    return lambda * (1 - exp(-(tau - tau_abs) / tau_rel)) * exp(-lambda * (tau - tau_abs + tau_rel * exp(-(tau - tau_abs) / tau_rel) - tau_rel))
end


# Function to compute the analytic CV
function compute_analytic_cv(lambda, tau_abs, tau_rel)
    mean_int, _ = quadgk(tau -> tau * rho(tau, lambda, tau_rel, tau_abs), tau_abs, Inf)
    second_moment_int, _ = quadgk(tau -> tau^2 * rho(tau, lambda, tau_rel, tau_abs), tau_abs, Inf)
    cv = sqrt(second_moment_int - mean_int^2) / mean_int
    return cv
end

# Initialize parameters
λ = 30.0 / 1000  # Hz to kHz for consistency with ms
τ_abs = 3.0  # ms, fixed absolute refractory period
τ_rel_values = 0:2:20  # ms, range of relative refractory periods
num_simulations = 1  # Number of simulations per τ_rel value
ind = 0

# Arrays to store results
τ_rel_array = collect(τ_rel_values)
cv_means = []
cv_analytic_array = []


for τ_rel in τ_rel_values
    global ind
    println(ind)
    ind += 1
    cvs = Float64[]  # To store CVs for each simulation
    for sim in 1:num_simulations
        spike_times = simulate_spike_times(λ * 1000, τ_abs, τ_rel, 1e6, 0.1)  # Adjust λ back to Hz
        cv = compute_cv(spike_times)
        push!(cvs, cv)
    end
    push!(cv_means, mean(cvs))  # Compute and store the average CV
    
    # Compute and store the analytic CV
    analytic_cv = compute_analytic_cv(λ, τ_abs, τ_rel)
    push!(cv_analytic_array, analytic_cv)
end

# Plotting both numerical and analytic CV as functions of τ_rel
p = plot(τ_rel_array, cv_means, xlabel="τ_rel (ms)", ylabel="Coefficient of Variation (CV)", 
         title="CV as a Function of τ_rel", label="Numerical CV", marker=:circle, legend=:topright, size=(800,600),left_margin=20mm,bottom_margin=20mm)

plot!(p, τ_rel_array, cv_analytic_array, label="Analytic CV", marker=:square, line=:dash, color=:red)

# Save the plot
savefig(p, "p3_z.png")
