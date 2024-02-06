using Distributions
using Plots
using Statistics
using Measures

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
        rate = λ * (1 - exp(-(t - (isempty(spike_times) ? 0 : spike_times[end])) / τ_rel))
        p_spike = 1 - exp(-rate * dt / 1000)  # Convert dt to seconds for this calculation if λ is in Hz
        if rand() < p_spike
            push!(spike_times, t)
            t += τ_abs  # Ensure τ_abs is in ms, consistent with t and dt
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

# Initialize parameters
λ = 30.0  # Hz
τ_rel = 3.0  # ms, fixed absolute refractory period
τ_abs_values = 0:1:20  # ms, range of relative refractory periods
num_simulations = 10  # Number of simulations per τ_rel value

# Arrays to store results
τ_abs_array = []
cv_means = []

for τ_abs in τ_abs_values
    cvs = Float64[]  # To store CVs for each simulation
    for sim in 1:num_simulations
        spike_times = simulate_spike_times(λ, τ_abs, τ_rel, 1e5,0.1)
        cv = compute_cv(spike_times)
        push!(cvs, cv)
    end
    push!(τ_abs_array, τ_abs)
    push!(cv_means, mean(cvs))  # Compute and store the average CV
end

# Plotting average CV as a function of τ_rel
p = plot(τ_abs_array, cv_means, xlabel="τ_abs (ms)", ylabel="Coefficient of Variation (CV)",
         title="Average CV as a Function of τ_abs", legend=false, marker=:circle,left_margin=20mm,bottom_margin=20mm,size=(800,600),left_margin)

# Save the plot
savefig(p, "p3_e.png")


