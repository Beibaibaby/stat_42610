using Distributions
using Plots
using Measures

function simulate_process(N_max, trials)
    λ_xi = 100
    λ_xc = 5
    variances = zeros(N_max)

    for N in 1:N_max
        Y_vars = zeros(trials)
        
        for trial in 1:trials
            xi = rand(Poisson(λ_xi), N)  # Generate x_i's
            xc = rand(Poisson(λ_xc))  # Generate x_c
            yi = xi .+ xc  # Generate y_i's by adding x_c to each x_i
            Y = sum(yi)/N  # Sum of y_i's for this trial
            Y_vars[trial] = Y  # Store Y for variance calculation
        end
        
        variances[N] = var(Y_vars)  # Compute variance of Y for this N
    end

    return variances
end

N_max = 200
trials = 500
variances = simulate_process(N_max, trials)

# Theoretical prediction for Var(Y) when N=1
theoretical_variances = [(105*(1/N+1/21*(N-1)/(N))) for N in 1:N_max]
# Plotting
p=plot(1:N_max, variances, label="Numerical", xlabel="N", ylabel="Var(Y)", title="Numerical vs. Theoretical Var(Y)",left_margin=20mm,bottom_margin=15mm)
plot!(1:N_max, theoretical_variances, label="Theoretical", linestyle=:dash)

# Define the filename and path for saving the plot
filename = "VarY_vs_Theoretical_Prediction.png"

# Save the plot
savefig(p, filename)

# This saves the plot in the current directory. You can specify a different path if needed.
