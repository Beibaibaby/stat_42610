using Random
using Plots
using LinearAlgebra
using Statistics
using Measures

# Time steps
T = 30000  

# Different N values
N_values = [1000,2500,5000,10000]

function simulate_coherence_for_p(N, p)
    alpha_c = 1 / (2 * log(N))
    simulations_coherence = zeros(5) # Array to store coherence values for each simulation

    for sim in 1:5 # Run each simulation 5 times
        xi_patterns = [rand([-1, 1], N) for _ in 1:p]
        J = zeros(N, N)

        for i in 1:N
            for j in 1:N
                if i != j
                    J[i, j] = sum(xi_patterns[μ][i] * xi_patterns[μ][j] for μ in 1:p) / N
                end
            end
        end

        S = copy(xi_patterns[1])

        for _ in 1:T
            i = rand(1:N)
            S[i] = sign(dot(J[i, :], S)) == 0 ? 1 : sign(dot(J[i, :], S))
        end

        coherences = [dot(S, xi_patterns[μ]) / N for μ in 1:p]
        average_coherence = mean(coherences)
        
        simulations_coherence[sim] = average_coherence
    end
    
    return mean(simulations_coherence) # Return the mean of the 5 simulation averages
end

# Initialize the plot
p = plot(xlabel="p", ylabel="<m^1(T)>", title="Average Coherence vs. Alpha", legend=:topright)

for N in N_values
    alpha_c = 1 / (2 * log(N))
    p_start = round(Int, 0.5 * 1000 / (2 * log(1000)))
    p_end = round(Int, 2 * 1000 / (2 * log(1000)))
    p_values = range(p_start, p_end, length=10)

    # Compute the average coherence for each p, averaged over 5 simulations
    average_coherences = [mean([simulate_coherence_for_p(N, Int(round(p))) for _ in 1:5]) for p in p_values]
    alpha_rel_c = p_values / (N * alpha_c)

    plot!(p, p_values, average_coherences, label="N=$N", marker=:circle)
end

# Adjusting the plot command for better visibility
plot!(left_margin=15mm, bottom_margin=5mm, top_margin=5mm, right_margin=5mm)

# Save the plot to a file
savefig("p2_e.png")