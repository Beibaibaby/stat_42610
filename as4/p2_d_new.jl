using Random
using Plots
using LinearAlgebra
using Statistics
using Measures
using ProgressMeter  # Import the ProgressMeter package

# Time steps
T = 30000  

# Different N values
N_values = [1000, 1500, 2000, 2500]

function simulate_coherence_for_p(N, p)
    alpha_c = 1 / (2 * log(N))
    n_ini = 25
    simulations_coherence = zeros(n_ini)  # Array to store coherence values for each simulation

    @showprogress 1 "Simulating for N=$N, p=$p..." for sim in 1:n_ini
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

        coherence = dot(S, xi_patterns[1]) / N
        simulations_coherence[sim] = coherence
    end
    
    return mean(simulations_coherence) 
end

# Initialize the plot
p = plot(xlabel="alpha / alpha_c", ylabel="<m^1(T)>", title="Average Coherence vs. Alpha", legend=:topright)

for N in N_values
    alpha_c = 1 / (2 * log(N))
    p_start = round(Int, 0.5 * N / (2 * log(N)))
    p_end = round(Int, 3.0 * N / (2 * log(N)))
    p_values = collect(range(p_start, p_end, length=6))  # Ensure it's a collection for progress iteration

    # Calculate average_coherences outside of the @showprogress loop
    average_coherences = [simulate_coherence_for_p(N, Int(round(p))) for p in p_values]
    alpha_rel_c = p_values / (N * alpha_c)

    plot!(p, alpha_rel_c, average_coherences, label="N=$N", marker=:circle)
end

# Adjusting the plot command for better visibility
plot!(left_margin=15mm, bottom_margin=5mm, top_margin=5mm, right_margin=5mm)

# Save the plot to a file
savefig("p2_d.png")
