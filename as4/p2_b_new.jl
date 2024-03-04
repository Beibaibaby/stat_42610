using Random
using Plots
using LinearAlgebra
using Statistics
using Measures
using ProgressMeter  # Import the ProgressMeter package

# Parameters
T = 30000  # Time steps
N = 1000  # Number of neurons
alpha_c = 1 / (2 * log(N))
p_start = round(Int, 0.5 * N / (2 * log(N)))  # Starting value of p, rounded to nearest integer
p_end = round(Int, 3.0 * N / (2 * log(N)))  # Ending value of p, rounded to nearest integer
p_values = collect(range(p_start, p_end, length=6))  # Sample p values linearly in the range and ensure it's a collection for progress iteration

function simulate_coherence_for_p(p)
    coherences = zeros(30)  # Store coherence for each run
    
    @showprogress 1 "Simulating for p=$p..." for simulation in 1:30
        # Generate p patterns
        xi_patterns = [rand([-1, 1], N) for _ in 1:p]

        # Compute the weight matrix J
        J = zeros(N, N)
        for i in 1:N
            for j in 1:N
                if i != j
                    J[i, j] = sum(xi_patterns[μ][i] * xi_patterns[μ][j] for μ in 1:p) / N
                end
            end
        end

        # Initialize the system with the first pattern
        S = copy(xi_patterns[1])

        # Update the system for T steps
        for _ in 1:T
            i = rand(1:N)
            S[i] = sign(dot(J[i, :], S)) == 0 ? 1 : sign(dot(J[i, :], S))
        end

        # Compute coherence relative to the first pattern
        coherence = dot(S, xi_patterns[1]) / N
        coherences[simulation] = coherence
    end
    
    return mean(coherences)  # Return the average coherence over the runs
end

# Prepare for progress display
p_len = length(p_values)
progress = Progress(p_len, desc="Computing coherences")

average_coherences = zeros(p_len)
for (i, p) in enumerate(p_values)
    average_coherences[i] = simulate_coherence_for_p(Int(round(p)))
    next!(progress)
end

# Plotting
alpha_rel_c = p_values ./ (N * alpha_c)  # Compute alpha relative to alpha_c
xticks = minimum(alpha_rel_c):0.1:maximum(alpha_rel_c)  # Define x-ticks

plot(alpha_rel_c, average_coherences, xlabel="alpha / alpha_c", ylabel="<m^1(T)>", title="Average Coherence vs. Alpha", legend=false, left_margin=15mm, bottom_margin=5mm, top_margin=5mm, right_margin=5mm)

# Save the plot to a file
savefig("p2_b.png")
