using Random
using Plots
using LinearAlgebra
using Measures

# Parameters
N = 1000  # Number of neurons
num_patterns = 300 # This can be changed as needed
updates = 50000  # Number of updates
alpha = 1 / (2 * log(N))
println("alpha = ", alpha)
CL = alpha * N
println("Greatest Storage = ", CL)

# Generate num_patterns random patterns
xi_patterns = [rand([-1, 1], N) for _ in 1:num_patterns]

# Compute weights J_ij, ensuring no self-connections
J = zeros(N, N)
for i in 1:N
    for j in 1:N
        if i != j
            J[i, j] = sum(xi_patterns[μ][i] * xi_patterns[μ][j] for μ in 1:num_patterns) / N
        end
    end
end

# Initial state S_i(0) set to the first pattern
S = copy(xi_patterns[1])

# Function to update a randomly chosen neuron
function update_neuron!(S, J)
    i = rand(1:N)
    h = dot(J[i, :], S)
    S[i] = sign(h) == 0 ? 1 : sign(h)
end

# Function to compute coherence m^mu(t)
coherence(S, pattern) = dot(S, pattern) / N

# Track the evolution of coherence for each pattern
m_evolution = zeros(num_patterns, updates)

# Simulation loop
for t in 1:updates
    update_neuron!(S, J)
    for μ in 1:num_patterns
        m_evolution[μ, t] = coherence(S, xi_patterns[μ])
    end
end


# Dynamic label creation for patterns
plot_obj = plot(xlabel="Time step", ylabel="Coherence m^mu(t)",
                title="Evolution of Coherence, num_patterns = " * string(num_patterns), left_margin=10mm, bottom_margin=5mm,
                top_margin=5mm, right_margin=5mm, legend=:topleft)

# Add each pattern to the plot. Label only "Pattern 1"
for μ in 1:num_patterns
    if μ == 1
        plot!(plot_obj, 1:updates, m_evolution[μ, :], label="Pattern 1", linewidth=2)
    else
        plot!(plot_obj, 1:updates, m_evolution[μ, :], label=false, linewidth=2)
    end
end

# Save the plot to a file
savefig(plot_obj, "p2_p=$(num_patterns).png")
