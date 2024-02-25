using Random
using Plots
using LinearAlgebra
using Measures

# Parameters
N = 1000  # Number of neurons
num_patterns = 2  # Number of patterns, renamed from p to avoid conflict with plot object
updates = 30000  # Number of updates
alpha = 1/(2*log(N))  
println("alpha = ", alpha)
CL=alpha*N
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
    i = rand(1:N)  # Randomly select a neuron
    h = dot(J[i, :], S)  # Compute the weighted input
    S[i] = sign(h) == 0 ? 1 : sign(h)  # Update state
end

# Function to compute coherence m^mu(t)
coherence(S, pattern) = dot(S, pattern) / N

# Track the evolution of coherence for each pattern
m_evolution = zeros(num_patterns, updates)

# Simulation loop
for t in 1:updates
    update_neuron!(S, J)  # Update a neuron
    for μ in 1:num_patterns  # Calculate coherence for each pattern
        m_evolution[μ, t] = coherence(S, xi_patterns[μ])
    end
end

# Plotting with different colors and corrected title
plot_obj = plot(1:updates, m_evolution', label=["Pattern 1" "Pattern 2"], xlabel="Time step", ylabel="Coherence m^mu(t)", title="Evolution of Coherence, p=" * string(num_patterns), left_margin=10mm, bottom_margin=5mm, top_margin=5mm, right_margin=5mm, color=[:orange :blue], legend=:topleft)

# Save the plot to a file
savefig(plot_obj, "p2_1.png")
