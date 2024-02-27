using Plots
using Distributions


# Function to initialize the network parameters and state variables
function initialize_network(N_E, N_I, V_rest, V_thresh, tau, J_EE, J_IE, J_EI, J_II, J_E0, J_I0, p)
    N = N_E + N_I  # Total number of neurons
    V = V_rest * ones(N)  # Membrane potentials
    spiking_neurons = []
    connections = rand(Bernoulli(p), N, N)  # Random connections between neurons
    J = zeros(N, N)  # Synaptic strengths
    K=N_E*p
    J[1:N_E, 1:N_E] .= J_EE / sqrt(K)  # EE connections
    J[N_E+1:N, 1:N_E] .= J_IE / sqrt(K)  # IE connections
    J[1:N_E, N_E+1:N] .= J_EI / sqrt(K)  # EI connections
    J[N_E+1:N, N_E+1:N] .= J_II / sqrt(K)  # II connections
    J[:, 1:N_E] .+= J_E0 / sqrt(K)  # External input to E population
    J[:, N_E+1:N] .+= J_I0 / sqrt(K)  # External input to I population
    return V, spiking_neurons, connections, J
end

# Function to update the network state using the forward Euler method
function update_network!(V, spiking_neurons, connections, J, dt, V_rest, V_thresh, tau)
    N = length(V)
    dV = (-V + connections * J) / tau * dt  # dV/dt update
    V .+= dV  # Update membrane potentials
    for i in 1:N
        if V[i] >= V_thresh  # Check if the neuron has fired
            push!(spiking_neurons, (i, t))  # Record the spike
            V[i] = V_rest  # Reset membrane potential
        end
    end
end

# Parameters
N_E = 1000
N_I = 1000
V_rest = 0.0
V_thresh = 1.0
tau = 15e-3  # 15 ms
J_EE = 1.0
J_IE = 2.0
J_EI = 3.0
J_II = 2.5
J_E0 = 1.2
J_I0 = 0.7
p = 0.2
dt = 1e-3  # Time step of 1 ms
T = 1.0  # Total simulation time of 1000 ms

# Initialize the network
V, spiking_neurons, connections, J = initialize_network(N_E, N_I, V_rest, V_thresh, tau, J_EE, J_IE, J_EI, J_II, J_E0, J_I0, p)

# Run the simulation
t = 0.0
while t < T
    update_network!(V, spiking_neurons, connections, J, dt, V_rest, V_thresh, tau)
    t += dt
end

# Generate the raster plot
plot()
for (neuron, spike_time) in spiking_neurons
    scatter!([spike_time], [neuron], color=:black, markersize=1, legend=false)
end
xlabel!("Time (ms)")
ylabel!("Neuron index")
title!("Raster Plot of Network Activity")

# Save the plot
savefig("raster_plot.png")
