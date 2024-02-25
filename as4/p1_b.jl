using Plots
using Measures
# Define parameters
tau_r = 1.0
tau_a = 100.0
theta = -0.1 # Updated theta value
K = 0.1
dt = 0.1 # Time step for the Euler method
T = 2000 # Total time for simulation, adjust if needed to ensure steady-state
time = 0:dt:T

# Sigmoid function
S(x) = 1 / (1 + exp(-(x - theta) / K))

# Grid of I and g values
I_values = 0:0.05:1
g_values = 0:0.05:1
rs_steady = zeros(length(I_values), length(g_values))

for (i, I) in enumerate(I_values)
    for (j, g) in enumerate(g_values)
        # Initial conditions
        r = 0.0
        a = 0.0
        
        for t in 1:length(time)
            # Euler method to update r and a
            drdt = (-r + S(I - g*a)) / tau_r
            dadt = (-a + r) / tau_a
            r += dt * drdt
            a += dt * dadt
        end
        
        # Assume r at the end of the simulation is the steady-state firing rate
        rs_steady[i, j] = r
    end
end

# Plot the steady-state firing rate over the (I, g) plane
heatmap(I_values, g_values, rs_steady', xlabel="I", ylabel="g", title="Steady-State Firing Rate (r_s) over (I, g) plane", c=:viridis,left_margin=15mm, bottom_margin=5mm, top_margin=5mm, right_margin=5mm)

# Save the plot
savefig("p1_b.png")
