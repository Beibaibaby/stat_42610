using Plots
using Statistics
using Measures

# Define parameters
tau_r = 1.0
tau_a = 100.0
theta = -0.1
K = 0.1
w = 0.45
g = 0.5
dt = 0.1
T = 1500
time = 0:dt:T
I_range = 0.0:0.01:1.0

# Sigmoid function
S(x) = 1 / (1 + exp(-(x - theta) / K))

# Function to simulate the network dynamics
function simulate_network_dynamics(I)
    r1, a1, r2, a2 = 0.1, 0.0, 0.0, 0.0
    r1_values = zeros(length(time))

    for t in 1:length(time)
        dr1dt = (-r1 + S(I - g * a1 - w * r2)) / tau_r
        da1dt = (-a1 + r1) / tau_a
        dr2dt = (-r2 + S(I - g * a2 - w * r1)) / tau_r
        da2dt = (-a2 + r2) / tau_a

        r1 += dt * dr1dt
        a1 += dt * da1dt
        r2 += dt * dr2dt
        a2 += dt * da2dt

        r1_values[t] = r1
    end
    return r1_values
end

# Arrays to store max and min values for each I
max_values = []
min_values = []

# Time indices for the range 600ms to 1100ms
start_index = Int(600 / dt) + 1
end_index = Int(1100 / dt)

# Main loop to compute max and min for each I value
for I in I_range
    r1_values = simulate_network_dynamics(I)
    # Extract the segment of interest
    segment = r1_values[start_index:end_index]
    # Compute max and min
    push!(max_values, maximum(segment))
    push!(min_values, minimum(segment))
end

# Plotting
p = plot(I_range, max_values, label="Max r1(t)", xlabel="I", ylabel="r1(t) Value", title="Max and Min of r1(t) vs. I", color=:orange)
plot!(p, I_range, min_values, label="Min r1(t)", color=:purple)
# Adding margins for better visualization
plot!(p, left_margin=10mm, bottom_margin=5mm, top_margin=5mm, right_margin=5mm)

# Save the plot to a file
savefig(p, "p1_d.png")
