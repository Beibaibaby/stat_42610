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

# Function to simulate the network dynamics with I2 = I + deltaI
function simulate_network_dynamics(I, deltaI)
    I2 = I + deltaI  # Define I2 based on I and deltaI
    r1, a1, r2, a2 = 0.1, 0.0, 0.0, 0.0
    r1_values = zeros(length(time))

    for t in 1:length(time)
        dr1dt = (-r1 + S(I - g * a1 - w * r2)) / tau_r
        da1dt = (-a1 + r1) / tau_a
        dr2dt = (-r2 + S(I2 - g * a2 - w * r1)) / tau_r
        da2dt = (-a2 + r2) / tau_a

        r1 += dt * dr1dt
        a1 += dt * da1dt
        r2 += dt * dr2dt
        a2 += dt * da2dt

        r1_values[t] = r1
    end
    return r1_values
end

deltaI_values = -0.2:0.01:0.2  # Range of deltaI values
range_of_I1 = []  # To store the range of I1 for each deltaI

for deltaI in deltaI_values
    max_values = []
    min_values = []

    # Compute max and min for each I value
    for I in I_range
        r1_values = simulate_network_dynamics(I, deltaI)
        start_index = Int(600 / dt) + 1
        end_index = Int(1100 / dt)
        segment = r1_values[start_index:end_index]
        push!(max_values, maximum(segment))
        push!(min_values, minimum(segment))
    end

    # Calculate the range of I1 as the max of max_values - min_values
    difference = max_values .- min_values
    push!(range_of_I1, maximum(difference))
end

# Plotting the range of I1 dependent on deltaI
plot(deltaI_values, range_of_I1, xlabel="ΔI", ylabel="Range of I1", title="Range of I1 vs. ΔI", label=false, color=:blue,left_margin=10mm, bottom_margin=5mm, top_margin=5mm, right_margin=5mm)
# Save the plot to a file
savefig("p1_e2.png")
