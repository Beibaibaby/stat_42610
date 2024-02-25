using Plots
using Measures

# Define parameters
tau_r = 1.0
tau_a = 100.0
theta = -0.1
K = 0.1
w = 0.3
g = 0.5
dt = 0.1 # Time step for the Euler method
T = 1000 # Total time for simulation
time = 0:dt:T
input_conditions = [(0.2, 0.7), (0.5, 0.5)]

# Sigmoid function
S(x) = 1 / (1 + exp(-(x - theta) / K))

# Initialize plots
plots = []

for (I1, I2) in input_conditions
    # Initial conditions
    r1 = 0.1
    a1 = 0.0
    r2 = 0.0
    a2 = 0.0

    # Arrays to store the time evolution of r1 and r2
    r1_values = zeros(length(time))
    r2_values = zeros(length(time))

    for t in 1:length(time)
        # Euler method to update r1, a1, r2, a2
        dr1dt = (-r1 + S(I1 - g*a1 - w*r2)) / tau_r
        da1dt = (-a1 + r1) / tau_a
        dr2dt = (-r2 + S(I2 - g*a2 - w*r1)) / tau_r
        da2dt = (-a2 + r2) / tau_a
        r1 += dt * dr1dt
        a1 += dt * da1dt
        r2 += dt * dr2dt
        a2 += dt * da2dt

        # Store the updated values
        r1_values[t] = r1
        r2_values[t] = r2
    end

    # Plot r1(t) and r2(t) for the current input condition
    p = plot(time, r1_values, label="r1(t) for I1=$(I1), I2=$(I2)", title="r1(t) and r2(t), w=0.3",left_margin=10mm, bottom_margin=5mm, top_margin=5mm, right_margin=5mm)
    plot!(p, time, r2_values, label="r2(t) for I1=$(I1), I2=$(I2)")
    push!(plots, p)
end

# Combine all plots into one figure
final_plot = plot(plots..., layout=(length(input_conditions), 1), size=(800, 600))

# Save the plot to a file
savefig(final_plot, "p1_c_1.png")



# Define parameters
tau_r = 1.0
tau_a = 100.0
theta = -0.1
K = 0.1
w = 0.0
g = 0.5
dt = 0.1 # Time step for the Euler method
T = 1000 # Total time for simulation
time = 0:dt:T
input_conditions = [(0.2, 0.7), (0.5, 0.5)]

# Sigmoid function
S(x) = 1 / (1 + exp(-(x - theta) / K))

# Initialize plots
plots = []

for (I1, I2) in input_conditions
    # Initial conditions
    r1 = 0.1
    a1 = 0.0
    r2 = 0.0
    a2 = 0.0

    # Arrays to store the time evolution of r1 and r2
    r1_values = zeros(length(time))
    r2_values = zeros(length(time))

    for t in 1:length(time)
        # Euler method to update r1, a1, r2, a2
        dr1dt = (-r1 + S(I1 - g*a1 - w*r2)) / tau_r
        da1dt = (-a1 + r1) / tau_a
        dr2dt = (-r2 + S(I2 - g*a2 - w*r1)) / tau_r
        da2dt = (-a2 + r2) / tau_a
        r1 += dt * dr1dt
        a1 += dt * da1dt
        r2 += dt * dr2dt
        a2 += dt * da2dt

        # Store the updated values
        r1_values[t] = r1
        r2_values[t] = r2
    end

    # Plot r1(t) and r2(t) for the current input condition
    p = plot(time, r1_values, label="r1(t) for I1=$(I1), I2=$(I2)", title="r1(t) and r2(t), w=0.0",left_margin=10mm, bottom_margin=5mm, top_margin=5mm, right_margin=5mm)
    plot!(p, time, r2_values, label="r2(t) for I1=$(I1), I2=$(I2)")
    push!(plots, p)
end

# Combine all plots into one figure
final_plot = plot(plots..., layout=(length(input_conditions), 1), size=(800, 600))

# Save the plot to a file
savefig(final_plot, "p1_c_2.png")