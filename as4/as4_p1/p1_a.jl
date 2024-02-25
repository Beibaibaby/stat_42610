using Plots
using Measures
# Define parameters
tau_r = 1.0
tau_a = 100.0
theta = -0.1
K = 0.1
dt = 0.1 # Time step for the Euler method
T = 1000 # Total time for simulation
time = 0:dt:T
conditions = [(0.2, 0.1), (0.7, 0.1), (0.2, 0.7), (0.7, 0.7)]

# Sigmoid function
S(x) = 1 / (1 + exp(-(x - theta) / K))

# Initialize an empty array to store plots
plots = []

for (I, g) in conditions
    # Initial conditions
    r = 0.0
    a = 0.0

    # Arrays to store the time evolution of r and a
    r_values = zeros(length(time))
    a_values = zeros(length(time))

    for t in 1:length(time)
        # Euler method to update r and a
        drdt = (-r + S(I - g*a)) / tau_r
        dadt = (-a + r) / tau_a
        r += dt * drdt
        a += dt * dadt

        # Store the updated values
        r_values[t] = r
        a_values[t] = a
    end

    # Plot r(t) and a(t) for the current condition
    p = plot(time, r_values, label="r(t)", title="I=$(I), g=$(g)", xlabel="Time", ylabel="Value",left_margin=5mm, bottom_margin=5mm, top_margin=5mm, right_margin=5mm)
    plot!(p, time, a_values, label="a(t)")
    push!(plots, p)
end

# Combine all plots into one figure
final_plot = plot(plots..., layout=(2, 2), size=(800, 600))

# Save the plot to a file
savefig(final_plot, "r_a_dynamics.png")
