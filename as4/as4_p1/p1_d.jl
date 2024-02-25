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
I_range = 0.0:0.02:1.0

# Sigmoid function
S(x) = 1 / (1 + exp(-(x - theta) / K))

# Function to simulate the network dynamics
function simulate_network_dynamics(I)
    r1, a1, r2, a2 = 0.1, 0.0, 0.0, 0.0
    r1_values, r2_values = zeros(length(time)), zeros(length(time))

    for t in 1:length(time)
        dr1dt = (-r1 + S(I - g * a1 - w * r2)) / tau_r
        da1dt = (-a1 + r1) / tau_a
        dr2dt = (-r2 + S(I - g * a2 - w * r1)) / tau_r
        da2dt = (-a2 + r2) / tau_a

        r1 += dt * dr1dt
        a1 += dt * da1dt
        r2 += dt * dr2dt
        a2 += dt * da2dt

        r1_values[t], r2_values[t] = r1, r2
    end
    return r1_values, r2_values
end

# Adjusted function to detect periodicity in r1's activity
function detect_periodicity(r1_values)
    # Find local maxima by checking if each point is greater than its neighbors
    peaks = [i for i in 2:length(r1_values)-1 if r1_values[i] > r1_values[i-1] && r1_values[i] > r1_values[i+1]]

    if length(peaks) > 10 # Ensure there are enough peaks to consider it periodic
        intervals = diff(peaks)
        mean_interval = mean(intervals)
        std_interval = std(intervals)
        # Check for regularity in peak intervals
        if std_interval < mean_interval * 0.5 # Adjust threshold for regularity as needed
            return true
        end
    end
    return false
end

periodic_I_values = [] # To store I values with periodic dynamics
example_plot_created = false # To check if the example plot is created

# Main loop to explore I values and detect periodic dynamics
# Main loop to explore I values and detect periodic dynamics
for I in I_range
    global example_plot_created # Add this line to declare the variable as global within the loop
    r1_values, r2_values = simulate_network_dynamics(I)
    if detect_periodicity(r1_values)
        push!(periodic_I_values, I)
        println("Periodic dynamics detected at I=$I")
        
            # Plot r1 and r2 for the first detected I value
            p = plot(time, r1_values, label="r1(t)", title="Periodic Dynamics for I=$I", xlabel="Time", ylabel="Firing Rate", color=:blue,left_margin=10mm, bottom_margin=5mm, top_margin=5mm, right_margin=5mm)
            plot!(p, time, r2_values, label="r2(t)", color=:red)
            savefig(p, "Periodic_Dynamics_I_$I.png")
            
        
    end
end


if !isempty(periodic_I_values)
    println("Periodic dynamics detected at I values: ", join(periodic_I_values, ", "))
else
    println("No periodic dynamics detected within the given range of I.")
end


