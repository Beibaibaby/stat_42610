using Statistics
using Plots
using Measures  # For specifying margins in mm

doplot = true

include("sim.jl")  # Ensure this sim.jl script is correctly defined in your working directory
T = 1100
times, ns, Ne, Ncells, T = sim(1000,T)

println("mean excitatory firing rate: ", mean(1000 * ns[1:Ne] / T), " Hz")
println("mean inhibitory firing rate: ", mean(1000 * ns[(Ne + 1):Ncells] / T), " Hz")

if doplot
    println("creating plot")
    p = plot(layout=(1,1), size=(600,400))  # Initialize an empty plot with specified size
    
    # Excitatory neurons
    for ci = 1:Ne
        vals = times[ci, 1:ns[ci]]
        scatter!(vals, fill(ci, length(vals)), label=false, markersize=0.3, markerstrokewidth=0, color=:red)
    end
    
    # Inhibitory neurons
    for ci = Ne+1:Ncells
        vals = times[ci, 1:ns[ci]]
        scatter!(vals, fill(ci, length(vals)), label=false, markersize=0.3, markerstrokewidth=0, color=:blue)
    end

    xlims!(0, T)
    ylims!(0, Ncells)
    ylabel!("Neuron")
    xlabel!("Time (ms)")
    title!("Neuronal Firing Raster Plot")

    plot!(margin=5mm)  # Set margins

    savefig("p1a.png")
end



# Assume `times` and `ns` are obtained from the `sim` function
# Gaussian kernel function
function gaussian_kernel(x, σ)
    return exp(-x^2 / (2σ^2)) / (σ * sqrt(2π))
end

σ_r = 3  # Standard deviation for the Gaussian kernel
dt = 0.1  # Time step used in simulation, for binning
T_end = T  # Total simulation time

# Initialize the firing rate array
time_bins = 0:dt:T
R_E_t = zeros(length(time_bins)-1)
# Precompute the Gaussian kernel over a range of time differences
max_time_diff = 3 * σ_r  # Typically, 3σ covers most of the relevant range
time_diffs = -max_time_diff:dt:max_time_diff
gaussian_kernel_values = exp.(-time_diffs.^2 / (2 * σ_r^2)) / (σ_r * sqrt(2π))

# Initialize the firing rate array
time_bins = 0:dt:T_end
R_E_t = zeros(length(time_bins)-1)

# Apply the precomputed Gaussian kernel to the spike trains
for ci = 1:Ne
    spike_times = times[ci, times[ci, :] .> 0]
    for spike_time in spike_times
        spike_index = min(max(1, Int(floor(spike_time / dt)) + 1), length(R_E_t))
        # Find the range of indices in R_E_t to update
        start_index = max(1, spike_index - length(time_diffs) ÷ 2)
        end_index = min(length(R_E_t), spike_index + length(time_diffs) ÷ 2)
        kernel_start_index = max(1, 1 + length(time_diffs) ÷ 2 - spike_index)
        kernel_end_index = kernel_start_index + end_index - start_index
        # Add the kernel values to the firing rate array
        R_E_t[start_index:end_index] .+= gaussian_kernel_values[kernel_start_index:kernel_end_index]
    end
end

R_E_t *= 1000 / Ne  # Convert to firing rate

plot_start = 50  # Start plotting from 50 ms
plot_end = T_end - 50  # Stop plotting 50 ms before the end
start_index = Int(ceil(plot_start / dt))
end_index = Int(floor(plot_end / dt))

# Compute M^E and CV^E
M_E = mean(R_E_t[start_index:end_index])
CV_E = std(R_E_t[start_index:end_index]) / M_E


if doplot
    # Define the time range to exclude first and last 50 ms from the plot
    plot_start = 50  # Start plotting from 50 ms
    plot_end = T_end - 50  # Stop plotting 50 ms before the end

    # Find the indices corresponding to the plot start and end times
    start_index = Int(ceil(plot_start / dt))
    end_index = Int(floor(plot_end / dt))

    # Plot R^E(t) within the specified range
    p2 = plot(time_bins[start_index:end_index], R_E_t[start_index:end_index], label="R^E(t)", color=:red,left_margin=10mm)
    title!("E Neuron Firing Rate Over Time")
    xlabel!("Time (ms)")
    ylabel!("Firing Rate (Hz)")
    plot!(margin=5mm)  # Set margins

    savefig("p1b.png")
end



########
#######
######Longer Stim for Statistics
T=10100
times, ns, Ne, Ncells, T = sim(1000,T)
println("mean excitatory firing rate: ", mean(1000 * ns[1:Ne] / T), " Hz")
println("mean inhibitory firing rate: ", mean(1000 * ns[(Ne + 1):Ncells] / T), " Hz")


σ_r = 3  # Standard deviation for the Gaussian kernel
dt = 0.1  # Time step used in simulation, for binning
T_end = T  # Total simulation time

# Initialize the firing rate array
time_bins = 0:dt:T
R_E_t = zeros(length(time_bins)-1)
# Precompute the Gaussian kernel over a range of time differences
max_time_diff = 3 * σ_r  # Typically, 3σ covers most of the relevant range
time_diffs = -max_time_diff:dt:max_time_diff
gaussian_kernel_values = exp.(-time_diffs.^2 / (2 * σ_r^2)) / (σ_r * sqrt(2π))

# Initialize the firing rate array
time_bins = 0:dt:T_end
R_E_t = zeros(length(time_bins)-1)

# Apply the precomputed Gaussian kernel to the spike trains
for ci = 1:Ne
    spike_times = times[ci, times[ci, :] .> 0]
    for spike_time in spike_times
        spike_index = min(max(1, Int(floor(spike_time / dt)) + 1), length(R_E_t))
        # Find the range of indices in R_E_t to update
        start_index = max(1, spike_index - length(time_diffs) ÷ 2)
        end_index = min(length(R_E_t), spike_index + length(time_diffs) ÷ 2)
        kernel_start_index = max(1, 1 + length(time_diffs) ÷ 2 - spike_index)
        kernel_end_index = kernel_start_index + end_index - start_index
        # Add the kernel values to the firing rate array
        R_E_t[start_index:end_index] .+= gaussian_kernel_values[kernel_start_index:kernel_end_index]
    end
end

R_E_t *= 1000 / Ne  # Convert to firing rate

plot_start = 50  # Start plotting from 50 ms
plot_end = T_end - 50  # Stop plotting 50 ms before the end
start_index = Int(ceil(plot_start / dt))
end_index = Int(floor(plot_end / dt))

# Compute M^E and CV^E
M_E = mean(R_E_t[start_index:end_index])
CV_E = std(R_E_t[start_index:end_index]) / M_E


println("M^E: $M_E Hz, CV^E: $CV_E")





