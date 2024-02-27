using Distributions, Plots, Measures, ProgressMeter

# Parameters
C = 1.0 # µF/cm^2
gL = 0.1 # mS/cm^2
VL = -65.0 # mV
VT = -60.0 # mV
Vreset = -70.0 # mV
Vpeak = -45.0 # mV
Delta = 2.0 # mV
lambdaE = 500 # Hz
lambdaI = 500 # Hz
T = 5000 # ms
dt = 0.2 # ms
N_realizations = 1000

# Synaptic input parameters for E and I
tau1E = 4.0 # ms
tau2E = 0.4 # ms
tau1I = 6.0 # ms
tau2I = 1.75 # ms

# Functions phi, x, generate_poisson_spikes remain unchanged

function phi(V)
    return gL * Delta * exp((V - VT) / Delta)
end

function x(t, tau1, tau2, t_hat)
    if t < t_hat
        return 0
    else
        return exp(-(t-t_hat)/tau1) - exp(-(t-t_hat)/tau2)
    end
end

function generate_poisson_spikes(rate, T, dt)
    lambda_dt = rate * dt / 1000
    return rand(Bernoulli(lambda_dt), length(0:dt:T))
end

# Adjusted simulate_neuron function to include aE, aI parameters
function simulate_neuron(aE, aI,dt,T)
    time = 0:dt:T
    V = VL * ones(length(time))
    spikesE = generate_poisson_spikes(lambdaE, T, dt)
    spikesI = generate_poisson_spikes(lambdaI, T, dt)
    IE = zeros(length(time))
    II = zeros(length(time))
    spike_counts = 0

    for i in 2:length(time)
        t = time[i]
        IE[i] = aE * sum([x(t, tau1E, tau2E, time[j]) for j in findall(spikesE .== 1)])
        II[i] = aI * sum([x(t, tau1I, tau2I, time[j]) for j in findall(spikesI .== 1)])
        I_total = IE[i] + II[i]
        dVdt = (-gL * (V[i-1] - VL) + phi(V[i-1]) + I_total) / C
        V[i] = V[i-1] + dVdt * dt

        if V[i] >= Vpeak
            V[i] = Vreset
            spike_counts += 1
        end
    end

    return V, IE, II, spike_counts / (T / 1000) # Return firing rate in Hz
end

# Global scope
optimal_aE = 0.0
optimal_aI = 0.0
closest_rate = 1000000  # Initialize with a high value to ensure it gets updated properly
optimal_rate = 0

# Loop through aE and aI values
for aE in 0.80:0.01:1.0
    aI = -aE * (tau1E - tau2E) / (tau1I - tau2I)
    _, _,_, firing_rate = simulate_neuron(aE, aI,dt, T)  # Ensure T and simulate_neuron are defined
    
    println("aE: $aE, Firing Rate: $firing_rate")
    println("Best Firing Rate: $closest_rate")
    # Update global variables if the condition is met
    if abs(firing_rate - 13) < abs(closest_rate - 13)
        global closest_rate = firing_rate
        global optimal_aE = aE
        global optimal_aI = aI
        global optimal_rate = firing_rate
    end
end

# After the loop, you can use optimal_aE, optimal_aI, and optimal_rate for further processing



println("Optimal aE: $optimal_aE, Optimal aI: $optimal_aI, Firing Rate: $optimal_rate Hz")
T = 1000 # ms
dt = 1 # ms # time step in ms

# Compute Fano Factor for the optimal values
spike_counts = zeros(N_realizations)
@showprogress "Computing Fano Factor..." for n in 1:N_realizations
    _, _,_, spike_count = simulate_neuron(optimal_aE, optimal_aI,dt,T)
    spike_counts[n] = spike_count
end


M = mean(spike_counts)
Var = var(spike_counts)
FF = Var / M

println("Fano Factor: $FF")


using Plots  # Ensure the Plots package is used

# Assuming the simulation loop is set up correctly
while true
    T = 1000 # ms
    dt = 0.1 # time step in ms
    time = 0:dt:T
    V, IE, II, spike_count = simulate_neuron(optimal_aE, optimal_aI, dt,T)
    println("Spike Count: $spike_count")

    # When spike_count equals 13, proceed to plot
    if spike_count == 13
        # Calculate total current for plotting
        I_total = IE + II

        # Plotting the currents (IE, II, I_total) in one subplot
        currents_plot = plot(time, IE, label="I_E", color=:red, linewidth=2,left_margin=15mm,bottom_margin=15mm)
        plot!(currents_plot, time, II, label="I_I", color=:green, linewidth=2)
        plot!(currents_plot, time, I_total, label="I_total", color=:black, linewidth=2)
        title!(currents_plot, "Synaptic Currents Over Time")
        xlabel!(currents_plot, "Time (ms)")
        ylabel!(currents_plot, "Current (µA/cm^2)")

        # Plotting the membrane potential V in another subplot
        potential_plot = plot(time, V, label="Membrane Potential V", color=:blue, linewidth=2,left_margin=15mm,bottom_margin=15mm)
        title!(potential_plot, "Membrane Potential Over Time")
        xlabel!(potential_plot, "Time (ms)")
        ylabel!(potential_plot, "Potential (mV)")

        # Combining the two plots into a single figure with subplots arranged vertically
        combined_plot = plot(currents_plot, potential_plot, layout=(2, 1), size=(1000, 800), legend=:outertop)

        savefig(combined_plot, "p1_5.png")

        break
    end

    # Optional: Include logic for adjusting parameters or conditions to ensure the loop can eventually meet the exit condition
end
