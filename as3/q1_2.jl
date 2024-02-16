using Distributions
using Plots
using Measures
# Parameters
aE = 1.0 # µA/cm^2 for excitatory synapses
lambdaE = 500 # Hz
lambdaI = 500 # Hz
tau1E = 4.0 # ms
tau2E = 0.4 # ms
tau1I = 6.0 # ms
tau2I = 1.75 # ms

# Time parameters
T = 5000 # total simulation time in ms
dt = 0.1 # time step in ms
time = 0:dt:T

# Compute aI using the balance condition
aI = -aE * (tau1E - tau2E) / (tau1I - tau2I)

# Function to generate Poisson spike trains
function generate_poisson_spikes(rate, T, dt)
    lambda_dt = rate * dt / 1000 # Convert rate to spikes per time step
    return rand(Bernoulli(lambda_dt), length(0:dt:T))
end

# Function to compute x(t) for a single synaptic input
function x(t, tau1, tau2, t_hat)
    if t < t_hat
        return 0
    else
        return exp(-(t-t_hat)/tau1) - exp(-(t-t_hat)/tau2)
    end
end

# Synaptic current computation
function compute_synaptic_current(spike_train, a, tau1, tau2, time)
    I = zeros(length(time))
    for t_hat_index in findall(spike_train .== 1)
        t_hat = time[t_hat_index]
        for t_index in t_hat_index:length(time)
            t = time[t_index]
            I[t_index] += a * x(t, tau1, tau2, t_hat)
        end
    end
    return I
end

# Generate spike trains
spikesE = generate_poisson_spikes(lambdaE, T, dt)
spikesI = generate_poisson_spikes(lambdaI, T, dt)

# Compute synaptic currents
IE = compute_synaptic_current(spikesE, aE, tau1E, tau2E, time)
II = compute_synaptic_current(spikesI, aI, tau1I, tau2I, time)

# Total current
I_total = IE + II

# Verify balance
mean_current = mean(I_total)
println("Mean total current: $mean_current µA/cm^2")

# Generate the plot
plot(time, IE, label="I_E", xlabel="Time (ms)", ylabel="Current (µA/cm^2)", title="Synaptic Currents Over Time", linewidth=2,size=(2000,500),left_margin=20mm,bottom_margin=15mm)
plot!(time, II, label="I_I", linewidth=2)
plot!(time, I_total, label="I_total", linewidth=2, linestyle=:dash)

# Display the plot (optional, depending on your environment)
#display(plot)

# Save the figure to a file
savefig("p1_2.png")

function simulate_synaptic_currents()
    # Parameters
    aE = 1.0 # µA/cm^2 for excitatory synapses
    lambdaE = 500 # Hz
    lambdaI = 500 # Hz
    tau1E = 4.0 # ms
    tau2E = 0.4 # ms
    tau1I = 6.0 # ms
    tau2I = 1.75 # ms

    # Time parameters
    T = 1000 # total simulation time in ms
    dt = 0.1 # time step in ms
    time = 0:dt:T

    # Compute aI using the balance condition
    aI = -aE * (tau1E - tau2E) / (tau1I - tau2I)

    # Function to generate Poisson spike trains
    function generate_poisson_spikes(rate, T, dt)
        lambda_dt = rate * dt / 1000 # Convert rate to spikes per time step
        return rand(Bernoulli(lambda_dt), length(0:dt:T))
    end

    # Function to compute x(t) for a single synaptic input
    function x(t, tau1, tau2, t_hat)
        if t < t_hat
            return 0
        else
            return exp(-(t-t_hat)/tau1) - exp(-(t-t_hat)/tau2)
        end
    end

    # Synaptic current computation
    function compute_synaptic_current(spike_train, a, tau1, tau2, time)
        I = zeros(length(time))
        for t_hat_index in findall(spike_train .== 1)
            t_hat = time[t_hat_index]
            for t_index in t_hat_index:length(time)
                t = time[t_index]
                I[t_index] += a * x(t, tau1, tau2, t_hat)
            end
        end
        return I
    end

    # Generate spike trains
    spikesE = generate_poisson_spikes(lambdaE, T, dt)
    spikesI = generate_poisson_spikes(lambdaI, T, dt)

    # Compute synaptic currents
    IE = compute_synaptic_current(spikesE, aE, tau1E, tau2E, time)
    II = compute_synaptic_current(spikesI, aI, tau1I, tau2I, time)

    # Total current
    I_total = IE + II

    # Return mean current
    return mean(I_total)
end

# Run simulation 50 times
mean_currents = [simulate_synaptic_currents() for _ in 1:100]

# Calculate average of the mean currents
average_mean_current = mean(mean_currents)

println("Average of the mean total currents over 100 simulations: $average_mean_current µA/cm^2")
