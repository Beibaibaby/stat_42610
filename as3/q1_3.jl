using Distributions, Plots, Measures, ProgressMeter

# Parameters
C = 1.0 # µF/cm^2
gL = 0.1 # mS/cm^2
VL = -65.0 # mV
VT = -60.0 # mV
Vreset = -70.0 # mV
Vpeak = -45.0 # mV
Delta = 2.0 # mV
aE = 0.2 # µA/cm^2, strength of excitatory input
lambdaE = 500 # Hz, rate of excitatory Poisson input
T = 1000 # ms, simulation time for each realization
dt = 1 # ms, time step
N_realizations = 1000 # Number of realizations

# Synaptic input as in part (b), inhibition set to zero, so we only simulate IE
tau1E = 4.0 # ms
tau2E = 0.4 # ms

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

function simulate_neuron()
    time = 0:dt:T
    V = VL * ones(length(time))
    spikes = generate_poisson_spikes(lambdaE, T, dt)
    IE = zeros(length(time))
    spike_counts = 0

    for i in 2:length(time)
        t = time[i]
        IE[i] = aE * sum([x(t, tau1E, tau2E, time[j]) for j in findall(spikes .== 1)])
        dVdt = (-gL * (V[i-1] - VL) + phi(V[i-1]) + IE[i]) / C
        V[i] = V[i-1] + dVdt * dt

        if V[i] >= Vpeak
            V[i] = Vreset
            spike_counts += 1
        end
    end

    return V, IE, spike_counts
end
using Distributions, Plots

# Run a single simulation for plotting
V, IE, _ = simulate_neuron()
time = 0:dt:T

# Create a plot with two subplots
p = plot(layout=(2, 1), xlabel="Time (ms)", legend=false)

# Plot membrane potential V(t) in the first subplot
plot!(p[1], time, V, label="V(t)", color=:blue, ylabel="V (mV)", title="Membrane Potential",left_margin=15mm,bottom_margin=15mm)

# Plot synaptic current I(t) in the second subplot
plot!(p[2], time, IE, label="I(t)", color=:red, ylabel="I (µA/cm²)", title="Synaptic Current",left_margin=15mm,bottom_margin=15mm)

# Display the combined plot
#display(p)

# Save the figure
savefig("p1_3.png")

# Run N_realizations simulations with a progress bar
spike_counts = zeros(N_realizations)
@showprogress for n in 1:N_realizations
    _, _, spike_count = simulate_neuron()
    spike_counts[n] = spike_count
end

# Compute Fano Factor
M = mean(spike_counts)
Var = var(spike_counts)
FF = Var / M

println("Fano Factor: $FF")