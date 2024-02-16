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
T = 1000 # ms
dt = 1 # ms
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
function simulate_neuron(aE, aI)
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

    return V, IE + II, spike_counts / (T / 1000) # Return firing rate in Hz
end

# Initialize variables for finding the optimal aE and aI
global optimal_aE = 0.0
global optimal_aI = 0.0
global closest_rate = 1000000
global optimal_rate = 0

# Loop through aE and aI values without @showprogress for troubleshooting
for aE in 0.1:0.05:1.0
    aI = -aE
    _, _, firing_rate = simulate_neuron(aE, aI)
    println(aE)
    println(firing_rate)

    global closest_rate, optimal_aE, optimal_aI, optimal_rate
    if abs(firing_rate - 13) < 2
        closest_rate = abs(firing_rate - 13)
        optimal_aE = aE
        optimal_aI = aI
        optimal_rate = firing_rate
    end
end


println("Optimal aE: $optimal_aE, Optimal aI: $optimal_aI, Firing Rate: $optimal_rate Hz")

# Compute Fano Factor for the optimal values
spike_counts = zeros(N_realizations)
@showprogress "Computing Fano Factor..." for n in 1:N_realizations
    _, _, spike_count = simulate_neuron(optimal_aE, optimal_aI)
    spike_counts[n] = spike_count
end

M = mean(spike_counts)
Var = var(spike_counts)
FF = Var / M

println("Fano Factor: $FF")
