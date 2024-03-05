using Statistics
using Plots
using Measures  # For specifying margins in mm
ENV["GKSwstype"] = "100"  # Use the off-screen GKS terminal
# Function to generate a spike based on a homogeneous Poisson process
function generate_homogeneous_poisson_spike(λ, dt)
    u = rand() # Uniform random variable U(0, 1)
    return u < λ * dt
end

function sim(N,T)
	println("setting up parameters")
	
	Ne = N
	Ni = N
	Ncells= Ne + Ni


	taue = 10 #membrane time constant for exc. neurons (ms)
	taui = 10
    
	p = 0.2
	#connection probabilities

	K = p * N #average number in-pop connections per neuron
	sqrtK = sqrt(K)

	jie = 1.25 /sqrtK
	jee = 1.0 /sqrtK

	jei = -1.25 /sqrtK
	jii = -1.0 /sqrtK

	jeX = 2.0 / sqrtK 
	jiX = 0.8 / sqrtK

	stimstart = 0
	stimend = T 


	vre = 0.0 #reset voltage

	threshe = 1 #threshold for exc. neurons
	threshi = 1

	dt = 0.1 #simulation timestep (ms)
	refrac = 0.0 #refractory period (ms)

	maxrate = 1000 #(Hz) maximum average firing rate.  if the average firing rate across the simulation for any neuron exceeds this value, some of that neuron's spikes will not be saved
    

	# Parameters for the λ(t) dynamics
    lambda_h = K / 100  # Base firing rate as before
    tau_lambda = 25 * taue  # Given τ_λ

    sigma_lambda = 3 * lambda_h / sqrt(K)  # Given σ_λ

    # Initial value of λ(t) for each neuron
    lambda_t = lambda_h

	mu = zeros(Ncells)
	thresh = zeros(Ncells)
	tau = zeros(Ncells)

	thresh[1:Ne] .= threshe
	thresh[(1+Ne):Ncells] .= threshi

	tau[1:Ne] .= taue
	tau[(1+Ne):Ncells] .= taui


	weights = zeros(Ncells,Ncells)

	#random connections
	weights[1:Ne,1:Ne] = jee*(rand(Ne,Ne) .< p)
	weights[1:Ne,(1+Ne):Ncells] = jei*(rand(Ne,Ni) .< p)
	weights[(1+Ne):Ncells,1:Ne] = jie*(rand(Ni,Ne) .< p)
	weights[(1+Ne):Ncells,(1+Ne):Ncells] = jii*(rand(Ni,Ni) .< p)
	

	for ci = 1:Ncells
		weights[ci,ci] = 0
	end

	maxTimes = round(Int,maxrate*T/1000)
	times = zeros(Ncells,maxTimes)
	ns = zeros(Int,Ncells)

	forwardInputsE = zeros(Ncells) #summed weight of incoming E spikes
	forwardInputsI = zeros(Ncells)
	forwardInputsEPrev = zeros(Ncells) #as above, for previous timestep
	forwardInputsIPrev = zeros(Ncells)

	xerise = zeros(Ncells) #auxiliary variables for E/I currents (difference of exponentials)
	xedecay = zeros(Ncells)
	xirise = zeros(Ncells)
	xidecay = zeros(Ncells)

	v = rand(Ncells) #membrane voltage 

	lastSpike = -1000*ones(Ncells) #time of last spike

	Nsteps = round(Int,T/dt)

	# Initialize arrays to track external spike times for each neuron

		
	println("starting simulation")
    lambda_ts=zeros(Nsteps)
	#begin main simulation loop
	for ti = 1:Nsteps
		if mod(ti,Nsteps/100) == 1  #print percent complete
			print("\r",round(Int,100*ti/Nsteps))
		end
		t = dt*ti
		forwardInputsE[:] .= 0
		forwardInputsI[:] .= 0
		dW = randn()  # N(0,1) for white noise
		lambda_t += (dt / tau_lambda) * (lambda_h - lambda_t) + sigma_lambda * sqrt(dt / tau_lambda) * dW
		
		lambda_ts[ti]=lambda_t

		for ci = 1:Ncells

			# Generate external spike based on updated λ(t)
			external_spike = generate_homogeneous_poisson_spike(lambda_t, dt)

			external_input = 0

			if external_spike

				if ci < Ne
					external_input = jeX
				else
					external_input = jiX
				end
			end

		    recurrent_input= forwardInputsEPrev[ci]+forwardInputsIPrev[ci]

			if t > (lastSpike[ci] + refrac)  #not in refractory period
				v[ci] += dt*((1/tau[ci])*(-v[ci]) )+recurrent_input+external_input

				if v[ci] > thresh[ci]  #spike occurred
					v[ci] = vre
					lastSpike[ci] = t
					ns[ci] = ns[ci]+1
					if ns[ci] <= maxTimes 
						times[ci,ns[ci]] = t
					end

					for j = 1:Ncells
						if weights[j,ci] > 0  #E synapse
							forwardInputsE[j] += weights[j,ci]
						elseif weights[j,ci] < 0  #I synapse
							forwardInputsI[j] += weights[j,ci]
						end
					end #end loop over synaptic projections
				end #end if(spike occurred)
			end #end if(not refractory)
		end #end loop over neurons

		forwardInputsEPrev = copy(forwardInputsE)
		forwardInputsIPrev = copy(forwardInputsI)
	end #end loop over time
	print("\r")

	times = times[:,1:maximum(ns)]

	return times,ns,Ne,Ncells,T,lambda_ts
end






doplot = true

T = 1100
times, ns, Ne, Ncells, T,lambda_ts = sim(1000,T)

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

    savefig("p2a.png")
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
max_time_diff = 10 * σ_r  # Typically, 3σ covers most of the relevant range
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

CV_lm=std(lambda_ts)/mean(lambda_ts)
println("CV lambda: $CV_lm") 
println("Size of Lambda_ts: $(size(lambda_ts))")
println("Size of R_E_t: $(size(R_E_t))")



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

    savefig("p2b.png")

	# Plot lambda(t) over time
    plot(time_bins[start_index:end_index], lambda_ts[start_index:end_index], label="λ(t) over time", xlabel="Time (ms)", ylabel="λ(t)", title="Evolution of λ(t)",left_margin=10mm)
    savefig("p2c.png")  # Save the plot to a file

end






########
#######
######Longer Stim for Statistics
T=3100
times, ns, Ne, Ncells, T,lambda_ts = sim(1000,T)
println("mean excitatory firing rate: ", mean(1000 * ns[1:Ne] / T), " Hz")
println("mean inhibitory firing rate: ", mean(1000 * ns[(Ne + 1):Ncells] / T), " Hz")


σ_r = 3  # Standard deviation for the Gaussian kernel
dt = 0.1  # Time step used in simulation, for binning
T_end = T  # Total simulation time

# Initialize the firing rate array
time_bins = 0:dt:T
R_E_t = zeros(length(time_bins)-1)
# Precompute the Gaussian kernel over a range of time differences
max_time_diff = 10 * σ_r  # Typically, 3σ covers most of the relevant range
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

CV_lm=std(lambda_ts)/mean(lambda_ts)
println("CV lambda: $CV_lm") 
println("Size of Lambda_ts: $(size(lambda_ts))")
println("Size of R_E_t: $(size(R_E_t))")


p4 = plot(time_bins[start_index:end_index], R_E_t[start_index:end_index], label="R^E(t)", color=:red,left_margin=10mm)
title!("E Neuron Firing Rate Over Time")
xlabel!("Time (ms)")
ylabel!("Firing Rate (Hz)")
plot!(margin=5mm)  # Set margins

savefig("p2d.png")


function simulate_lambda(T, dt, lambda_h, tau_lambda, sigma_lambda)
    Nsteps = round(Int, T/dt)
    lambda_t = zeros(Nsteps)
    lambda = lambda_h  # Initial condition
    
    for ti = 1:Nsteps
        dW = randn() 
        # Update lambda using the Euler-Maruyama method for the SDE
        lambda += (dt / tau_lambda) * (lambda_h - lambda) + sigma_lambda * sqrt(dt / tau_lambda) * dW
        lambda_t[ti] = lambda
    end
    
    return lambda_t
end

function compute_CV_X(lambda_t)
    mean_lambda = mean(lambda_t)
    std_lambda = std(lambda_t)
    CV_X = std_lambda / mean_lambda
    return CV_X, mean_lambda, std_lambda
end

# Parameters
T = 20000.0  # Total simulation time in ms
dt = 0.1    # Time step in ms
lambda_h = 2.0  # Base firing rate
tau_lambda = 250.0  # Time constant
sigma_lambda = 3.0 * lambda_h /sqrt(200) # Noise intensity

# Simulate lambda(t)
lambda_t = simulate_lambda(T, dt, lambda_h, tau_lambda, sigma_lambda)

# Compute CV_X
CV_X, mean_lambda, std_lambda = compute_CV_X(lambda_t)

println("CV_X: $CV_X, Mean: $mean_lambda, Standard Deviation: $std_lambda")


