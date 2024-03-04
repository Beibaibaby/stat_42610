using Plots
using Measures

function sim(N,weights_same,pertub_spike,pertub_v)
	println("setting up parameters")
    weights =copy(weights_same)
	
	Ne = N
	Ni = N
	Ncells= Ne + Ni

	T = 1000 #simulation time (ms)

	taue = 15 #membrane time constant for exc. neurons (ms)
	taui = 15
    
	p = 0.2
	#connection probabilities

	K = p * N #average number in-pop connections per neuron
	sqrtK = sqrt(K)

	jie = 2.0 /sqrtK
	jee = 1.0 /sqrtK

	jei = -3.0 /sqrtK
	jii = -2.5 /sqrtK


	stimstr_E = 1.2 * sqrtK 
	stimstr_I = 0.7 * sqrtK 

	stimstart = 0
	stimend = T 

	tauerise = 1
	tauedecay = 3
	tauirise = 1
	tauidecay = 2


	vre = 0.0 #reset voltage

	threshe = 1 #threshold for exc. neurons
	threshi = 1

	dt = 0.1 #simulation timestep (ms)
	refrac = 0.0 #refractory period (ms)

	maxrate = 1000 #(Hz) maximum average firing rate.  if the average firing rate across the simulation for any neuron exceeds this value, some of that neuron's spikes will not be saved

	mu = zeros(Ncells)
	thresh = zeros(Ncells)
	tau = zeros(Ncells)

	thresh[1:Ne] .= threshe
	thresh[(1+Ne):Ncells] .= threshi

	tau[1:Ne] .= taue
	tau[(1+Ne):Ncells] .= taui





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

	v = zeros(Ncells) #membrane voltage 

	lastSpike = -1000*ones(Ncells) #time of last spike

	Nsteps = round(Int,T/dt)
		
	println("starting simulation")

	#begin main simulation loop
	for ti = 1:Nsteps
		if mod(ti,Nsteps/100) == 1  #print percent complete
			print("\r",round(Int,100*ti/Nsteps))
		end
		t = dt*ti
		forwardInputsE[:] .= 0
		forwardInputsI[:] .= 0
		
		for ci = 1:Ncells
		   
			xerise[ci] += -dt*xerise[ci]/tauerise + forwardInputsEPrev[ci]
			xedecay[ci] += -dt*xedecay[ci]/tauedecay + forwardInputsEPrev[ci]
			xirise[ci] += -dt*xirise[ci]/tauirise + forwardInputsIPrev[ci]
			xidecay[ci] += -dt*xidecay[ci]/tauidecay + forwardInputsIPrev[ci]

			recurrent_input = (xedecay[ci] - xerise[ci])/(tauedecay - tauerise) + (xidecay[ci] - xirise[ci])/(tauidecay - tauirise)

            
			external_input = 0

			if (ci < Ne) && (t > stimstart) && (t < stimend) 
				external_input = stimstr_E;	
			elseif (ci >= Ne) && (t > stimstart) && (t < stimend)
				external_input = stimstr_I;
			end
            
		

			if t > (lastSpike[ci] + refrac)  #not in refractory period
				v[ci] += dt*((1/tau[ci])*(-v[ci]) + recurrent_input + external_input/tau[ci])
                
                if pertub_spike == true
                    if ci == 1
                        if ti == round(Int, 500 / dt)
                          v[ci] = 1.1
                        end
                    end
                end

                if pertub_v == true
                    if ci == 1
                        if ti == round(Int, 500 / dt)
                          v[ci] = v[ci] + 0.5
                        end
                    end
                end
    

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

	return times,ns,Ne,Ncells,T
end

N = 1000
Ncells = 2N
Ne=N
Ni=N

p = 0.2 #connection probabilities

K = p * N #average number in-pop connections per neuron
sqrtK = sqrt(K)

jie = 2.0 /sqrtK
jee = 1.0 /sqrtK

jei = -3.0 /sqrtK
jii = -2.5 /sqrtK


stimstr_E = 1.2 * sqrtK 
stimstr_I = 0.7 * sqrtK 


weights_same = zeros(Ncells,Ncells)

#random connections
weights_same[1:Ne,1:Ne] = jee*(rand(Ne,Ne) .< p)
weights_same[1:Ne,(1+Ne):Ncells] = jei*(rand(Ne,Ni) .< p)
weights_same[(1+Ne):Ncells,1:Ne] = jie*(rand(Ni,Ne) .< p)
weights_same[(1+Ne):Ncells,(1+Ne):Ncells] = jii*(rand(Ni,Ni) .< p)

for ci = 1:Ncells
    weights_same[ci,ci] = 0
end


times1, ns1, Ne, Ncells, T = sim(1000,weights_same,false,false)
times2, ns2, Ne, Ncells, T = sim(1000,weights_same,true,false)



using Plots

function plot_rasters(times1, ns1, times2, ns2, Ne, Ncells, T)
    p = plot(layout = (1, 1), size = (800, 600), legend = false,left_margin=25mm, bottom_margin=25mm)

    # Plot spikes from the first simulation in green
    for ci = 1:Ncells
        vals1 = times1[ci, 1:ns1[ci]]
        scatter!(vals1, fill(ci, length(vals1)), markersize=0.6, markerstrokewidth=0, color=:green)
    end
    
    # Plot spikes from the second simulation in purple
    for ci = 1:Ncells
        vals2 = times2[ci, 1:ns2[ci]]
        scatter!(vals2, fill(ci, length(vals2)), markersize=0.7, markerstrokewidth=0, color=:black)
    end

    xlims!(0, T)
    ylims!(0, Ncells)
    ylabel!("Neuron ID")
    xlabel!("Time (ms)")
    title!("Neuronal Firing Raster Plot")

    savefig(p, "p2_4.png")
end



# Ensure you have the results from both simulations: times1, ns1 for the first and times2, ns2 for the second
plot_rasters(times1, ns1, times2, ns2, Ne, Ncells, T)

function compute_firing_rates(times, ns, Ne, Ncells, T, window1 = (50, 500), window2 = (550, 1000))
    # Initialize counters for spikes within the specified windows
    spikes_window1_E = 0
    spikes_window2_E = 0
    spikes_window1_I = 0
    spikes_window2_I = 0

    for ci = 1:Ncells
        for spk = 1:ns[ci]
            spike_time = times[ci, spk]
            if spike_time >= window1[1] && spike_time <= window1[2]
                if ci <= Ne
                    spikes_window1_E += 1
                else
                    spikes_window1_I += 1
                end
            elseif spike_time >= window2[1] && spike_time <= window2[2]
                if ci <= Ne
                    spikes_window2_E += 1
                else
                    spikes_window2_I += 1
                end
            end
        end
    end

    # Calculate firing rates (spikes/second)
    duration1 = (window1[2] - window1[1]) / 1000.0  # Convert ms to seconds
    duration2 = (window2[2] - window2[1]) / 1000.0
    firing_rate_window1_E = spikes_window1_E / (Ne * duration1)
    firing_rate_window2_E = spikes_window2_E / (Ne * duration2)
    firing_rate_window1_I = spikes_window1_I / ((Ncells - Ne) * duration1)
    firing_rate_window2_I = spikes_window2_I / ((Ncells - Ne) * duration2)

    return firing_rate_window1_E, firing_rate_window2_E, firing_rate_window1_I, firing_rate_window2_I
end


# For the first simulation
firing_rate_window1_E_1, firing_rate_window2_E_1, firing_rate_window1_I_1, firing_rate_window2_I_1 = compute_firing_rates(times1, ns1, Ne, Ncells, T)
println("Firing rates for the first simulation:")
println("Window 1 (50-500 ms): Excitatory = ", firing_rate_window1_E_1, " Hz, Inhibitory = ", firing_rate_window1_I_1, " Hz")
println("Window 2 (550-1000 ms): Excitatory = ", firing_rate_window2_E_1, " Hz, Inhibitory = ", firing_rate_window2_I_1, " Hz")



# For the second simulation
firing_rate_window1_E_2, firing_rate_window2_E_2, firing_rate_window1_I_2, firing_rate_window2_I_2 = compute_firing_rates(times2, ns2, Ne,Ncells,T)
println("Firing rates for the second simulation:")
println("Window 1 (50-500 ms): Excitatory = ", firing_rate_window1_E_2, " Hz, Inhibitory = ", firing_rate_window1_I_2, " Hz")
println("Window 2 (550-1000 ms): Excitatory = ", firing_rate_window2_E_2, " Hz, Inhibitory = ", firing_rate_window2_I_2, " Hz")




