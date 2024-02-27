#this file is part of litwin-kumar_doiron_cluster_2012
#Copyright (C) 2014 Ashok Litwin-Kumar
#see README for more information

function sim(N)
	println("setting up parameters")
	
	Ne = N
	Ni = N
	Ncells= Ne + Ni

	T = 2500 #simulation time (ms)

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


