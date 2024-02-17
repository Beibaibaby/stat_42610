#this file is part of litwin-kumar_doiron_cluster_2012
#Copyright (C) 2014 Ashok Litwin-Kumar
#see README for more information

function sim()
	println("setting up parameters")
	Ncells = 2000
	Ne = 1000
	Ni = 1000
	T = 1000 #simulation time (ms)

	taue = 15 #membrane time constant for exc. neurons (ms)
	taui = 15

	#connection probabilities
	pei = 0.2
	pie = 0.2
	pii = 0.2
	pee = 0.2

	K = pei*Ne #average number of E->E connections per neuron
	
	sqrtK = sqrt(K)

	#jie = 4. / (taui*sqrtK)
	#jei = -16. * 1.2 /(taue*sqrtK)
	#jii = -16. / (taui*sqrtK)
	#jee = 1 / (sqrtk*taue)

	jie = 2.0 / (sqrtK)
	jei = -3.0 /(sqrtK)
	jii = -2.5 / (sqrtK)
	jee = 1.0 / (sqrtK)


	#stimulation
	Nstim = 400 #number of neurons to stimulate (indices 1 to Nstim will be stimulated)
	stimstr_E = 1.2/(sqrtK)
	stimstr_I = 0.7/(sqrtK)
	stimstart = 0
	stimend = T 

	#constant bias to each neuron type
	muemin = 1.1
	muemax = 1.2
	muimin = 1
	muimax = 1.05

	vre = 0. #reset voltage

	threshe = 1 #threshold for exc. neurons
	threshi = 1

	dt = .1 #simulation timestep (ms)
	refrac = 0.0 #refractory period (ms)

	maxrate = 150 #(Hz) maximum average firing rate.  if the average firing rate across the simulation for any neuron exceeds this value, some of that neuron's spikes will not be saved

	mu = zeros(Ncells)
	thresh = zeros(Ncells)
	tau = zeros(Ncells)

	mu[1:Ne] .= 0 #(muemax-muemin)*rand(Ne) .+ muemin
	mu[(Ne+1):Ncells] .= 0 #(muimax-muimin)*rand(Ni) .+ muimin

	thresh[1:Ne] .= threshe
	thresh[(1+Ne):Ncells] .= threshi

	tau[1:Ne] .= taue
	tau[(1+Ne):Ncells] .= taui


	weights = zeros(Ncells,Ncells)

	#random connections
	weights[1:Ne,1:Ne] = jee*(rand(Ne,Ne) .< pee)
	weights[1:Ne,(1+Ne):Ncells] = jei*(rand(Ne,Ni) .< pei)
	weights[(1+Ne):Ncells,1:Ne] = jie*(rand(Ni,Ne) .< pie)
	weights[(1+Ne):Ncells,(1+Ne):Ncells] = jii*(rand(Ni,Ni) .< pii)


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

	lastSpike = -100*ones(Ncells) #time of last spike

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
			#xerise[ci] += -dt*xerise[ci]/tauerise + forwardInputsEPrev[ci]
			#xedecay[ci] += -dt*xedecay[ci]/tauedecay + forwardInputsEPrev[ci]
			#xirise[ci] += -dt*xirise[ci]/tauirise + forwardInputsIPrev[ci]
			#xidecay[ci] += -dt*xidecay[ci]/tauidecay + forwardInputsIPrev[ci]

			synInput = forwardInputsEPrev[ci]+forwardInputsIPrev[ci]#(xedecay[ci] - xerise[ci])/(tauedecay - tauerise) + (xidecay[ci] - xirise[ci])/(tauidecay - tauirise)

			if (ci < Ne) && (t > stimstart) && (t < stimend) 
				synInput += stimstr_E;
			elseif (ci >= Ne) && (t > stimstart) && (t < stimend)
				synInput += stimstr_I;
			end

			if t > (lastSpike[ci] + refrac)  #not in refractory period
				v[ci] += dt*((1/tau[ci])*(mu[ci]-v[ci]) + synInput)

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
