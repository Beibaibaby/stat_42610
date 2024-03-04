using Statistics, Plots, Measures

# Include the sim.jl file
include("sim.jl")
jie = 2.0 
jee = 1.0 

jei = 3.0
jii = 2.5 


stimstr_E = 1.2 
stimstr_I = 0.7 

# Define theoretical rates
rE_theory = ( stimstr_E*jii- stimstr_I *jei) / (jei*jie - jee*jii) * 1000/15
rI_theory = ( stimstr_E*jie- stimstr_I *jee)  / (jei*jie - jee*jii) * 1000/15
println("Theoretical E Neuron Rate: ", rE_theory)
println("Theoretical I Neuron Rate: ", rI_theory)
# Arrays to hold simulation results
ls = 1:7
Ns = 100 .* ls.^2
Ks = 0.2 .* Ns
ratesE = []
ratesI = []

for l in ls
    N = 100 * l^2
    println("Simulating for N = $N")
    times, ns, Ne, Ncells, T = sim(N)
    
    # Filter spikes that occur after 500ms
    filtered_spike_countsE = sum(times[1:Ne, :] .> 500)
    filtered_spike_countsI = sum(times[(Ne + 1):Ncells, :] .> 500)
    
    # Compute firing rates excluding the first 500ms
    # Note: Adjust the total time (T) to account for the exclusion of the first 500ms
    adjusted_T = T - 500
    rateE = 1000 * filtered_spike_countsE / (Ne * adjusted_T)
    rateI = 1000 * filtered_spike_countsI / ((Ncells - Ne) * adjusted_T)
    
    push!(ratesE, rateE)
    push!(ratesI, rateI)
end

println("E Neuron Rates: ", ratesE)
println("I Neuron Rates: ", ratesI)

using Plots
using Measures # Make sure this package is included for margin adjustments

# Assuming Ns, ratesE, ratesI, rE_theory, and rI_theory are already defined

# Plotting both E and I neuron firing rates on the same plot
p = plot(Ns, ratesE, label="E Neurons", color=:red, marker=:circle, 
         title="Firing Rate of Neurons", xlabel="N", ylabel="Firing Rate (Hz)", 
         left_margin=10mm, bottom_margin=10mm, ylim=(0, maximum([ratesE; ratesI; rE_theory; rI_theory]) * 1.1))

plot!(Ns, ratesI, label="I Neurons", color=:blue, marker=:square)

# Adding theoretical firing rate lines for both E and I neurons
hline!([rE_theory], label="rE_theory", color=:red, linestyle=:dash)
hline!([rI_theory], label="rI_theory", color=:blue, linestyle=:dash)

# Adjust layout if necessary, here it's a single plot so layout adjustment isn't needed
# plot!(layout=(1,1)) # This line is not necessary for a single plot

# Save the combined plot to a file
savefig("p2_2.png")



# Plotting
p3 = plot(Ns[5:7], ratesE[5:7], label = "E Neurons", color = :red, marker = :circle, title = "Firing Rate of E Neurons",left_margin = 10mm, bottom_margin = 10mm)
hline!([rE_theory], label = "rE_theory", color = :red, linestyle = :dash)
xlabel!("N")
ylabel!("Firing Rate (Hz)")

p4 = plot(Ns[5:7], ratesI[5:7], label = "I Neurons", color = :blue, marker = :square, title = "Firing Rate of I Neurons",left_margin = 10mm, bottom_margin = 10mm)
hline!([rI_theory], label = "rI_theory", color = :blue, linestyle = :dash)
xlabel!("N)")
ylabel!("Firing Rate (Hz)")

plot(p3, p4, layout = (2,1), size = (600, 400))
savefig("p2_2_zoomed.png")

