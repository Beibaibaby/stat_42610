using Statistics
using Plots
using Measures  # For specifying margins in mm

doplot = true

include("sim.jl")  # Ensure this sim.jl script is correctly defined in your working directory

times, ns, Ne, Ncells, T = sim()

println("mean excitatory firing rate: ", mean(1000 * ns[1:Ne] / T), " Hz")
println("mean inhibitory firing rate: ", mean(1000 * ns[(Ne + 1):Ncells] / T), " Hz")

if doplot
    println("creating plot")
    p = plot(layout=(1,1), size=(600,400))  # Initialize an empty plot with specified size
    
    # Excitatory neurons
    for ci = 1:Ne
        vals = times[ci, 1:ns[ci]]
        scatter!(vals, fill(ci, length(vals)), label=false, markersize=1, markerstrokewidth=0, color=:red)
    end
    
    # Inhibitory neurons
    for ci = Ne+1:Ncells
        vals = times[ci, 1:ns[ci]]
        scatter!(vals, fill(ci, length(vals)), label=false, markersize=1, markerstrokewidth=0, color=:blue)
    end

    xlims!(0, T)
    ylims!(0, Ncells)
    ylabel!("Neuron")
    xlabel!("Time (ms)")
    title!("Neuronal Firing Raster Plot")

    plot!(margin=5mm)  # Set margins

    savefig("p2_1.png")
end
