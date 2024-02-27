using Plots
using Measures  # For margin adjustments

# Constants
tau1E = 4.0 # ms
tau2E = 0.4 # ms
tau1I = 6.0 # ms
tau2I = 1.75 # ms

# Define a range for aE values
aE_values = 0:0.1:2 # For example, from 0 to 2 µA/cm²

# Calculate corresponding aI values for the balance line
aI_values = -aE_values * (tau1E - tau2E) / (tau1I - tau2I)

# Plot the balance line
p = plot(aE_values, aI_values, label="aI vs. aE (Balance Line)", xlabel="aE (µA/cm²)", ylabel="aI (µA/cm²)", title="Balance Line in (aE, aI) Space", left_margin=15mm, bottom_margin=15mm)

# Points to highlight
highlight_aE = [0.99, 0.92, 0.86]
highlight_aI = -highlight_aE * (tau1E - tau2E) / (tau1I - tau2I)

# Add highlighted points to the plot
scatter!(highlight_aE, highlight_aI, label="Highlighted Points", color=:red, marker=:circle, markersize=6)

# Save the plot
savefig(p, "p1_4.png")
