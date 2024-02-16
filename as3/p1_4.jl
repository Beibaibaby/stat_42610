using Plots, Measures

# Define a range for aE values
aE_values = 0:0.1:5 # For example, from 0 to 1 µA/cm²

# Assuming a simple balance relation where aI = -aE for illustration
aI_values = -aE_values

# Plot the line
plot(aE_values, aI_values, label="aI vs. aE (Balance Line)", xlabel="aE (µA/cm²)", ylabel="aI (µA/cm²)", title="Balance Line in (aE, aI) Space",left_margin=15mm,bottom_margin=15mm)

# Save the plot
savefig("p1_4.png")