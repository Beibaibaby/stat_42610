using Plots


# Parameters
beta = -1
I = 2
V_T = 1

# Functions for T
function T_case_1(alpha, Delta)
    return (1/sqrt(Delta)) * log(((beta + sqrt(Delta)) * V_T + 2 * I) / ((beta - sqrt(Delta)) * V_T + 2 * I))
end

function T_case_3(alpha, Delta)
    return (2/sqrt(-Delta)) * (atan((2 * alpha * V_T + beta) / sqrt(-Delta)) - atan(beta / sqrt(-Delta)))
end

# Generate the plot data
alpha_values = range(-1, stop=2, length=400)
T_values = [begin
                Delta = beta^2 - 4 * a * I
                Delta > 0 ? T_case_1(a, Delta) : T_case_3(a, Delta)
            end for a in alpha_values]

# Plotting save the plot
p=plot(alpha_values, T_values, label="T(α)", xlabel="α", ylabel="T", title="Curve of T")

#save the plot
savefig(p, "T_plot.png")