using DifferentialEquations, Plots
using ProgressMeter

# Model parameters
C = 2.0
gL = 2.0
gK = 20.0
gNa = 20.0
EK = -100.0
ENa = 50.0
EL = -70.0
beta_m = -1.2
gamma_m = 18.0
gamma_n = 10.0
phi = 0.15


# Steady-state activation functions and time scale function
x_inf(V, beta_x, gamma_x) = 0.5 * (1.0 + tanh((V - beta_x) / gamma_x))
tau_n(V, beta_x, gamma_x) = 1.0 / (phi * cosh((V - beta_x) / gamma_x))

# Model equations
function neuron_model(du, u, p, t)
    I0, I1, f, beta_n = p
    I = I0 + I1 * sin(2 * pi * f * t)

    m_inf = x_inf(u[1], beta_m, gamma_m)
    n_inf = x_inf(u[1], beta_n, gamma_n)

    du[1] = (1 / C) * (-gNa * m_inf * (u[1] - ENa) - gK * u[2] * (u[1] - EK) - gL * (u[1] - EL) + I)
    du[2] = (n_inf - u[2]) / tau_n(u[1], beta_n, gamma_n)
end

# Simulation function
function simulate_neuron(beta_n, I0, I1, f, T, dt)
    V0 = -70.0
    n0 = x_inf(V0, beta_m, gamma_m)

    u0 = [V0, n0]
    tspan = (0.0, T)
    p = [I0, I1, f, beta_n]

    prob = ODEProblem(neuron_model, u0, tspan, p)
    #sol = solve(prob, Tsit5(), dt=dt) 
    sol = solve(prob, alg_hints = [:stiff], saveat = dt, abstol = 1e-8, reltol = 1e-6)

    return sol
end

# Class I excitability parameters
beta_n_class_I = 0.0
I0_class_I = 36.70
T = 400.0  # Total simulation time
dt = 0.005  # Time step
transient = 200.0  # Time to exclude for transient dynamics

# Parameter sweep range for I1 and f
I1_range = 0.4:0.02:2.0
f_range = 0.04:0.005:0.17

# Initialize result matrix
results_Class_I = zeros(length(I1_range), length(f_range))

# Perform simulations for Class I
total_iterations = length(I1_range) * length(f_range)

# Set up the progress bar
p = Progress(total_iterations, 1, "Computing... ", 50)

# Perform simulations for Class I
for (i, I1) in enumerate(I1_range)
    for (j, f) in enumerate(f_range)
        sol = simulate_neuron(beta_n_class_I, I0_class_I, I1, f, T, dt)
        # Exclude transient dynamics and check for spiking
        spiking = any([V > -40 for (t, V) in zip(sol.t, sol[1, :]) if t > transient])
        results_Class_I[i, j] = spiking ? 1 : 0

        # Update the progress bar
        next!(p)
    end
end

# Plotting the results
p = contour(f_range, I1_range, results_Class_I, title="Class I: Spiking Boundary", xlabel="Frequency f", ylabel="I1", levels=[0.5], colors=:blue)

#save plot
savefig(p, "Class_I.png")

heatmap_plot = heatmap(f_range, I1_range, results_Class_I,
                       title="Spiking Behavior",
                       xlabel="Frequency f",
                       ylabel="I1",
                       color=:blues,
                       aspect_ratio=:auto)

# Save the plot
savefig(heatmap_plot, "Class_I_heatmap.png")




#then complete it same as class I but change all the class I to class II
# Class II excitability parameters
beta_n_class_II = -13.0
I0_class_II = 45.75

p = Progress(total_iterations, 1, "Computing... ", 50)

# Perform simulations for Class II
for (i, I1) in enumerate(I1_range)
    for (j, f) in enumerate(f_range)
        sol = simulate_neuron(beta_n_class_II, I0_class_II, I1, f, T, dt)
        # Exclude transient dynamics and check for spiking
        spiking = any([V > -40 for (t, V) in zip(sol.t, sol[1, :]) if t > transient])
        results_Class_II[i, j] = spiking ? 1 : 0

        # Update the progress bar
        next!(p)
    end
end

# Plotting the results
p = contour(f_range, I1_range, results_Class_II, title="Class II: Spiking Boundary", xlabel="Frequency f", ylabel="I1", levels=[0.5], colors=:blue)

#save plot
savefig(p, "Class_II.png")

heatmap_plot = heatmap(f_range, I1_range, results_Class_II,
                       title="Spiking Behavior",
                       xlabel="Frequency f",
                       ylabel="I1",
                       color=:blues,
                       aspect_ratio=:auto)

# Save the plot
savefig(heatmap_plot, "Class_II_heatmap.png")




