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

function count_spikes(voltage_trace, threshold, time_points, transient)
    num_spikes = 0
    is_above_threshold = false

    for (t, V) in zip(time_points, voltage_trace)
        if t > transient
            if V > threshold && !is_above_threshold
                num_spikes += 1
                is_above_threshold = true
            elseif V < threshold
                is_above_threshold = false
            end
        end
    end

    return num_spikes
end

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
    #sol = solve(prob, Euler(), dt=dt)
    sol = solve(prob, alg_hints = [:stiff], saveat = dt, abstol = 1e-8, reltol = 1e-6)
    #sol = solve(prob, Tsit5(), dt=dt) 
    return sol
end


# Modified Simulation Function
function simulate_and_plot_neuron(beta_n, I0, I1, f, T, dt, filename)
    sol = simulate_neuron(beta_n, I0, I1, f, T, dt)

    # Plotting
    p1 = plot(sol.t, sol[1,:], label="Voltage V(t)", xlabel="Time (ms)", ylabel="Voltage (mV)")
    p2 = plot(sol.t, sol[2,:], label="Gating variable n(t)", xlabel="Time (ms)", ylabel="n")
    
    plot(p1, p2, layout=(2,1), legend=true)
    
    # Save the plot
    savefig(filename)

    # Printing the steady-state values (the last values in the solution)
    V_ss = sol[1,end]
    n_ss = sol[2,end]
    println("Steady-state Voltage (V_ss): $V_ss mV")
    println("Steady-state Gating variable (n_ss): $n_ss")

    return sol
end

# Set parameters
beta_n_class_I = 0.0
I0_class_I = 36.70

I1 = 0.0
f = 0.0  # Frequency is zero since I1 is zero
T = 1000.0
dt = 0.005

beta_n_class_II = -13.0
I0_class_II = 45.75

# Run the simulation and plot
sol_I = simulate_and_plot_neuron(beta_n_class_I, I0_class_I, I1, f, T, dt,"SS_class_I.png")
sol_II = simulate_and_plot_neuron(beta_n_class_II, I0_class_II, I1, f, T, dt,"SS_class_II.png")




