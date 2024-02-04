using DifferentialEquations, Plots, Measures

# Define the system of differential equations, including dy/dt
function bursting_model!(du, u, p, t)
    v, w, y = u
    delta = 0.006368
    k = -0.22
    m_inf = 0.5 * (1 + tanh((v + 0.01) / 0.15))
    w_inf = 0.5 * (1 + tanh((v - 0.1) / 0.145))
    tau = cosh((v - 0.1) / 0.29)
    
    du[1] = y - 0.5*(v + 0.5) - 2*w*(v + 0.7) - m_inf*(v - 1)  # dv/dt
    du[2] = 1.15*(w_inf - w)*tau  # dw/dt
    du[3] = delta * (k - v)  # dy/dt
end

# Set initial conditions for v, w, and y
u0 = [0.0, 0.0, 0.0]  # Initial conditions

# Define the time span for the simulation
tspan = (0.0, 1000.0)

# Create an ODE problem and solve it
prob = ODEProblem(bursting_model!, u0, tspan)
sol = solve(prob, Tsit5(), saveat=0.1)

# Plot v(t) and y(t) from the solution
p1 = plot(sol.t, getindex.(sol.u, 1), label="v(t)", xlabel="Time", ylabel="v", title="v over Time", left_margin=10mm)
p2 = plot(sol.t, getindex.(sol.u, 3), label="y(t)", xlabel="Time", ylabel="y", title="y over Time", left_margin=10mm)

# Create a combined plot
plot(p1, p2, layout=(2, 1), size=(600, 400))

# Save the combined plot to a file
savefig("d_combined_v_y_time_plot.png")
