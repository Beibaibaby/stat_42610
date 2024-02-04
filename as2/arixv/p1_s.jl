using DifferentialEquations, Plots, Measures

# Define the system of differential equations
function bursting_model!(du, u, p, t)
    v, w = u
    y = p
    m_inf = 0.5 * (1 + tanh((v + 0.01) / 0.15))
    w_inf = 0.5 * (1 + tanh((v - 0.1) / 0.145))
    tau = cosh((v - 0.1) / 0.29)
    
    du[1] = y - 0.5*(v + 0.5) - 2*w*(v + 0.7) - m_inf*(v - 1) # dv/dt
    du[2] = 1.15*(w_inf - w)*tau # dw/dt
end

# Set the parameters for the simulation
y = 0.069  # Example value of y
u0 = [0.034528054780827946, 0.28842180302784154]  # Initial conditions for v and w
tspan = (0.0, 1000.0)  # Time span for the simulation

# Create an ODE problem and solve it
prob = ODEProblem(bursting_model!, u0, tspan, y)
sol = solve(prob, Tsit5(), saveat=0.1)

# Plot v over time
plot(sol.t, getindex.(sol.u, 1), label="v(t)", xlabel="Time", ylabel="v", title="v over Time for y=$y",left_margin=10mm)
savefig("v_time_plot.png")  # Save the plot to a file
