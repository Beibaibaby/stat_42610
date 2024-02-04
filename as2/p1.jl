using DifferentialEquations, CSV, DataFrames, Plots, ProgressMeter

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

# Function to count the number of times v crosses a certain threshold
function count_crossings(v_sol, threshold)
    crossings = 0
    for i in 2:length(v_sol)
        if (v_sol[i-1] < threshold && v_sol[i] > threshold) || (v_sol[i-1] > threshold && v_sol[i] < threshold)
            crossings += 1
        end
    end
    return crossings
end

# Initialize arrays to hold the simulation results
y_values = -0.1:0.0005:0.2
v_initials = -2:0.5:2
w_initials = -2:0.5:2  # Span for w
results = []

@showprogress for y in y_values
    for v_init in v_initials
        for w_init in w_initials  # Iterate over w_initials
            # Initial conditions and parameters
            u0 = [v_init, w_init] # Initial [v, w]
            tspan = (0.0, 100.0)
            
            # Solve the differential equations
            prob = ODEProblem(bursting_model!, u0, tspan, y)
            sol = solve(prob, Tsit5(), saveat=0.1)
            
            # Check for oscillatory behavior by counting crossings
            v_sol = [u[1] for u in sol.u]
            threshold = 0 # Adjust this threshold as needed
            crossings = count_crossings(v_sol, threshold)
            
            # Determine if the solution is oscillatory based on the number of crossings
            if crossings > 20 # Adjust as needed based on your criteria for oscillatory behavior
                push!(results, (y, v_init, w_init, "Oscillatory"))
            else
                push!(results, (y, v_init, w_init, "Non-Oscillatory"))
            end
        end
    end
end

# Convert results to a DataFrame and save to a CSV file
df = DataFrame(results, [:Y, :V_initial, :W_initial, :Behavior])
CSV.write("bursting_model_results.csv", df)
