using NLsolve
using DifferentialEquations
using Plots
using LinearAlgebra
using Measures
using ProgressMeter
using CSV
using DataFrames

# Define the system dynamics and Jacobian matrix
function system_dynamics!(F, x, y)
    v, w = x
    F[1] = y - 0.5 * (v + 0.5) - 2 * w * (v + 0.7) - (1/2 * (1 + tanh((v + 0.01) / 0.15))) * (v - 1)
    F[2] = 1.15 * ((1/2 * (1 + tanh((v - 0.1) / 0.145))) - w) * cosh((v - 0.1) / 0.29)
end

function jacobian_matrix(v, w)
    J11 = -2*w - (3.33333333333333 - 3.33333333333333*tanh(6.66666666666667*v + 0.066666666666667)^2)*(v - 1) - tanh(6.66666666666667*v + 0.066666666666667)/2 - 1.0
    J12 = -2*v - 1.4
    J21 = (3.96551724137931 - 3.96551724137931*tanh(6.89655172413793*v - 0.689655172413793)^2)*cosh(3.44827586206897*v - 0.344827586206897) + 3.44827586206897*(-1.15*w + 0.575*tanh(6.89655172413793*v - 0.689655172413793) + 0.575)*sinh(3.44827586206897*v - 0.344827586206897)
    J22 = -1.15*cosh(3.44827586206897*v - 0.344827586206897)
    return [J11 J12; J21 J22]
end

function is_stable(J)
    eigenvalues = eigvals(J)
    return all(real.(eigenvalues) .< 0)
end

# Prepare vectors to store y values and corresponding v values for stable and unstable points
stable_y = Float64[]
stable_v = Float64[]
stable_w = Float64[]
unstable_y = Float64[]
unstable_v = Float64[]
unstable_w = Float64[]

# Solve for fixed points and check their stability across a range of y values
#ys = -0.1:0.001:0.2
ys = -0.1:0.001:0.2

for y in ys
    # Range of initial guesses for v and w
    for v_guess in -2:0.5:2
        for w_guess in -2:0.5:2
            sol = nlsolve((F, x) -> system_dynamics!(F, x, y), [v_guess, w_guess])
            if converged(sol)
                v, w = sol.zero
                J = jacobian_matrix(v, w)
                if is_stable(J)
                    push!(stable_y, y)
                    push!(stable_v, v)
                    push!(stable_w, w)
                
                else
                    push!(unstable_y, y)
                    push!(unstable_v, v)
                    push!(unstable_w, w)

                end
            end
        end
    end
end



# Your existing functions system_dynamics!, jacobian_matrix, and is_stable go here

# Define the system of differential equations for the model
function bursting_model!(du, u, p, t)
    v, w = u
    y = p
    m_inf = 0.5 * (1 + tanh((v + 0.01) / 0.15))
    w_inf = 0.5 * (1 + tanh((v - 0.1) / 0.145))
    tau = cosh((v - 0.1) / 0.29)
    
    du[1] = y - 0.5*(v + 0.5) - 2*w*(v + 0.7) - m_inf*(v - 1) # dv/dt
    du[2] = 1.15*(w_inf - w)*tau # dw/dt
end


function count_crossings(values, threshold)
    crossings = 0
    for i in 2:length(values)
        if (values[i-1] < threshold && values[i] > threshold) || (values[i-1] > threshold && values[i] < threshold)
            crossings += 1
        end
    end
    return crossings
end

# Prepare for simulation
min_v_values = Float64[]
max_v_values = Float64[]
associated_y_min_v = Float64[]
associated_y_max_v = Float64[]


spk=zip(-0.1:0.01:0.2, -1:0.000001:1, -1:0.000001:1)

y_range = 0.05:0.001:0.12
v_range = -0.3:0.01:1  # Example with coarser granularity
w_range = -0.5:0.1:1  # Example with coarser granularity

total_iterations = length(y_range) * length(v_range) * length(w_range)
progress = Progress(total_iterations, desc="Processing", barlen=20)

#for (y, v, w) in spk# zip(stable_y, stable_v, stable_w)
for y in y_range
    for v in v_range
        for w in w_range

    next!(progress)
    initial_conditions = [v, w]
    tspan = (0.0, 3000.0)  # Adjust simulation time as necessary
    prob = ODEProblem(bursting_model!, initial_conditions, tspan, y)
    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

    # Find min and max of v for this simulation
    v_values = sol[1, 1000:end]
    if count_crossings(v_values, v) > 5
        min_v, max_v = minimum(v_values), maximum(v_values)
        push!(min_v_values, min_v)
        push!(max_v_values, max_v)
        push!(associated_y_min_v, y)
        push!(associated_y_max_v, y) 
    end
end
end
end





#=
total_iterations = length(y_range) * length(v_range) * length(w_range)
progress = Progress(total_iterations, desc="Processing", barlen=20)

# Repeat for unstable points
#for (y, v, w) in spk #zip(unstable_y, unstable_v, unstable_w)
for y in y_range
    for v in v_range
        for w in w_range
        next!(progress)
    
        initial_conditions = [v, w]
        tspan = (0.0, 3000.0)  # Adjust simulation time as necessary
        prob = ODEProblem(bursting_model!, initial_conditions, tspan, y)
        sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

        v_values = sol[1, 1000:end]
        if count_crossings(v_values, v) > 5
            min_v, max_v = minimum(v_values), maximum(v_values)
            push!(min_v_values, min_v)
            push!(max_v_values, max_v)
            push!(associated_y_min_v, y)
            push!(associated_y_max_v, y) 
        end
    
    
end
end
end
=#

# Initialize dictionaries to store the overall min and max v values for each unique y
min_v_for_y = Dict{Float64, Float64}()
max_v_for_y = Dict{Float64, Float64}()

# Process minimum v values
for (i, y) in enumerate(associated_y_min_v)
    min_v = min_v_values[i]
    # If y already exists in the dictionary, update it if the current min_v is smaller
    if haskey(min_v_for_y, y)
        min_v_for_y[y] = min(min_v_for_y[y], min_v)
    else
        min_v_for_y[y] = min_v
    end
end

# Process maximum v values
for (i, y) in enumerate(associated_y_max_v)
    max_v = max_v_values[i]
    # If y already exists in the dictionary, update it if the current max_v is larger
    if haskey(max_v_for_y, y)
        max_v_for_y[y] = max(max_v_for_y[y], max_v)
    else
        max_v_for_y[y] = max_v
    end
end

# Convert the dictionaries back to lists if needed, or use them as is depending on your requirements.
# For example, to get lists of unique y values and their associated min/max v values:
unique_y_values = collect(keys(min_v_for_y))
associated_min_v = [min_v_for_y[y] for y in unique_y_values]
associated_max_v = [max_v_for_y[y] for y in unique_y_values]


# Plotting
p = scatter(stable_y, stable_v, label="Stable", color=:blue, markersize=4, markerstrokecolor=:blue,size=(800, 600),left_margin=20mm)
scatter!(p, unstable_y, unstable_v, label="Unstable", color=:red, markersize=4, markerstrokecolor=:red, linestyle=:dash)

xlabel!(p, "y")
ylabel!(p, "Fixed Point v*")
title!(p, "Stability of Fixed Points")

# Plot min and max v values with their associated y values
scatter!(p, unique_y_values, min_v_values, label="Min v Values", color=:green, markersize=4, markerstrokecolor=:green)
scatter!(p, unique_y_values, max_v_values, label="Max v Values", color=:purple, markersize=4, markerstrokecolor=:purple)


savefig("p1_b.png")

# Create a DataFrame from the processed data
data = DataFrame(y=unique_y_values, min_v=min_v_values, max_v=max_v_values)

# Save the DataFrame to a CSV file
CSV.write("processed_v_values.csv", data)
