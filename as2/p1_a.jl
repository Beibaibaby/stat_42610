using NLsolve
using Plots
using LinearAlgebra  # Add this line to use eigvals
using Measures

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

# Plotting
p = scatter(stable_y, stable_v, label="Stable", color=:blue, markersize=4, markerstrokecolor=:blue,size=(800, 600),left_margin=20mm)
scatter!(p, unstable_y, unstable_v, label="Unstable", color=:red, markersize=4, markerstrokecolor=:red, linestyle=:dash)

xlabel!(p, "y")
ylabel!(p, "Fixed Point v*")
title!(p, "Stability of Fixed Points")
savefig("p1_a.png")
