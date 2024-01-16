
using DifferentialEquations, Plots

############################################################
# Constant variables in global scope:
# (Global scope gives ODE solver access without passing as parameters) 

const C = 1.0           # Membrance capacitance (uF/cm^2)
const gNabar = 120.0    # Max Na conductance (mS/cm^2)
const gKbar = 36.0      # Max K conductance (mS/cm^2)
const gLbar = 0.3       # Max leakage conductance (mS/cm^2)
const ENa = 45.0        # Na reversal potential (mV)
const EK = -82.0        # K reversal potential (mV)
const EL = -59.0        # Leakage reversal potential (mV)

############################################################

#################### Auxilliary Functions ##########################

function alpham(v)
    theta = (v + 45) / 10
    if theta == 0
        return 1.0
    else
        return 1.0 * theta / (1 - exp(-theta))
    end
end

function betam(v)
    return 4.0 * exp(-(v + 70) / 18)
end

function alphah(v)
    return 0.07 * exp(-(v + 70) / 20)
end

function betah(v)
    return 1.0 / (1 + exp(-(v + 40) / 10))
end

function alphan(v)
    theta = (v + 60) / 10
    if theta == 0
        return 0.1
    else
        return 0.1 * theta / (1 - exp(-theta))
    end
end

function betan(v)
    return 0.125 * exp(-(v + 70) / 80)
end

# This function will be applied to to ODE solver
function fxn!(du, u, p, t)
    du[1] = alpham(u[4]) * (1 - u[1]) - betam(u[4]) * u[1]
    du[2] = alphan(u[4]) * (1 - u[2]) - betan(u[4]) * u[2]
    du[3] = alphah(u[4]) * (1 - u[3]) - betah(u[4]) * u[3]
    gNa = gNabar * (u[1]^3) * u[3]
    gK = gKbar * (u[2]^4)

    iapp = p # Unpack parameters (for readability)
    du[4] = -(1 / C) * (gNa * (u[4] - ENa) + gK * (u[4] - EK) + gLbar * (u[4] - EL)) + iapp / C
end

##################################################################

###################### Main function ########################
#
# Michael Schwemmer, modified by Brent Doiron and Satchal Postlewaite
# Hodgkin-Huxley Model
# ====================
#
# simulates the response of the squid axon to a pulse of appiled current.
#
# parameters as in 
#   Hoppensteadt and Peskin
#   Modeling and Simulation in Medicine and the Life Sciences
#
# to run, call function from command line while in interactive session or from another script: 
#    
#    HH(vstart,iapp,T)

#    vstart is the initial voltage. 
#    iapp is the amplitude of the pulse. 
#    T simulation run time.
#
# voltages in mV, current in uA/cm^2, conductance in mS/cm^2, time is msec
#


function HH(vstart, iapp, T)
    # Initial Values for m, n, and h are set to rest state
    v = vstart
    m = alpham(v)/(alpham(v) + betam(v))
    n = alphan(v)/(alphan(v) + betan(v))
    h = alphah(v)/(alphah(v) + betah(v))

    # Set initial conditions for ODE solver
    u0 = [m, n, h, v]
    tspan = (0.0, T)
    params = iapp

    # Solve ODE
    dt = 0.01
    ode_prob = ODEProblem(fxn!, u0, tspan, params)
    sol = solve(ode_prob, alg_hints = [:stiff], saveat = dt, abstol = 1e-8, reltol = 1e-6)

    # sol.u returns as an array of arrays (# time steps, # state variables)
    # Unpack each state variable into array for plotting
    u = zeros(length(sol.u), length(sol.u[1]))
    for i = 1:length(sol.u)
        u[i,:] = sol.u[i]
    end

    # Plot and save figure
    plot_obj = plot(plot(sol.t, u[:,4], xlabel = "time", ylabel = "Voltage"),
                    plot(sol.t, u[:,1:3], label = ["m" "n" "h"], xlabel = "time"), layout=(3, 1), legend = false)
    filename = "HH.png"
    savefig(plot_obj, filename)

end

############################################################## 
#HH(-70, 6, 10)
##############################################################
#write a motified function based on HH named HH_planar_reduced by fixxng the value of m to be alpham(v)/(alpham(v) + betam(v)) and n to be 0.8-n
function fxn_r!(du, u, p, t)
    u[1] = alpham(u[4]) / (alpham(u[4]) + betam(u[4]))
    du[2] = alphan(u[4]) * (1 - u[2]) - betan(u[4]) * u[2]
    u[3] = 0.8 - u[2]
    gNa = gNabar * (u[1]^3) * u[3]
    gK = gKbar * (u[2]^4)

    iapp = p # Unpack parameters (for readability)
    du[4] = -(1 / C) * (gNa * (u[4] - ENa) + gK * (u[4] - EK) + gLbar * (u[4] - EL)) + iapp / C
end

function HH_planar_reduced(vstart, iapp, T)
   # Initial Values for m, n, and h are set to rest state
   v = vstart
   m = alpham(v)/(alpham(v) + betam(v))
   n = alphan(v)/(alphan(v) + betan(v))
   h = alphah(v)/(alphah(v) + betah(v))

   # Set initial conditions for ODE solver
   u0 = [m, n, h, v]
   tspan = (0.0, T)
   params = iapp

   # Solve ODE
   dt = 0.01
   ode_prob = ODEProblem(fxn_r!, u0, tspan, params)
   sol = solve(ode_prob, alg_hints = [:stiff], saveat = dt, abstol = 1e-8, reltol = 1e-6)

   # sol.u returns as an array of arrays (# time steps, # state variables)
   # Unpack each state variable into array for plotting
   u = zeros(length(sol.u), length(sol.u[1]))
   for i = 1:length(sol.u)
       u[i,:] = sol.u[i]
   end

   # Plot and save figure
   plot_obj = plot(plot(sol.t, u[:,4], xlabel = "time", ylabel = "Voltage"),
                   plot(sol.t, u[:,1:3], label = ["m" "n" "h"], xlabel = "time"), layout=(3, 1), legend = false)
   filename = "HH_reduced.png"
   savefig(plot_obj, filename)
end

HH_planar_reduced(-70, 6, 10)
HH(-70, 6, 10)