using QuadGK

# Define the parameters
lambda = 30/1000 # Example value, adjust as needed
tau_abs = 3.0 # Example value, adjust as needed
tau_rel = 3.0 # Example value, adjust as needed

# Define the corrected ρ(τ) function
function rho(tau, lambda, tau_rel, tau_abs)
    return lambda * (1 - exp(-(tau - tau_abs) / tau_rel)) * exp(-lambda * (tau - tau_abs + tau_rel * exp(-(tau - tau_abs) / tau_rel) - tau_rel))
end

# Numerical integration of ρ(τ) from tau_abs to infinity
function integrate_mean(lambda, tau_abs, tau_rel)
    integral, error = quadgk(tau -> tau*rho(tau, lambda, tau_rel, tau_abs), tau_abs, 10000)
    return integral
end

function integrate_second(lambda, tau_abs, tau_rel)
    integral, error = quadgk(tau -> tau^2*rho(tau, lambda, tau_rel, tau_abs), tau_abs, 10000)
    return integral
end



# Evaluate the integral
result_mean = integrate_mean(lambda, tau_abs, tau_rel)
println("Result of the integral (Mean): ", result_mean)

#second
result_second = integrate_second(lambda, tau_abs, tau_rel)
println("Result of the integral (Second): ", result_second)

#Compute the coeﬃcient of variance
cv = sqrt(result_second - result_mean^2) / result_mean
println("Coefficient of variance: ", cv)