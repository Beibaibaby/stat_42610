from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import sympy as sp
# Define the ODE system for numerical integration
def ode_system(vars, t, y_val):
    V, W = vars
    dv_dt_num = f_dv_dt(V, W, y_val)
    dw_dt_num = f_dw_dt(V, W, y_val)
    return [dv_dt_num, dw_dt_num]

# Define time span for the simulation
t_span = np.linspace(0, 100, 2000)  # Adjust as necessary for your system

# Simulate the system and find the min and max of V for oscillatory solutions
oscillatory_solutions = []

for y_val in y_range:
    # Use an initial guess near a possible limit cycle (if known)
    initial_conditions = [0, 0]  # Adjust based on the system's behavior
    sol = odeint(ode_system, initial_conditions, t_span, args=(y_val,))
    V_sol, W_sol = sol.T
    
    # Heuristic check for oscillations: variance in V over time exceeds a threshold
    if np.var(V_sol) > threshold:  # Define a suitable threshold for your system
        V_min, V_max = np.min(V_sol), np.max(V_sol)
        oscillatory_solutions.append((y_val, V_min, V_max))

# Unpack the oscillatory solutions for plotting
if oscillatory_solutions:
    osc_y, osc_min, osc_max = zip(*oscillatory_solutions)
else:
    osc_y, osc_min, osc_max = [], [], []

# Add to the existing plot
plt.figure(figsize=(10, 6))
plt.plot(stable_y, stable_v, 'bo', label='Stable Fixed Points', markersize=3)
plt.plot(unstable_y, unstable_v, 'ro', label='Unstable Fixed Points', markersize=3)
plt.plot(osc_y, osc_min, 'g.', label='Min Oscillatory Solutions', markersize=3)
plt.plot(osc_y, osc_max, 'm.', label='Max Oscillatory Solutions', markersize=3)
plt.xlabel('y')
plt.ylabel('v')
plt.title('Stable Oscillatory Solutions in the Fast Subsystem')
plt.legend()
plt.grid(True)
plt.savefig('oscillatory_solutions.png')