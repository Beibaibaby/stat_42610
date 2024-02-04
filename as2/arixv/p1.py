import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import sympy as sp
from scipy.integrate import odeint
import pandas as pd

# Define symbols
v, w, y = sp.symbols('v w y', real=True)

# Define functions for m_inf, w_inf, and tau using the given expressions
m_inf = sp.Rational(1, 2) * (1 + sp.tanh((v + 0.01) / 0.15))
w_inf = sp.Rational(1, 2) * (1 + sp.tanh((v - 0.1) / 0.145))
tau = sp.cosh((v - 0.1) / 0.29)

# Define the system of equations
dv_dt = y - 0.5 * (v + 0.5) - 2 * w * (v + 0.7) - m_inf * (v - 1)
dw_dt = 1.15 * (w_inf - w) * tau

# Create sympy expressions for the equations
equations = sp.Matrix([dv_dt, dw_dt])

# Compute the Jacobian matrix using sympy
J = equations.jacobian([v, w])
print(J)
# Convert sympy expressions to lambda functions for numerical evaluation
f_dv_dt = sp.lambdify((v, w, y), dv_dt, 'numpy')
f_dw_dt = sp.lambdify((v, w, y), dw_dt, 'numpy')
f_J = sp.lambdify((v, w, y), J, 'numpy')

# Define a function for fsolve to find the fixed point
def find_fixed_point(y_val, guess_v, guess_w):
    # Find the root
    sol, infodict, ier, mesg = fsolve(lambda vars: [f_dv_dt(vars[0], vars[1], y_val),
                                                     f_dw_dt(vars[0], vars[1], y_val)], [guess_v, guess_w], full_output=True)
    # Check if the solution converged
    if ier == 1:
        v_star, w_star = sol
        # Compute the Jacobian at the fixed point
        J_val = f_J(v_star, w_star, y_val)
        # Calculate eigenvalues
        eigenvals = np.linalg.eigvals(J_val)
        # Determine stability: stable if all real parts of eigenvalues are negative
        is_stable = np.all(np.real(eigenvals) < 0)
        return v_star, is_stable
    else:
        return None, None

# Evaluate stability across the range of y values
y_range = np.linspace(-0.1, 0.2, 500)
v_stable = []
v_unstable = []

# Check each value of y
for y_val in y_range:
    # Try different initial guesses for v and w
    for guess_v in np.linspace(-2, 2, 5):
        for guess_w in np.linspace(-1, 1, 5):
            try:
                v_star, stability = find_fixed_point(y_val, guess_v, guess_w)
                if stability is not None:  # A fixed point was found
                    if stability:
                        v_stable.append((y_val, v_star))
                    else:
                        v_unstable.append((y_val, v_star))
            except Exception as e:
                # Handle exceptions raised by fsolve
                print(f"Exception at y={y_val}, guess_v={guess_v}, guess_w={guess_w}: {e}")

# Unpack the results for plotting
stable_y, stable_v = zip(*v_stable) if v_stable else ([], [])
unstable_y, unstable_v = zip(*v_unstable) if v_unstable else ([], [])

# Plot the fixed points v* as a function of y and save the figure
plt.figure(figsize=(10, 6))
plt.plot(stable_y, stable_v, 'bo', label='Stable', markersize=3)
plt.plot(unstable_y, unstable_v, 'ro', label='Unstable', markersize=3)
plt.xlabel('y')
plt.ylabel('Fixed Point v*')
plt.title('Fixed Points and their Stability')
plt.legend()
plt.grid(True)
plt.savefig('fixed_points_stability.png')





 # Assuming the CSV file is named 'bursting_model_results.csv'
# and contains columns 'Y', 'V_initial', and 'Behavior' where 'Behavior' indicates "Oscillatory" or "Non-Oscillatory"
df = pd.read_csv('bursting_model_results.csv')

# Filter rows where behavior is 'Oscillatory'
oscillatory_df = df[df['Behavior'] == 'Oscillatory']

# We will plot the range of V_initial for each Y value where behavior is 'Oscillatory'
oscillatory_y = oscillatory_df['Y'].unique()
oscillatory_min = []
oscillatory_max = []

for y in oscillatory_y:
    v_vals = oscillatory_df[oscillatory_df['Y'] == y]['V_initial']
    oscillatory_min.append((y, v_vals.min()))
    oscillatory_max.append((y, v_vals.max()))

# Unpack the min and max values for plotting
osc_y, osc_min = zip(*oscillatory_min) if oscillatory_min else ([], [])
_, osc_max = zip(*oscillatory_max) if oscillatory_max else ([], [])

# Now integrate this with your existing plotting code

plt.figure(figsize=(10, 6))
# Plot stable and unstable points (assuming you have these from your previous script)
plt.plot(stable_y, stable_v, 'bo', label='Stable Fixed Points', markersize=3)
plt.plot(unstable_y, unstable_v, 'ro', label='Unstable Fixed Points', markersize=3)
# Plot oscillatory min and max values
plt.plot(osc_y, osc_min, 'g^', label='Min Oscillatory Solutions', markersize=3)
plt.plot(osc_y, osc_max, 'yv', label='Max Oscillatory Solutions', markersize=3)

plt.xlabel('y')
plt.ylabel('v* / V_initial')
plt.title('Fixed Points, Stability, and Oscillatory Solutions')
plt.legend()
plt.grid(True)
plt.savefig('complete_analysis.png')

