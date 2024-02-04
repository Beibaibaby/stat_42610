import sympy as sp

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

# Display the Jacobian matrix
sp.pprint(J, use_unicode=True)