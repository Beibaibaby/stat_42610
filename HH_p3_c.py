import sympy as sp

# Define the symbols
V, n = sp.symbols('V n', real=True)

# Constants
C = 2.0
gNa = 20.0
gK = 20.0
gL = 2.0
ENa = 50.0
EK = -100.0
EL = -70.0
beta_m = -1.2
beta_n = 0  # Replace depending on type of class
gamma_m = 18.0
gamma_n = 10.0
phi = 0.15
I = 36.70  # Replace with actual current value if needed

# Define the functions
def x_inf(V, beta_x, gamma_x):
    return 0.5 * (1 + sp.tanh((V - beta_x) / gamma_x))

def tau_n(V, beta_x, gamma_x, phi):
    return 1.0 / (phi * sp.cosh((V - beta_x) / gamma_x))

# Compute derivatives
dx_inf_dV = sp.diff(x_inf(V, beta_m, gamma_m), V)
dtau_n_dV = sp.diff(tau_n(V, beta_n, gamma_n, phi), V)

# Define m_inf and n_inf for use in the Jacobian
m_inf = x_inf(V, beta_m, gamma_m)
n_inf = x_inf(V, beta_n, gamma_n)

# Define F and G functions from the neuron model
F = (1 / C) * (-gNa * m_inf * (V - ENa) - gK * n * (V - EK) - gL * (V - EL) + I)
G = (n_inf - n) / tau_n(V, beta_n, gamma_n, phi)

# Partial derivatives for the Jacobian
dF_dV = sp.diff(F, V)
dF_dn = sp.diff(F, n)
dG_dV = sp.diff(G, V)
dG_dn = sp.diff(G, n)

# Jacobian matrix
J = sp.Matrix([[dF_dV, dF_dn], [dG_dV, dG_dn]])

# Print the Jacobian matrix
print("Jacobian Matrix J:")
sp.pprint(J)

V_val = -42.00440593572987   # Example value for V
n_val = 0.0002246187963222547 # Example value for n
J_substituted = J.subs({V: V_val, n: n_val})
print("\nJacobian Matrix with V =", V_val, "and n =", n_val, ":")
sp.pprint(J_substituted)


# Compute and print the eigenvalues
eigenvalues = J_substituted.eigenvals()
print("\nEigenvalues of the Jacobian Matrix(Class I):")
sp.pprint(eigenvalues)


#######################################



V, n = sp.symbols('V n', real=True)

# Constants
C = 2.0
gNa = 20.0
gK = 20.0
gL = 2.0
ENa = 50.0
EK = -100.0
EL = -70.0
beta_m = -1.2
beta_n = -13  # Replace depending on type of class
gamma_m = 18.0
gamma_n = 10.0
phi = 0.15
I = 45.75  # Replace with actual current value if needed

# Define the functions
def x_inf(V, beta_x, gamma_x):
    return 0.5 * (1 + sp.tanh((V - beta_x) / gamma_x))

def tau_n(V, beta_x, gamma_x, phi):
    return 1.0 / (phi * sp.cosh((V - beta_x) / gamma_x))

# Compute derivatives
dx_inf_dV = sp.diff(x_inf(V, beta_m, gamma_m), V)
dtau_n_dV = sp.diff(tau_n(V, beta_n, gamma_n, phi), V)

# Define m_inf and n_inf for use in the Jacobian
m_inf = x_inf(V, beta_m, gamma_m)
n_inf = x_inf(V, beta_n, gamma_n)

# Define F and G functions from the neuron model
F = (1 / C) * (-gNa * m_inf * (V - ENa) - gK * n * (V - EK) - gL * (V - EL) + I)
G = (n_inf - n) / tau_n(V, beta_n, gamma_n, phi)

# Partial derivatives for the Jacobian
dF_dV = sp.diff(F, V)
dF_dn = sp.diff(F, n)
dG_dV = sp.diff(G, V)
dG_dn = sp.diff(G, n)

# Jacobian matrix
J = sp.Matrix([[dF_dV, dF_dn], [dG_dV, dG_dn]])

# Print the Jacobian matrix
print("Jacobian Matrix J:")
sp.pprint(J)

V_val = -35.795555450866665   # Example value for V
n_val = 0.010362849744594825 # Example value for n
J_substituted = J.subs({V: V_val, n: n_val})
print("\nJacobian Matrix with V =", V_val, "and n =", n_val, ":")
sp.pprint(J_substituted)


# Compute and print the eigenvalues
eigenvalues = J_substituted.eigenvals()
print("\nEigenvalues of the Jacobian Matrix (Class II):")
sp.pprint(eigenvalues)