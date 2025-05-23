import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# Parameters from the C++ code
nu = 3.99       # for Methane
YF_b = 0.1
YO_oo = 0.21
S_IC = nu * YF_b / YO_oo

q_dim = 50.15e6  # Heat of combustion (J/kg)
cp_tot = 1937.35  # Specific heat capacity (J/kg-K)
Too = 298.0       # Ambient temperature (K)
Q_IC = (q_dim * YF_b) / (cp_tot * Too)
H_IC = (S_IC + 1.0) * 1.0 / Q_IC + 1.0

# Domain parameters
y_or = 0.18       # Origin y-coordinate
x_or = 0.1       # Origin x-coordinate
y_sf = 0.4       # Scaling factor for y
V_asymp = 10.0    # Asymptotic velocity
Pe_asymp = 55.0   # Asymptotic Peclet number

# Create grid
x = np.linspace(0, 0.7, 100)
y = np.linspace(0.01, 1.4, 100)  # Avoid division by zero
X, Y = np.meshgrid(x, y)

# Initialize temperature and Z arrays
T = np.zeros_like(X)
Z = np.zeros_like(X)
Z_plot = np.zeros_like(X)


# Calculate temperature and Z fields
for i in range(len(x)):
    for j in range(len(y)):
        if Y[j, i] < y_or:
            T[j, i] = Too  # prob_parm.T_mean in original code
            Z[j, i] = 0.0  # No fuel below y_or
        else:
            # Calculate unnormalized Z
            term1 = np.sqrt(y_sf / (Y[j, i] - y_or))
            exponent = - (V_asymp * Pe_asymp / 4.0) * ((X[j, i] - x_or)**2) / (Y[j, i] - y_or)
            Z_unnormalized = term1 * np.exp(exponent)
            
            # Normalize Z
            Z_max_value = np.sqrt(y_sf / (Y[j, i] - y_or))
            Z[j, i] = Z_unnormalized / Z_max_value
            Z_plot[j,i] = Z_unnormalized
            # Calculate temperature
            if Z_unnormalized > 1:
                T[j, i] = 298.0 * (Q_IC / (S_IC + 1.0)) * (H_IC - (Z_unnormalized - 1.0) / S_IC)
            else:
                T[j, i] = 298.0 * (Q_IC / (S_IC + 1.0)) * (H_IC - (-Z_unnormalized + 1.0))

# Plotting with 1:1 aspect ratio
plt.figure(figsize=(8, 8))
contour = plt.contourf(X, Y, T, levels=50, cmap=cm.jet)
plt.colorbar(contour, label='Temperature (K)')

# Add Z=1 contour line
CS = plt.contour(X, Y, Z_plot, levels=[1.0], colors='white', linewidths=2, linestyles='dashed')
#plt.clabel(CS, inline=True, fmt='Z=1', fontsize=10)

plt.xlabel('x coordinate')
plt.ylabel('y coordinate')
plt.title('Temperature Distribution in Methane Combustion (Z=1 Contour)')
plt.gca().set_aspect('equal')  # 1:1 aspect ratio
plt.grid(True)
plt.show()