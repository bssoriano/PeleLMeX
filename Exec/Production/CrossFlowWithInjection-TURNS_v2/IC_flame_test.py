import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# Constants
CH4_ID = 0  # Just for reference
O2_ID = 1
N2_ID = 2

# Define the mass fraction functions
#Z = (yf /y)^(1/2)*exp(−(V Pe/4) x^2/y) asymptotic solution for flameshape
def massfrac_CH4(x, y):
    return np.sqrt(1.8 / (y-0.1)) * np.exp(-(10.0 * 5.0 / 4.0) * (x * x) / (y-0.1))

def massfrac_O2(x, y):
    return 0.233

def massfrac_N2(x, y):
    return 1 - massfrac_CH4(x, y) - massfrac_O2(x, y)

# Create grid
x = np.linspace(0, 0.7, 100)
y = np.linspace(0, 1.4, 100)
X, Y = np.meshgrid(x, y)

# Calculate mass fractions
CH4 = massfrac_CH4(X, Y)
O2 = massfrac_O2(X, Y) * np.ones_like(X)  # Constant value
N2 = massfrac_N2(X, Y)

# Ensure mass fractions are physically reasonable (between 0 and 1)
#CH4 = np.clip(CH4, 0, 1)
#N2 = np.clip(N2, 0, 1)

# Plotting
fig, ax = plt.subplots(figsize=(10, 8))

# Plot CH4 mass fraction
contour = ax.contourf(X, Y, CH4, levels=20, cmap=cm.viridis)
cbar = fig.colorbar(contour, ax=ax)
cbar.set_label('CH4 Mass Fraction')

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('Mass Fraction of CH4')

plt.tight_layout()
plt.show()

# Optional: Plot all three species in subplots
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))

# CH4
cont1 = ax1.contourf(X, Y, CH4, levels=20, cmap=cm.viridis)
fig.colorbar(cont1, ax=ax1)
ax1.set_title('CH4 Mass Fraction')

# O2
cont2 = ax2.contourf(X, Y, O2, levels=20, cmap=cm.plasma)
fig.colorbar(cont2, ax=ax2)
ax2.set_title('O2 Mass Fraction')

# N2
cont3 = ax3.contourf(X, Y, N2, levels=20, cmap=cm.inferno)
fig.colorbar(cont3, ax=ax3)
ax3.set_title('N2 Mass Fraction')

plt.tight_layout()
plt.show()