import matplotlib.pyplot as plt
import numpy as np

# Original data
API_gravity = np.array([20, 25, 30, 35, 40, 45, 50])
viscosity_100F = np.array([60, 20.62, 9.96, 4.5, 2.51, 1.316, 0.91])

# Generate more API_gravity points using linear spacing
API_gravity_smooth = np.linspace(API_gravity.min(), API_gravity.max(), 100)

# Interpolate viscosity values for the smoother API_gravity points
viscosity_100F_smooth = np.interp(API_gravity_smooth, API_gravity, viscosity_100F)

# Plot
plt.figure()

plt.plot(API_gravity_smooth, viscosity_100F_smooth, marker='', label='T=100F', color='black')  # No marker for smoother appearance

# Logarithmic scale for y-axis
plt.yscale('log')

# Labels and title
plt.xlabel('Oil API gravity')
plt.ylabel('Dead oil viscosity, Cp')

plt.grid(True, which="both", ls="--", linewidth=0.5)
plt.legend()

# Set y-axis limits
plt.ylim([0.01, 100])

# Show plot
plt.show()
