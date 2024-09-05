import matplotlib.pyplot as plt
import numpy as np

# Data for the plot
API_gravity = np.array([20, 25, 30, 35, 40, 45, 50])
viscosity_100F = np.array([60, 20.62, 9.96, 5.5, 2.51, 1.316, 0.91])
viscosity_150F = np.array([31.62, 10, 3.16, 1, 0.316, 0.1, 0.0316])
viscosity_200F = np.array([10, 3.16, 1, 0.316, 0.1, 0.0316, 0.01])
viscosity_250F = np.array([3.16, 1, 0.316, 0.1, 0.0316, 0.01, 0.00316])
viscosity_300F = np.array([1, 0.316, 0.1, 0.0316, 0.01, 0.00316, 0.001])

# Create the plot
plt.figure()

# Plot each temperature line
plt.plot(API_gravity, viscosity_100F, marker='s', label='T=100F', color='black')
plt.plot(API_gravity, viscosity_150F, marker='v', label='T=150F', color='black')
plt.plot(API_gravity, viscosity_200F, marker='^', label='T=200F', color='black')
plt.plot(API_gravity, viscosity_250F, marker='o', label='T=250F', color='black')
plt.plot(API_gravity, viscosity_300F, marker='s', label='T=300F', color='white', markerfacecolor='white')

# Logarithmic scale for y-axis
plt.yscale('log')

# Labels and title
plt.xlabel('Oil API gravity')
plt.ylabel('Dead oil viscosity, Cp')

# Grid
plt.grid(True, which="both", ls="--", linewidth=0.5)

# Add legend
plt.legend()

# Set y-axis limits
plt.ylim([0.01, 100])

# Show plot
plt.show()
