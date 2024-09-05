import matplotlib.pyplot as plt
import numpy as np

# Data for the plot
API_gravity = np.array([20, 25, 30, 35, 40, 45, 50])
viscosity_100F = np.array([60, 20.62, 9.96, 4.5, 2.51, 1.316, 0.91])

plt.figure()

plt.plot(API_gravity, viscosity_100F, marker='s', label='T=100F', color='black')

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
