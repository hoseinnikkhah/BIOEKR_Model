import numpy as np
import matplotlib.pyplot as plt

# Sample data points
P = np.array([0, 0.01, 0.02, 0.07, 0.14, 0.23, 0.26])
B = np.array([0, 0.01, 0.04, 0.09, 0.17, 0.25, 0.30])
R = np.array([0, 0.03, 0.05, 0.07, 0.15, 0.23, 0.25])
Bl = np.array([0, 0.02, 0.04, 0.08, 0.14, 0.21, 0.24])

# Sample data points
epsilon = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])

# Fit a second-degree polynomial (adjust the degree as needed)
p = np.polyfit(epsilon, P, 2)
b = np.polyfit(epsilon, B, 2)
r = np.polyfit(epsilon, R, 2)
bl = np.polyfit(epsilon, Bl, 2)

# Generate x-values for missing points or extrapolation (replace with desired range)
new_epsilon = np.linspace(min(epsilon), max(epsilon), 100)

# Use the fitted polynomial to estimate epsilon for new x-values
new_P = np.polyval(p, new_epsilon)
new_B = np.polyval(b, new_epsilon)
new_R = np.polyval(r, new_epsilon)
new_Bl = np.polyval(bl, new_epsilon)

# Plot original data and fitted curves
plt.plot(epsilon, P, 'o', label='Data Points (Model 1)')
plt.plot(new_epsilon, new_P, 'm-', label='Model 1')
plt.plot(epsilon, B, 'o', label='Data Points (Model 2)')
plt.plot(new_epsilon, new_B, 'k-', label='Model 2')
plt.plot(epsilon, R, 'o', label='Data Points (Model 3)')
plt.plot(new_epsilon, new_R, 'r-', label='Model 3')
plt.plot(epsilon, Bl, 'o', label='Data Points (Model 4)')
plt.plot(new_epsilon, new_Bl, 'b-', label='Model 4')

# Set legends
plt.legend()

plt.show()
