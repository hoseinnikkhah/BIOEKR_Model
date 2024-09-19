import numpy as np

# Define the slope
a = 0.0252

# Define the range of x values
x_values = np.arange(1, 37)

# Assuming b = 0 (no intercept)
b = 0

# Calculate y values
y_values = a * x_values + b

# Print the y values
print(y_values)
