import numpy as np
import matplotlib.pyplot as plt

# Define the frequency range for the plot, using a logarithmic scale
frequency = np.logspace(2, 11, 500)  # Frequencies ranging from 10^2 to 10^11 Hz

# Define the relative permittivity as a function of frequency
# This is an arbitrary example curve, meant to simulate the behavior of crude oil
permittivity = 2.15 + 0.15 / (1 + (frequency / 1e8)**1)  # Simple model for permittivity

# Initialize the plot with a specified figure size
plt.figure(figsize=(8, 6))

# Plot the frequency vs. permittivity, using a blue line ('b')
plt.plot(frequency, permittivity, 'b')

# Set the x-axis to a logarithmic scale to better represent the frequency range
plt.xscale('log')

# Label the axes to indicate what is being plotted
plt.xlabel('Frequency (Hz)')
plt.ylabel('Relative Permittivity')

# Define the limits of the axes, based on the expected range of the data
plt.ylim([2.1, 2.35])
plt.xlim([10**2, 10**11])

# Add a title that describes the content of the plot
plt.title('Typical Permittivity Spectrum of Crude Oil')

# Enable grid lines for better readability, with a dashed line style
plt.grid(True, which="both", linestyle="--")

# Display the plot to the user
plt.show()
