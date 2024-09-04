import numpy as np
import matplotlib.pyplot as plt

# Frequency range (logarithmic scale)
frequency = np.logspace(2, 11, 500)  # from 10^2 to 10^11 Hz

# Relative permittivity (arbitrary example curve)
permittivity = 2.3 - 0.1 / (1 + (frequency / 1e7)**2)  # Example function

# Create the plot
plt.figure(figsize=(8, 6))
plt.plot(frequency, permittivity, 'b')  # 'b' stands for blue line

# Set the scale to logarithmic
plt.xscale('log')

# Set axis labels
plt.xlabel('Frequency (Hz)')
plt.ylabel('Relative permittivity')

# Set axis limits similar to your image
plt.ylim([2.1, 2.35])
plt.xlim([10**2, 10**11])

# Add the title or figure caption
plt.title('Typical permittivity spectrum of a crude oil')

# Add grid
plt.grid(True, which="both", ls="--")

# Show the plot
plt.show()
