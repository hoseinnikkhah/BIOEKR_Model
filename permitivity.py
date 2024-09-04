import matplotlib.pyplot as plt
import numpy as np

# Define frequency and permittivity data (adjust as needed)
frequency = np.logspace(2, 10, 100)  # 100 points logarithmically spaced between 10^2 and 10^10
permittivity = [2.3, 2.28, 2.26, 2.24, 2.22, 2.2, 2.18, 2.16, 2.14, 2.12]  # Example data (replace with your actual values)

# Create the plot
plt.figure(figsize=(8, 6))
plt.plot(frequency, permittivity, color='blue')

# Set labels and title
plt.xlabel('Frequency (Hz)')
plt.ylabel('Relative permittivity')
plt.title('Typical Permittivity Spectrum of a Crude Oil')

# Set x-axis scale to logarithmic
plt.xscale('log')

# Gridlines
plt.grid(True)

# Show the plot
plt.show()