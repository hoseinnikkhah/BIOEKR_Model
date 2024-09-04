import numpy as np
import matplotlib.pyplot as plt
import csv 

frequency = np.logspace(2, 11, 500)  
permittivity = 2.15 + 0.15 / (1 + (frequency / 1e8)**1) 

# Save the data to a CSV file
with open('permittivity.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Frequency (Hz)', 'Relative Permittivity'])
    for f, p in zip(frequency, permittivity):
        writer.writerow([f, p])

# Plotting code
plt.figure(figsize=(8, 6))
plt.plot(frequency, permittivity, 'b')
plt.xscale('log')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Relative Permittivity')
plt.ylim([2.1, 2.35])
plt.xlim([10**2, 10**11])
plt.title('Typical Permittivity Spectrum of Crude Oil')
plt.grid(True, which="both", linestyle="--")
plt.show()
