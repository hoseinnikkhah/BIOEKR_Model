import matplotlib.pyplot as plt
import numpy as np
import csv

gravity_API = np.array([20, 25, 30, 35, 40, 45, 50])
viscosity_at_100F = np.array([60, 20.62, 9.96, 4.5, 2.51, 1.316, 0.91])

# Fit a polynomial of degree 2 (you can adjust the degree if needed)
poly_coeff = np.polyfit(gravity_API, np.log(viscosity_at_100F), 2)  # log for log scale fitting
poly_fit = np.polyval(poly_coeff, gravity_API)

# Create a smooth API gravity range for plotting
gravity_API_smooth = np.linspace(np.min(gravity_API), np.max(gravity_API), 100)

# Calculate the viscosity using the polynomial fit
smoothed_viscosity = np.exp(np.polyval(poly_coeff, gravity_API_smooth))  # Exponentiate to undo log

# Write the data to CSV
csv_filename = 'oil_viscosity.csv'
with open(csv_filename, mode='w', newline='') as csvfile:
    csv_writer = csv.writer(csvfile)
    csv_writer.writerow(["API_gravity", "Viscosity_100F"])  # Header row
    csv_writer.writerows(zip(gravity_API_smooth, smoothed_viscosity))

# Plot the results
plt.figure()

plt.plot(gravity_API_smooth, smoothed_viscosity, label='Polyfit Viscosity', color='blue')
plt.scatter(gravity_API, viscosity_at_100F, label='Original Data', color='red', zorder=5)
plt.yscale('log')
plt.xlabel('API Gravity of Oil')
plt.ylabel('Viscosity (Cp)')
plt.grid(True, which="both", linestyle="--", linewidth=0.7)
plt.legend()

plt.ylim(0.01, 100)

plt.show()
