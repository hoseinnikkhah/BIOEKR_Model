import matplotlib.pyplot as plt
import numpy as np
import csv

gravity_API = np.array([20, 25, 30, 35, 40, 45, 50])
viscosity_at_100F = np.array([60, 20.62, 9.96, 4.5, 2.51, 1.316, 0.91])

gravity_API_smooth = np.linspace(np.min(gravity_API), np.max(gravity_API), 100)
smoothed_viscosity = np.interp(gravity_API_smooth, gravity_API, viscosity_at_100F)

csv_filename = 'oil_viscosity.csv'
with open(csv_filename, mode='w', newline='') as csvfile:
    csv_writer = csv.writer(csvfile)
    csv_writer.writerow(["API_gravity", "Viscosity_100F"])  # Header row
    csv_writer.writerows(zip(gravity_API_smooth, smoothed_viscosity))

plt.figure()

plt.plot(gravity_API_smooth, smoothed_viscosity, label='Viscosity at 100F', color='blue')
plt.yscale('log')
plt.xlabel('API Gravity of Oil')
plt.ylabel('Viscosity (Cp) at 100F')
plt.grid(True, which="both", linestyle="--", linewidth=0.7)
plt.legend()

plt.ylim(0.01, 100)

plt.show()
