import numpy as np
import matplotlib.pyplot as plt
import csv
import os

# Simulated time data spanning 35 days
days = np.arange(1, 36)

# Given microbial concentration data (mg/L) over the time period
# This has no intercept (not correct as there is a medium for growth phase)
'''
microbial_concentration_up = np.array([
     1.02552020405620, 1.05169168892748, 1.07853107543312, 1.10605540855913, 1.13428216828303,
     1.16322928067492, 1.19291512928190, 1.22335856680290, 1.25457892706162, 1.28659603728484,
     1.31943023069425, 1.35310235941949, 1.38763380774081, 1.42304650566964, 1.45936294287580,
     1.49660618297005, 1.53479987815122, 1.57396828422707, 1.61413627601854, 1.65532936315706,
     1.69757370628505, 1.74089613366989, 1.78532415824180, 1.83088599506660, 1.87761057926434,
     1.92552758438526, 1.97466744125462, 2.02506135729857, 2.07674133636317, 2.12974019903911,
     2.18409160350528, 2.23983006690418, 2.29699098726279, 2.35561066597300, 2.41572633084560
])
'''
microbial_concentration = np.array([
       3.32011692, 3.25246043, 3.18618262, 3.1212554 , 3.05765126,
       2.99534322, 2.93430488, 2.87451037, 2.81593433, 2.75855194,
       2.70233887, 2.6472713 , 2.59332588, 2.54047975, 2.4887105 ,
       2.4379962 , 2.38831533, 2.33964685, 2.29197013, 2.24526494,
       2.19951151, 2.15469042, 2.11078269, 2.0677697 , 2.02563322,
       1.98435538, 1.9439187 , 1.90430602, 1.86550056, 1.82748586,
       1.79024582, 1.75376465, 1.71802688, 1.68301737, 1.64872127
])




# Calculate the natural logarithm of the concentration data
log_concentration = np.log(microbial_concentration)

# Plotting the log-transformed concentration data over time
plt.figure(figsize=(8, 5))
plt.plot(days, log_concentration, color='purple', marker='o', linestyle='--', label='ln(Concentration)')
plt.plot(days, microbial_concentration, color='green', marker='x', linestyle='-', label='Concentration')
plt.xlabel('Time (days)')
plt.ylabel('lnC')
plt.title('Logarithmic Growth of Microbial Concentration')
plt.grid(True, which='both', linestyle=':')
plt.legend()
plt.tight_layout()
plt.show()

# Perform linear regression using numpy's polyfit function
linear_fit = np.polyfit(days, log_concentration, 1)
growth_rate_k = linear_fit[0]
y_intercept = linear_fit[1]

# Validate if the file already exists
output_file = 'growth_results.csv'
if os.path.exists(output_file):
    print(f"Warning: '{output_file}' already exists and will be overwritten.")

# Save the growth rate and intercept to a CSV file
with open(output_file, mode='w', newline='') as file:
    csv_writer = csv.writer(file)
    csv_writer.writerow(['Growth Rate (k)', 'Y-Intercept'])
    csv_writer.writerow([growth_rate_k, y_intercept])

# Display calculated growth rate and intercept
print(f"Calculated Growth Rate (k): {growth_rate_k:.4f}")
print(f"Calculated Y-Intercept: {y_intercept:.4f}")
